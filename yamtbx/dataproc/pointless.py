"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
import os
import re
from yamtbx import util
from libtbx.utils import null_out
from cctbx import sgtbx
from cctbx import crystal
from yamtbx.dataproc.xds import xds_ascii

pointless_comm = "pointless"
re_cell = re.compile(r"([0-9]+\.[0-9][0-9]) *([0-9]+\.[0-9][0-9]) *([0-9]+\.[0-9][0-9]) +([0-9]+\.[0-9][0-9]) *([0-9]+\.[0-9][0-9]) *([0-9]+\.[0-9][0-9])")

"""
Pointless's cell format
scala_util.cpp:
  std::string FormatCell(const Scell& cell, const int w, const int p)
  {
    clipper::String line;
    for (int i=0;i<3;i++) line += StringUtil::ftos(cell[i],w,p);
    line += "   ";
    for (int i=0;i<3;i++) line += StringUtil::ftos(cell[i+3],w,p);
    return line;
  }

w=7 and p=2 is default (defined in hkl_datatypes.hh; Scell::format())
"""

def parse_pointless_output_for_symm(output):
    lines = output.splitlines()
    ret = {}

    best_symm, best_cell = None, None
    read_flag = False
    for l in lines:
        if l.startswith("Best Solution:"):
            group_type = l.split()[2]
            assert group_type in ("point", "space") # How can we support Laue?
            best_symm = l[l.index("group")+6:].strip()
            read_flag = True
        elif read_flag and "Reindex operator:" in l:
            tmp = l.split(":")[1]
            tmp = tmp[tmp.index("[")+1:tmp.index("]")]
            ret["reindex_op"] = sgtbx.change_of_basis_op(tmp)
        elif read_flag and "Laue group probability:" in l:
            ret["laue_prob"] = util.safe_float(l.split(":")[1])
        elif read_flag and "Unit cell:" in l:
            r = re_cell.search(l)
            if r:
                best_cell = list(map(float, r.groups()))

    if None not in (best_symm, best_cell):
        ret["symm"] = crystal.symmetry(unit_cell=best_cell, space_group_symbol=best_symm,
                                       assert_is_compatible_unit_cell=False)
    return ret
# parse_pointless_output_for_symm()

def parse_pointless_output_for_input_files(output):
    lines = output.splitlines()
    ret = []

    for i, l in enumerate(lines):
        if "Reading XDS ascii file from file" in l:
            j = i + 1
            filename = ""
            while True:
                if "Header lines:" in lines[j]: break
                filename += lines[j].rstrip()
                j += 1
            ret.append(filename)
    return ret
# parse_pointless_output_for_runs()

def parse_pointless_output_for_runs(output):
    lines = output.splitlines()

    # Possibly line is broken..?
    re_run = re.compile(r"Run number: *([0-9]+) *consists of batches *([0-9]+) *to *([0-9]+)")
    ret = []

    for i, l in enumerate(lines):
        if "Run number:" in l:
            r = re_run.search(l)
            if r:
                run, sbatch, ebatch = list(map(int, r.groups()))
                ret.append((run, sbatch, ebatch))
    return ret
# parse_pointless_output_for_runs()


class Pointless:
    def __init__(self):#, hklin=None, xdsin=None, hklref=None, hklout=None, logout=None):
        pass
        #if hklref is None:
        #    self.run_for_symm(hklin, xdsin, logout)
    # __init__()

    def input_str(self, hklin=None, xdsin=None):
        filein = ""
        if hklin is not None:
            if isinstance(hklin, str):
                filein = "hklin %s" % os.path.basename(hklin)
                wdir = os.path.dirname(hklin)
            else:
                filein = "\n".join(["hklin %s"%x for x in hklin])
                wdir = os.path.dirname(hklin[0])
        else:
            if isinstance(xdsin, str):
                filein = "xdsin %s" % os.path.basename(xdsin)
                wdir = os.path.dirname(xdsin)
            else:
                filein = "\n".join(["xdsin %s"%x for x in xdsin])
                wdir = os.path.dirname(xdsin[0])
        return filein+"\n", wdir
    # input_str()

    def run_for_symm(self, hklin=None, xdsin=None, logout=None, tolerance=None, d_min=None, choose_laue=None, xdsin_to_p1=False, setting="c2"):
        assert (hklin, xdsin).count(None) == 1
        assert setting.lower().startswith(("c2", "cell", "symm"))
        tmpfile = None

        out = open(logout, "w") if logout else null_out()
        
        # Merge in P1 before pointless for quick evaluation
        if xdsin and xdsin_to_p1:
            out.write("Merging in P1 symmetry before Pointless..\n")
            xac = xds_ascii.XDS_ASCII(xdsin)
            xac.remove_rejected()
            i_obs = xac.i_obs(anomalous_flag=True)
            i_obs.show_summary(out, "  ")
            i_obs = i_obs.apply_change_of_basis(change_of_basis="to_niggli_cell", out=out)[0]
            i_obs = i_obs.customized_copy(space_group_info=sgtbx.space_group_info("P1"))
            i_obs = i_obs.merge_equivalents(use_internal_variance=False).array()
            tmpfile = util.get_temp_filename(prefix="", suffix=os.path.basename(xdsin) + "_p1.mtz")
            i_obs.as_mtz_dataset(column_root_label="I").mtz_object().write(tmpfile)
            xdsin, hklin = None, tmpfile
            out.write("\n")
        
        filein, wdir = self.input_str(hklin, xdsin)

        # Need centre keyword?
        stdin = """\
SETTING %s
%s
""" % (setting, filein)

        if tolerance is not None:
            stdin += "tolerance %f\n" % tolerance
        if d_min is not None:
            stdin += "resolution high %f\n" % d_min
        if choose_laue is not None:
            stdin += "choose lauegroup %s\n" % choose_laue

        exitcode, output, err = util.call(pointless_comm, 
                                          stdin=stdin, wdir=wdir)
        
        parsed = parse_pointless_output_for_symm(output)
        out.write(output)

        if tmpfile and os.path.isfile(tmpfile):
            os.remove(tmpfile)
        
        return parsed
    # run_for_symm()

    def run_copy(self, hklout, wdir, hklin=None, xdsin=None, logout=None, tolerance=None):
        assert (hklin, xdsin).count(None) == 1
        stdin, junk = self.input_str(hklin, xdsin)
        

        if tolerance is not None:
            stdin += "tolerance %f\n" % tolerance
        
        exitcode, output, err = util.call(pointless_comm, "-copy hklout %s" % hklout,
                                          stdin=stdin, wdir=wdir)
        if logout is not None:
            open(logout, "w").write(output)
            return parse_pointless_output_for_runs(output)
    # run_copy()

# class Pointless
