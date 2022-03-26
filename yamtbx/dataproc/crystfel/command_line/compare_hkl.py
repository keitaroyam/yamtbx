from __future__ import print_function
from __future__ import unicode_literals
from yamtbx.dataproc import crystfel
import iotbx.phil
from cctbx.array_family import flex
import math

master_params_str = """\
pdb = None
 .type = path
 .help = Symmetry source
datout = shells_all.dat
 .type = path
 .help = Output dat file

dmin = None
 .type = float
dmax = None
 .type = float
nshells = 20
 .type = int

fom = *cc *ccano *rsplit
 .type = choice(multi=True)
"""

def calc_cc(a1, a2):
    correlation = flex.linear_correlation(a1.data(), a2.data())
    if correlation.is_well_defined(): return correlation.coefficient()
    else: return float("nan")
# calc_cc()

def calc_ccano(a1, a2):
    if not (a1.anomalous_flag() and a2.anomalous_flag()):
        return float("nan")

    a1, a2 = [x.anomalous_differences() for x in (a1,a2)]
    return calc_cc(a1, a2)
# calc_ccano()

def calc_rsplit(a1, a2):
    num = flex.sum(flex.abs(a1.data() - a2.data()))
    den = flex.sum(a1.data() + a2.data()) / 2.

    return 1./math.sqrt(2.) * num / den
# clac_rsplit()


def run(hklfiles, params):
    arrays = [crystfel.hkl.HKLfile(symm_source=params.pdb, hklin=x) for x in hklfiles]
    for a in arrays: a.set_resolution(d_min=params.dmin, d_max=params.dmax)
    
    ofs = open(params.datout, "w")
    ofs.write("  dmax   dmin  nref  cmpl   red1   red2 %s\n" % " ".join(["%7s"%x for x in params.fom]))

    a1, a2 = arrays[0].array.common_sets(arrays[1].array)
    r1, r2 = arrays[0].redundancies.common_sets(arrays[1].redundancies)
    r1, r2 = [x.as_double() for x in (r1, r2)]

    binner = a1.setup_binner(n_bins=params.nshells)
    
    for i_bin in binner.range_used():
        sel = binner.selection(i_bin)
        d_max, d_min = binner.bin_d_range(i_bin)

        r1s, r2s = r1.select(sel), r2.select(sel)
        r1sm = flex.mean(r1s.data()) if r1s.size()>0 else float("nan")
        r2sm = flex.mean(r2s.data()) if r2s.size()>0 else float("nan")

        a1s, a2s = a1.select(sel), a2.select(sel)
        ofs.write("%6.2f %6.2f %5d %5.1f %6.1f %6.1f " % (d_max, d_min, a1s.size(),
                                                          a1s.completeness(d_max=d_max)*100.,
                                                          r1sm, r2sm))

        if "cc" in params.fom: ofs.write("% 7.4f " % calc_cc(a1s, a2s))
        if "ccano" in params.fom: ofs.write("% 7.4f " % calc_ccano(a1s, a2s))
        if "rsplit" in params.fom: ofs.write("% 7.4f " % calc_rsplit(a1s, a2s))
        ofs.write("\n")

    ofs.write("#     overall %5d %5.1f %6.1f %6.1f " % (a1.size(),
                                                        a1.completeness(d_max=params.dmax)*100.,
                                                        flex.mean(r1.data()), flex.mean(r2.data())))
    if "cc" in params.fom: ofs.write("% 7.4f " % calc_cc(a1, a2))
    if "ccano" in params.fom: ofs.write("% 7.4f " % calc_ccano(a1, a2))
    if "rsplit" in params.fom: ofs.write("% 7.4f " % calc_rsplit(a1, a2))
    ofs.write("\n")
# run()

if __name__ == "__main__":
    import sys
    import os

    if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
        print("All parameters:\n")
        iotbx.phil.parse(master_params_str).show(prefix="  ", attributes_level=1)
        quit()

    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args
    hklfiles = []

    for arg in args:
        if os.path.isfile(arg):
            if ".hkl" in arg: hklfiles.append(arg)
            elif params.pdb is None: params.pdb = arg

    if params.pdb is None:
        print("Give pdb file")
        quit()
    if len(hklfiles) != 2:
        print("Give two hkl files.")
        quit()

    run(hklfiles, params)
