#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2017. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals

from yamtbx.util import call
from yamtbx.util.xtal import format_unit_cell
import iotbx.file_reader
import iotbx.scalepack.merge
import iotbx.shelx.hklf
from cctbx import miller
from cctbx.array_family import flex
from iotbx import crystal_symmetry_from_any
from mmtbx.scaling.matthews import p_vm_calculator
from mmtbx.scaling.absolute_scaling import ml_aniso_absolute_scaling
import os

master_phil = """
hklin = None
  .type = path
pdbin = None
  .type = path
anisotropy_correction = False
  .type = bool
wdir = anode
  .type = str
"""

def check_symm(xs_hkl, xs_pdb):
    pdb_laue = xs_pdb.space_group().build_derived_reflection_intensity_group(False)
    hkl_laue = xs_hkl.space_group().build_derived_reflection_intensity_group(False)
    if pdb_laue != hkl_laue or not xs_pdb.unit_cell().is_similar_to(xs_hkl.unit_cell()):
        print("WARNING! Incompatible symmetry (symmetry mismatch or too different unit cell)")
        print(" hklin:")
        xs_hkl.show_summary(prefix="  ")
        print(" pdbin:")
        xs_pdb.show_summary(prefix="  ")
# check_symm()

def pha2mtz(phain, xs, mtzout):
    hkl, f, fom, phi, sigf = [], [], [], [], []
    
    for l in open(phain):
        sp = l.split()
        if len(sp) != 7: break
        hkl.append(tuple(map(int, sp[:3])))
        f.append(float(sp[3]))
        fom.append(float(sp[4]))
        phi.append(float(sp[5]))
        sigf.append(float(sp[6]))

    if not hkl:
        return
        
    f_array  = miller.array(miller.set(xs, flex.miller_index(hkl)),
                            data=flex.double(f),
                            sigmas=flex.double(sigf))
    mtz_ds = f_array.as_mtz_dataset(column_root_label="ANOM", column_types="FQ") # To open with Coot, column type F is required (not D)
    mtz_ds.add_miller_array(f_array.customized_copy(data=flex.double(phi),
                                                    sigmas=None),
                            column_root_label="PANOM",
                            column_types="P")
    mtz_ds.add_miller_array(f_array.customized_copy(data=flex.double(fom),
                                                    sigmas=None),
                            column_root_label="FOM",
                            column_types="W")

    mtz_ds.mtz_object().write(mtzout)
# pha2mtz()

def run(hklin, pdbin, wdir, anisotropy_correction=False):
    arrays = iotbx.file_reader.any_file(hklin).file_server.miller_arrays
    i_arrays = [x for x in arrays if x.is_xray_intensity_array() and x.anomalous_flag()]
    f_arrays = [x for x in arrays if x.is_xray_amplitude_array() and x.anomalous_flag()]

    if not i_arrays and not f_arrays:
        print("No anomalous observation data")
        return

    if os.path.exists(wdir):
        print("%s already exists. quiting." % wdir)
        return

    os.mkdir(wdir)

    xs = crystal_symmetry_from_any.extract_from(pdbin)

    sh_out = open(os.path.join(wdir, "run_anode.sh"), "w")
    sh_out.write("#!/bin/sh\n\n")
    sh_out.write("shelxc anode <<+ > shelxc.log 2>&1\n")
    sh_out.write("cell %s\n" % format_unit_cell(xs.unit_cell()))
    sh_out.write("spag %s\n" % str(xs.space_group_info()).replace(" ",""))

    if i_arrays:
        obs_array = i_arrays[0]
        infile = "%s.hkl" % os.path.splitext(os.path.basename(hklin))[0]
        in_opt = "%s" % infile
        print("Using intensity array:", obs_array.info().label_string())
    else:
        obs_array = f_arrays[0]
        infile = "%s_f.hkl" % os.path.splitext(os.path.basename(hklin))[0]
        in_opt ="-f %s" % infile
        print("No intensity arrays. Using amplitude arrays instead:", obs_array.info().label_string())
        
    sh_out.write("! data from %s : %s\n" % (os.path.abspath(hklin), obs_array.info().label_string()))
    obs_array.crystal_symmetry().show_summary(sh_out, prefix="! ")
    check_symm(obs_array.crystal_symmetry(), xs)

    n_org = obs_array.size()
    obs_array = obs_array.eliminate_sys_absent()
    n_sys_abs = n_org - obs_array.size()
    if n_sys_abs > 0:
        print("  %d systematic absences removed." % n_sys_abs)

    if anisotropy_correction:
        print("Correcting anisotropy..")
        n_residues = p_vm_calculator(obs_array, 1, 0).best_guess
        abss = ml_aniso_absolute_scaling(obs_array, n_residues=n_residues)
        abss.show()
        tmp = -2. if i_arrays else -1.
        b_cart = [x*tmp for x in abss.b_cart]
        obs_array = obs_array.apply_debye_waller_factors(b_cart=b_cart)
        
        
    sh_out.write("sad %s\n" % in_opt)
    iotbx.shelx.hklf.miller_array_export_as_shelx_hklf(obs_array, open(os.path.join(wdir, infile), "w"),
                                                       normalise_if_format_overflow=True)
    sh_out.write("+\n\n")
    sh_out.write('ln -s "%s" anode.pdb\n\n' % os.path.relpath(pdbin, wdir))
    sh_out.write("anode anode\n")
    sh_out.close()

    call(cmd="sh", arg="./run_anode.sh", wdir=wdir)

    pha_file = os.path.join(wdir, "anode.pha")
    if os.path.isfile(pha_file):
        pha2mtz(pha_file, xs, os.path.join(wdir, "anode.pha.mtz"))
    
    print("Done. See %s/" % wdir)

    fa_file = os.path.join(wdir, "anode_fa.hkl")
    if os.path.isfile(fa_file):
        r = iotbx.shelx.hklf.reader(open(fa_file))
        fa_array = r.as_miller_arrays(crystal_symmetry=xs)[0]
        print("\nData stats:")
        print(" # Cmpl.o = Anomalous completeness in original data")
        print(" # Cmpl.c = Anomalous completeness in shelxc result (rejections)")
        print(" # SigAno = <d''/sigma> in shelxc result")
        print(" d_max d_min Cmpl.o Cmpl.c SigAno")
        binner = obs_array.setup_binner(n_bins=12)
        for i_bin in binner.range_used():
            d_max_bin, d_min_bin = binner.bin_d_range(i_bin)
            obs_sel = obs_array.resolution_filter(d_max_bin, d_min_bin)
            obs_sel_ano = obs_sel.anomalous_differences()
            fa_sel = fa_array.resolution_filter(d_max_bin, d_min_bin)
            cmplset = obs_sel_ano.complete_set(d_max=d_max_bin, d_min=d_min_bin).select_acentric()
            n_acentric = cmplset.size()
            sigano = flex.mean(fa_sel.data()/fa_sel.sigmas()) if fa_sel.size() else float("nan")
            print(" %5.2f %5.2f %6.2f %6.2f %6.2f" % (d_max_bin, d_min_bin,
                                                      100.*obs_sel_ano.size()/n_acentric,
                                                      100.*fa_sel.size()/n_acentric,
                                                      sigano))

    lsa_file = os.path.join(wdir, "anode.lsa")
    if os.path.isfile(lsa_file):
        print("")
        flag = False
        for l in open(lsa_file):
            if "Strongest unique anomalous peaks" in l:
                flag  = True
            elif "Reflections written to" in l:
                flag = False
            if flag:
                print(l.rstrip())

    if os.path.isfile(("anode_fa.res")):
        x = iotbx.shelx.cctbx_xray_structure_from(file=open("anode_fa.res"))
        open("anode_fa.pdb", "w").write(x.as_pdb_file())

        
# run()

if __name__ == "__main__":
    import sys
    import iotbx.phil
    
    cmdline = iotbx.phil.process_command_line_with_files(args=sys.argv[1:],
                                                         master_phil_string=master_phil,
                                                         reflection_file_def="hklin",
                                                         pdb_file_def="pdbin")
    params = cmdline.work.extract()
    run(params.hklin, params.pdbin, params.wdir, params.anisotropy_correction)
