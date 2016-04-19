"""
To check if peak of DANO map is due to model bias or not.

Pass output mtz file of phenix.refine to this script.
And then run FFT.
Can you still see peak of anomalous scatterers?
"""

import sys
import iotbx.mtz
from cctbx.array_family import flex
import numpy

def get_cc(a1, a2, title):
    a1, a2 = a1.common_sets(a2)
    corr = flex.linear_correlation(a1.data(), a2.data())
    assert corr.is_well_defined()
    #print "%30s CC= %.4f" %(title,
    cc =  corr.coefficient()
    return cc**2 / (2 - cc**2)

def get_peak(data, phase, site):
    data, phase = data.common_sets(phase)
    m = data.phase_transfer(phase_source=phase, deg=True).fft_map().apply_sigma_scaling().real_map_unpadded()
    return m.eight_point_interpolation(site)

def run(mtzin, pdbin):
    arrays = iotbx.mtz.object(mtzin).as_miller_arrays()
    #f_sel = filter(lambda x:x.is_amplitude_array() and x.anomalous_flag(), arrays)[0]
    i_obs = filter(lambda x:x.is_xray_intensity_array() and x.anomalous_flag(), arrays)[0]
    f_model = filter(lambda x:"F-model" in x.info().label_string() and x.anomalous_flag(), arrays)[0]
    dano_model = f_model.as_amplitude_array().anomalous_differences()
    phi_model = f_model.as_non_anomalous_array().phases(deg=True)
    phi_model = phi_model.customized_copy(data=phi_model.data()-90)

    print "  I-obs:", i_obs.info()
    print "F-model:", f_model.info()

    # Original DANO
    f_obs = i_obs.french_wilson()
    dano = f_obs.anomalous_differences()
    dano = dano.customized_copy(data=dano.data())

    dano = dano_model

    dano_max = flex.max(flex.abs(dano.data()))
    print "MAX dano=", dano_max

    # Really random DANO
    #random = dano.customized_copy(data=flex.random_int_gaussian_distribution(size=dano.size(), mu=0, sigma=10), sigmas=None)
    sigma_dano = flex.sum(flex.pow2(dano.data()-flex.mean(dano.data()))) / dano.size()
    print "SIGMA of DANO:", sigma_dano
    random = dano.customized_copy(data=flex.double(numpy.random.normal(size=dano.size(), loc=0, scale=sigma_dano)), sigmas=None)

    # DANO + noise
    dano_noises = []
    for p in (.01, .1, .2, .4, .5, .6, .7, .8, .9, 1):
        dano_noises.append((p, dano.customized_copy(data=flex.double(numpy.random.normal(loc=numpy.array(dano.data()), scale=p*dano_max)), sigmas=None)))

    mtz_ds = dano.as_mtz_dataset(column_root_label="DANO")
    mtz_ds.add_miller_array(random, column_root_label="RANDOM", column_types="F")

    for p, d in dano_noises:
        mtz_ds.add_miller_array(d, column_root_label="DANO_r%d"%(int(p*100)), column_types="F")

    mtz_ds.add_miller_array(dano_model, column_root_label="DANO-model", column_types="F")
    mtz_ds.add_miller_array(phi_model, column_root_label="PHI-model-90", column_types="P")

    mtz_ds.mtz_object().write(file_name="dano_noise_test.mtz")

    # Read Hg site from PDB
    from iotbx import file_reader
    #pdbin = "../LRE11_pLRE_Mg_MPD_refine_001_HG.pdb"
    pdb_in = file_reader.any_file(pdbin, force_type="pdb")
    pdb_symm = pdb_in.file_object.crystal_symmetry()
    hg_site = pdb_in.file_object.xray_structure_simple(crystal_symmetry=pdb_symm).sites_frac()[0]

    cc = get_cc(dano, dano, "model vs data")
    peak = get_peak(dano, phi_model, hg_site)
    print "model vs original DANO CC= %.5f map= %.2f"%(cc, peak)

    cc = get_cc(dano, random, "model vs random")
    peak = get_peak(random, phi_model, hg_site)
    print "model vs random DANO CC= %.5f map= %.2f"%(cc, peak)

    for p, d in dano_noises:
        cc = get_cc(dano, d, "")
        peak = get_peak(d, phi_model, hg_site)
        print "model vs DANO+%d%%noise CC= %.5f map= %.2f"%(int(p*100), cc, peak)


    # FFT commands
    print """\
fft hklin dano_noise_test.mtz mapout dano.ccp4 <<+
labin F1=DANO PHI=PHI-model-90
xyzlim asu
+

fft hklin dano_noise_test.mtz mapout random.ccp4 <<+
labin F1=RANDOM PHI=PHI-model-90
xyzlim asu
+

"""
    for p, d in dano_noises:
        print """\
fft hklin dano_noise_test.mtz mapout dano_r%(p)d.ccp4 <<+
labin F1=DANO_r%(p)d PHI=PHI-model-90
xyzlim asu
+
""" % dict(p=int(p*100))


if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2])
