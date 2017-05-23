"""
(c) RIKEN 2017. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from mmtbx.scaling import absolute_scaling
from mmtbx.scaling.matthews import p_vm_calculator
from cctbx.array_family import flex
from cctbx import miller
from yamtbx.util.maths import weighted_correlation_coefficient
import numpy

def axis_label(v, orth_mat):
    h = numpy.dot(v, orth_mat)
    h /= numpy.linalg.norm(h)

    labs = ("a*", "b*", "c*")
    hstr = map(lambda x: "%+.2f%s" % (x[1],labs[x[0]]), enumerate(h))

    hstr = filter(lambda x: "0.00" not in x, hstr)
    hstr = map(lambda x: x.replace("1.00",""), hstr)
    hstr = "".join(hstr)
    if hstr[0]=="+": hstr = hstr[1:]
    return hstr
# axis_label()

def calc_principal_vectors(miller_array, log_out):
    """
    Reference: Aimless program
    Triclinic & Monoclinic: from B_cart
    Orthorhombic: along a*, b*, c*
    Trigonal, Hexagonal, Tetragonal: c* and a*b*-plane
    Cubic: no anisotropy

    return type: list of (vector, label, in-plane or not)
    """

    n_residues = p_vm_calculator(miller_array, 1, 0).best_guess

    aniso_scale_and_b = absolute_scaling.ml_aniso_absolute_scaling(miller_array=miller_array,
                                                                   n_residues=n_residues,
                                                                   n_bases=0)

    if aniso_scale_and_b.eigen_values[0] == 0:
        print >>log_out, "Error! Cannot determine B_cart"
        return

    b_cart = aniso_scale_and_b.b_cart
    print >>log_out, """ML estimate of overall B_cart:
 /  %5.2f %5.2f %5.2f \\
 |  %11.2f %5.2f |
 \\  %17.2f /
""" % (b_cart[0], b_cart[3], b_cart[4],
       b_cart[1], b_cart[5], b_cart[2])

    ev = aniso_scale_and_b.eigen_vectors

    orth_mat = numpy.array(miller_array.unit_cell().orthogonalization_matrix()).reshape(3,3)

    print >>log_out, "Eigenvalues/vectors:"
    for i, eg in enumerate(aniso_scale_and_b.eigen_values):
        v = ev[3*i:3*(i+1)]
        vs = ", ".join(map(lambda x: "% .4f"%x, v))
        vl = axis_label(v, orth_mat)
        print >>log_out, " %8.3f (%s) %s" % (eg, vs, vl)
    cs = miller_array.space_group().crystal_system()

    if cs == "Cubic":
        print >>log_out, "No anisotropy in this symmetry."
        return []
    elif cs == "Orthorhombic":
        return ([1.,0.,0.], "a*", False), ([0.,1.,0.], "b*", False), ([0.,0.,1.], "c*", False)
    elif cs in ("Triclinic", "Monoclinic"):
        return ((ev[:3], axis_label(ev[:3], orth_mat), False),
                (ev[3:6], axis_label(ev[3:6], orth_mat), False),
                (ev[6:], axis_label(ev[6:], orth_mat), False))
    else: # in "Tetragonal", "Hexagonal", "Trigonal", "Cubic"
        return ([0.,0.,1.], "c*", False), ([0.,0.,1.], "a*b*", True)
        
# calc_principal_vectors()

def calc_weighted_cc_half(unmerged, vstar, cos_limit):
    split_datasets = miller.split_unmerged(unmerged_indices=unmerged.indices(),
                                           unmerged_data=unmerged.data(),
                                           unmerged_sigmas=unmerged.sigmas(),
                                           seed=0)
    data_1 = split_datasets.data_1
    data_2 = split_datasets.data_2
    indices = split_datasets.indices

    if indices.size() == 0:
        return float("nan")

    if 0:
        return flex.linear_correlation(data_1, data_2).coefficient()

    pstar = unmerged.unit_cell().reciprocal_space_vector(indices).each_normalize()
    cos_theta = flex.abs(pstar.dot(vstar))
    weights = flex.abs(cos_theta - cos_limit)/(1 - cos_limit)

    return weighted_correlation_coefficient(data_1, data_2, weights)
# calc_weighted_cc_half()

def calc_stats_along_axes(miller_array, vstars, binner, angle=20, kind="cchalf"):
    assert kind in ("cchalf", "ioversigma")

    pstar = miller_array.unit_cell().reciprocal_space_vector(miller_array.indices()).each_normalize()

    ret = []

    for i, (vstar, label, plane_normal) in enumerate(vstars):
        vstar = tuple(vstar / numpy.linalg.norm(vstar))
        cos_theta = flex.abs(pstar.dot(vstar))

        if plane_normal:
            cos_limit = numpy.sin(numpy.deg2rad(angle))
            sel_in_cone = cos_theta <= cos_limit
        else:
            cos_limit = numpy.cos(numpy.deg2rad(angle))
            sel_in_cone = cos_theta >= cos_limit

        # DEBUG
        if 0:
            junk = miller_array.customized_copy(sigmas=miller_array.data()*0)
            junk.sigmas().set_selected(sel_in_cone, miller_array.data().select(sel_in_cone)/w_t.select(sel_in_cone))
            junk.as_mtz_dataset("JUNK").mtz_object().write("junk_%.3d.mtz"%i)
            continue

        in_cone = miller_array.select(sel_in_cone)

        if kind == "cchalf":
            overall = calc_weighted_cc_half(in_cone, vstar, cos_limit)
        else:
            overall = flex.mean(in_cone.data()/in_cone.sigmas()) if in_cone.size() > 0 else float("nan")

        ret.append([vstar, label, plane_normal, overall, [], []]) # nref, stat
        for i_bin in binner.range_used():
            d_max, d_min = binner.bin_d_range(i_bin)
            a_sel = in_cone.resolution_filter(d_max=d_max, d_min=d_min)

            if kind == "cchalf":
                val = calc_weighted_cc_half(a_sel, vstar, cos_limit)
            else:
                val = flex.mean(a_sel.data()/a_sel.sigmas()) if a_sel.size() > 0 else float("nan")

            ret[-1][-1].append(val)
            ret[-1][-2].append(a_sel.size())
    return ret
# calc_stats_along_axes()

if __name__ == "__main__":
    import sys
    import iotbx.file_reader

    f = iotbx.file_reader.any_file(sys.argv[1])
    try:
        array = f.file_server.get_xray_data(file_name=None,
                                            labels=None,
                                            ignore_all_zeros=True,
                                            parameter_scope="",
                                            prefer_anomalous=False,
                                            prefer_amplitudes=False)
    except:
        array = f.file_server.miller_arrays[0].as_amplitude_array()

    #for i, (v,pn) in enumerate(calc_principal_vectors(array)):
    calc_stats_along_axes(array, calc_principal_vectors(array, sys.stdout), None)
