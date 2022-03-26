from __future__ import print_function
from __future__ import unicode_literals
from cctbx import miller
from cctbx import crystal
from cctbx.array_family import flex

def read_x(xfile, take_full=False):
    """
    according to http://www.hkl-xray.com/denzo-output

Regular (film output file, Denzo _ip) imaging plate output format (i.e. the format you would you to read in in FORTRAN) is:

format (3i4,i2,2f8.0,f7.0,f6.0,f6.0,2f7.0,f6.0,f8.0)
h,k,l
Flag 0 - full 1 - partial
Intensity (F**2) by profile fitting
s of intensity (F**2)
c2 of profile fitting
Intensity (F**2) by profile summation
Cosine of incidence angle at detector
Predicted pixel position of spot centroid (slow, fast) directions
Lorentz, polarization, obliquity combined factor
Strength of averaged profile in arbitrary units
    """

    lines = open(xfile).readlines()

    indices = flex.miller_index()
    data = flex.double()
    sigmas = flex.double()
    full_flags = flex.bool()
    cell, sg_str = None, None

    read_hkl = True

    for l in lines[5:]:
        if l.startswith(" 999") and len(l.split()) == 3:
            read_hkl = False

        if read_hkl:
            hkl = list(map(int, (l[0:4], l[4:8], l[8:12])))
            ispart = int(l[13])
            iobs, sigma = list(map(float, (l[14:22], l[22:30])))

            indices.append(hkl)
            data.append(iobs)
            sigmas.append(sigma)
            full_flags.append(ispart==0)
        else:
            if l.startswith("unit cell"):
                cell = list(map(float, l.split()[2:]))
            elif l.startswith("space group"):
                sg_str = l.split()[-1]

    #sg_str="p1"
    symm = crystal.symmetry(unit_cell=cell, space_group=sg_str)

    array = miller.array(miller_set=miller.set(crystal_symmetry=symm,
                                               indices=indices,
                                               anomalous_flag=True),
                         data=data, sigmas=sigmas).set_observation_type_xray_intensity()

    if take_full:
        array = array.select(full_flags)

    return array
# read_x()

def run(x1, x2, n_bins):
    hemispheres = []
    merged_data = []

    for x in (x1, x2):
        print("Processing %s" % x)
        print("====================")
        data = read_x(x, take_full=False) # if you need only full, change this to True.
        data.crystal_symmetry().show_summary()
        merge = data.merge_equivalents(use_internal_variance=False)
        merge.show_summary()
        print()
        array = merge.array()
        array = array.select(array.sigmas() > 0)
        merged_data.append(array)

        # separate + and -
        matches = array.match_bijvoet_mates()[1] # returns asu and matches
        sel_p = matches.pairs_hemisphere_selection("+")
        sel_p.extend(matches.singles("+"))
        sel_m = matches.pairs_hemisphere_selection("-")
        sel_m.extend(matches.singles("-"))

        hemispheres.append([array.select(sel_p, anomalous_flag=False),
                            array.select(sel_m, anomalous_flag=False).map_to_asu()])

    print("x1: merged=%d (+)=%d (-)=%d" % (merged_data[0].size(),
                                           hemispheres[0][0].size(), hemispheres[0][1].size()))
    print("x2: merged=%d (+)=%d (-)=%d" % (merged_data[1].size(),
                                           hemispheres[1][0].size(), hemispheres[1][1].size()))

    # for sigma calculation (new_sigma^2 = sigma1^2 + sigma2^2)
    additive_sigmas = lambda x, y: flex.sqrt(flex.pow2(x.sigmas()) + flex.pow2(y.sigmas()))

    # calculate data1(+) - data2(-)
    """ # for debug
    for i, x in enumerate(hemispheres[0][0]):
        print x
        if i > 10: break
    print
    for i, x in enumerate(hemispheres[1][0]):
        print x
        if i > 10: break
    """
    h1p, h2m = hemispheres[0][0].common_sets(hemispheres[1][1])
    h1p_h2m = h1p.customized_copy(data=h1p.data() - h2m.data(),
                                  sigmas=additive_sigmas(h1p, h2m))
    print(h1p_h2m.size())
    #for x in h1p_h2m: print x

    # calculate data2(+) - data1(-)
    h2p, h1m = hemispheres[1][0].common_sets(hemispheres[0][1])
    h2p_h1m = h2p.customized_copy(data=h2p.data() - h1m.data(),
                                  sigmas=additive_sigmas(h2p, h1m))
    print(h2p_h1m.size())
    print() 
    #for x in h2p_h1m: print x

    # concatenate data1(+) - data2(-) and data2(+) - data1(-)
    dano_tmp = h1p_h2m.concatenate(h2p_h1m)
    merge = dano_tmp.merge_equivalents(use_internal_variance=False)
    print("Merging stats of (+)-(-) data")
    print("=============================")
    merge.show_summary()
    dano = merge.array()
    
    print("num_dano=", dano.size())
    print()

    # process with binning
    dano.setup_binner(n_bins=n_bins)
    binner = dano.binner()

    print("Result:")
    print("   dmax    dmin  nrefs  dano")
    for i_bin in binner.range_used():
        # selection for this bin. sel is flex.bool object (list of True of False)
        sel = binner.selection(i_bin)
        count = binner.count(i_bin)
        # take mean of absolute value of anomalous differences in a bin
        if count > 0:
            bin_mean = flex.mean(flex.abs(dano.select(sel).data()))
        else:
            bin_mean = float("nan")
        d_max, d_min = binner.bin_d_range(i_bin)
        print("%7.2f %7.2f %6d %.2f" % (d_max, d_min, count, bin_mean))
# run()

if __name__ == "__main__":
    import sys
    x1, x2 = sys.argv[1:3]
    n_bins = int(sys.argv[3])

    run(x1, x2, n_bins)
