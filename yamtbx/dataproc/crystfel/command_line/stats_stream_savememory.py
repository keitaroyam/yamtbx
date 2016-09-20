"""
R
d<-read.table("stats.log",h=T)
library(ggplot2)
library(reshape2)
d.m<-melt(d,id=c("nframes","red"))
ggplot(d.m, aes(x=red, y=value)) + geom_point() + geom_line() + facet_grid(variable~.,scale="free")
"""

from yamtbx.dataproc import crystfel
from yamtbx.dataproc.crystfel.command_line import split_stats

from iotbx import crystal_symmetry_from_any
import iotbx.file_reader
import iotbx.phil
from cctbx import miller
from cctbx import crystal
from cctbx.array_family import flex
import time
import random
import sys
import re
import os
import bz2

master_params_str = """
dmin = 2.1
 .type = float
dmax = None
 .type = float
nsplit = 200
 .type = int
adu_cutoff = None
 .type = float
 .help = upper limit value of peak
anomalous = True
 .type = bool
nshells = 20
 .type = int
stop_after = None
 .type = int
start_at = 0
 .type = int
halve_method = *alternate random
 .type = choice(multi=False)
random_seed = 1234
 .type = int

#frame_scaling = False
# .type = bool
# .help = Perform frame scaling against merged data

reference {
  anohkl = None
   .type = path
  anolabel = None
   .type = str

  hkl = None
   .type = path
  label = None
   .type = str
}

output_prefix = None
 .type = str
"""

def get_reference_data(refdata, reflabel, assert_anomalous=True):
    arrays = iotbx.file_reader.any_file(refdata).file_server.miller_arrays
    array = filter(lambda x: x.info().label_string()==reflabel, arrays)[0]
    if assert_anomalous:
        assert array.anomalous_flag()
    return array.as_intensity_array()
# get_reference_data()

def merge_obs(indices, iobs, sel, symm, anomalous_flag, d_min, d_max):
    #indices = flex.miller_index(indices)
    #iobs = flex.double(iobs)
    if sel is not None and len(sel) > 0:
        sel = flex.bool(sel)
        indices = indices.select(sel)
        iobs = iobs.select(sel)

    miller_set = miller.set(crystal_symmetry=symm,
                            indices=indices,
                            anomalous_flag=anomalous_flag)

    array = miller.array(miller_set=miller_set,
                         data=iobs,
                         ).set_observation_type_xray_intensity()
    array = array.resolution_filter(d_min=d_min, d_max=d_max)

    ## New way
    #array = array.map_to_asu()
    #
    #merger = yamtbx_dataproc_crystfel_ext.merge_equivalents_crystfel()
    #merger.add_observations(array.indices(), array.data())
    #merger.merge()

    return array.merge_equivalents(algorithm="crystfel") # if sigmas is None, merge_equivalents_real() is used which simply averages.
# merge_obs()

def read_stream(stream, start_at=0):
    if stream.endswith(".bz2"):
        fin = bz2.BZ2File(stream)
    else:
        fin = open(stream)
    line = fin.readline()
    format_ver = re.search("CrystFEL stream format ([0-9\.]+)", line).group(1)
    print "# format version:", format_ver
    assert float(format_ver) >= 2.2 # TODO support other version
    chunk = None
    count = 0
    read_flag = False
    for l in fin:
        if "----- Begin chunk -----" in l:
            read_flag = True
            chunk = crystfel.stream.Chunk()
        elif "----- End chunk -----" in l:
            read_flag = False
            if chunk.indexed_by is not None:
                if start_at <= count: yield chunk
                count += 1
        elif read_flag:
            chunk.parse_line(l)

def show_split_stats(stream, nindexed, symm, params, anoref=None, ref=None, out_prefix="out", start_at=0):
    random.seed(params.random_seed)
    nstep = nindexed // params.nsplit

    out = open(out_prefix + ".dat", "w")

    out_byres = open(out_prefix + "_byres.dat", "w")
    print >>out_byres, ("%"+str(len(str(nindexed)))+"s")%"nframes",
    print >>out_byres, "dmin  dmax   red  cmpl acmpl  rsplit    rano   cc1/2   ccano  snr  rano/split CCanoref CCref"

    #formatstr_for_eachout = out_prefix+"_to%."+str(len(str(nindexed)))+"d_byres.dat"
    #indices1, iobs1, sel1 = [], [], []
    #indices2, iobs2, sel2 = [], [], []

    indices1, indices2 = flex.miller_index(), flex.miller_index()
    iobs1, iobs2 = flex.double(), flex.double()
    sel1, sel2 = flex.bool(), flex.bool()

    print >>out, "   nframes     red  Rsplit    Rano   CC1/2   CCano  snr  Rano/Rsplit CCanoref CCref"

    chunks = read_stream(stream, start_at)

    for i in xrange(params.nsplit):
        e = (i+1)*nstep
        if i == params.nsplit-1: e = nindexed

        read_all = (i == params.nsplit-1)

        slc = []
        for j, c in enumerate(chunks):
            sys.stderr.write("\rprocessed: %d" % (i*nstep+j))
            sys.stderr.flush()
            slc.append(c)
            if not read_all and j >= nstep-1:
                break
            if nindexed <= i*nstep+j:
                break

        if params.halve_method == "alternate":
            slc1 = slc[0::2]
            slc2 = slc[1::2]
        elif params.halve_method == "random":
            perm = range(len(slc))
            random.shuffle(perm)
            nhalf = len(slc)//2
            slc1 = [slc[r] for r in perm[:nhalf]]
            slc2 = [slc[r] for r in perm[nhalf:]]
        else:
            raise "Not-supported:", params.halve_method

        i1, o1, s1 = [], [], []
        for x in slc1:
            i1.extend(x.indices)
            o1.extend(x.iobs)
            if params.adu_cutoff is not None:
                s1.extend((a <= params.adu_cutoff for a in x.peak))

        i2, o2, s2 = [], [], []
        for x in slc2:
            i2.extend(x.indices)
            o2.extend(x.iobs)
            if params.adu_cutoff is not None:
                s2.extend((a <= params.adu_cutoff for a in x.peak))

        # Concatenate and sort
        indices1 = indices1.concatenate(flex.miller_index(i1))
        iobs1 = iobs1.concatenate(flex.double(o1))
        sel1 = sel1.concatenate(flex.bool(s1))
        indices2 = indices2.concatenate(flex.miller_index(i2))
        iobs2 = iobs2.concatenate(flex.double(o2))
        sel2 = sel2.concatenate(flex.bool(s2))
        perm1 = flex.sort_permutation(data=miller.index_span(indices1).pack(indices1), reverse=False)
        perm2 = flex.sort_permutation(data=miller.index_span(indices2).pack(indices2), reverse=False)
        indices1 = indices1.select(perm1)
        indices2 = indices2.select(perm2)
        iobs1 = iobs1.select(perm1)
        iobs2 = iobs2.select(perm2)
        if params.adu_cutoff is not None:
            sel1 = sel1.select(perm1)
            sel2 = sel2.select(perm2)

        # Merge
        m = merge_obs(indices1.concatenate(indices2), iobs1.concatenate(iobs2), sel1.concatenate(sel2), symm, anomalous_flag=params.anomalous, d_min=params.dmin, d_max=params.dmax)

        m1 = merge_obs(indices1, iobs1, sel1, symm, anomalous_flag=params.anomalous, d_min=params.dmin, d_max=params.dmax)
        m2 = merge_obs(indices2, iobs2, sel2, symm, anomalous_flag=params.anomalous, d_min=params.dmin, d_max=params.dmax)

        a1, a2 = m1.array().common_sets(m2.array())

        red = flex.sum(m.redundancies().data()) / m.array().data().size() if m.array().data().size()>0 else 0
        rsplit = split_stats.calc_rsplit(a1, a2)
        rano = split_stats.calc_rano(a1, a2) if params.anomalous else float("nan")
        cc = split_stats.calc_cc(a1, a2)
        ccano = split_stats.calc_ccano(a1, a2) if params.anomalous else float("nan")
        ccanoref = split_stats.calc_ccano(m.array(), anoref, take_common=True) if params.anomalous and anoref is not None else float("nan")
        ccref = split_stats.calc_cc(m.array(), ref, take_common=True) if ref is not None else float("nan")
        snr = flex.mean(m.array().data()/m.array().sigmas())
        print >>out, "%10d %7.2f %7.4f %7.4f % .4f % .4f %.2f %.4f % .4f % .4f" % (e, red, rsplit, rano, cc, ccano, snr, rano/rsplit, ccanoref, ccref)
        out.flush()

        byresolution_stats(params, a1, a2, m, anoref, ref, ("%"+str(len(str(nindexed)))+"d")%e, out_byres)
        print >>out_byres, "#"
        #open(formatstr_for_eachout%e, "w"))

    # Oveall stats by redundancy

# show_split_stats()

def byresolution_stats(params, a1, a2, m, anoref, ref, prefix, out):
    binner = a1.setup_binner(n_bins=params.nshells)

    for i_bin in binner.range_used():
        dmax, dmin = binner.bin_d_range(i_bin)
        sel1 = a1.select(binner.bin_indices() == i_bin)
        sel2 = a2.select(binner.bin_indices() == i_bin)
        selm = m.array().resolution_filter_selection(d_max=dmax, d_min=dmin)
        m_sel = m.array().select(selm)
        red = flex.sum(m.redundancies().select(selm).data()) / m_sel.data().size() if m_sel.data().size()>0 else 0
        cmpl = m_sel.completeness(d_max=dmax)*100.
        cmplano = m_sel.anomalous_completeness(d_max=dmax)*100. if params.anomalous else float("nan")
        rsplit = split_stats.calc_rsplit(sel1, sel2)
        rano = split_stats.calc_rano(sel1, sel2) if params.anomalous else float("nan")
        cc = split_stats.calc_cc(sel1, sel2)
        ccano = split_stats.calc_ccano(sel1, sel2) if params.anomalous else float("nan")
        ccanoref = split_stats.calc_ccano(m_sel, anoref, take_common=True) if params.anomalous and anoref is not None else float("nan")
        ccref = split_stats.calc_cc(m_sel, ref, take_common=True) if ref is not None else float("nan")
        snr = flex.mean(m_sel.data()/m_sel.sigmas())
        print >>out, "%s %5.2f %5.2f %5d %5.1f %5.1f %7.4f %7.4f % .4f % .4f %.2f %7.4f % .4f % .4f" % (prefix, dmax, dmin, red, cmpl, cmplano, rsplit, rano, cc, ccano, snr, rano/rsplit, ccanoref, ccref)
# byresolution_stats()

def run(streamin, symm_source, params):
    symm = crystal_symmetry_from_any.extract_from(symm_source)

    ref, anoref = None, None
    if None not in (params.reference.anohkl,params.reference.anolabel):
        anoref = get_reference_data(params.reference.anohkl, params.reference.anolabel)

    if None not in (params.reference.hkl,params.reference.label):
        ref = get_reference_data(params.reference.hkl, params.reference.label,
                                 assert_anomalous=False) # no need to be anomalous..
        if not params.anomalous and ref.anomalous_flag():
            ref = ref.average_bijvoet_mates()
    # how many frames indexed?
    nindexed = 0
    t = time.time()

    for l in bz2.BZ2File(streamin) if streamin.endswith(".bz2") else open(streamin):
        if l.startswith("indexed_by =") and l[l.index("=")+1:].strip() != "none":
            nindexed += 1
        if params.stop_after is not None and params.stop_after <= (nindexed-params.start_at):
            break

    print >> sys.stderr, "# nframes checked (%d). time:" % nindexed, time.time() - t

    nindexed -= params.start_at
    
    if params.output_prefix is None:
        params.output_prefix = "stats_%s" % os.path.splitext(os.path.basename(streamin))[0]

    show_split_stats(streamin, nindexed, symm, params, anoref=anoref, ref=ref, out_prefix=params.output_prefix, start_at=params.start_at)
# run()

if __name__ == "__main__":
    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    streamin, symm_source = cmdline.remaining_args

    run(streamin, symm_source, params)
