import os
import cPickle as pickle
import math
from yamtbx.dataproc import crystfel
import iotbx.phil
from cctbx.array_family import flex
from cctbx import uctbx
from mmtbx.scaling.absolute_scaling import ml_iso_absolute_scaling


master_params_str = """\
streamin = None
 .type = path
 .help = Input file.
pklout = None
 .type = path
 .help = Output pickle file, which is needed by sort_stream.py
datout = None
 .type = path
 .help = Output dat file, which is useful for plotting

space_group = P222
 .type = str
 .help = used for merging data
ref_cell = 48.07, 77.45, 84.77, 90., 90., 90.
 .type = floats(size=6)
 .help = used for comparing unit cells
n_residues = 308
 .type = int
 .help = used for Wilson B calculation

stats = *reslimit *ioversigma *resnatsnr1 *pr *wilsonb *abdist
 .type = choice(multi=True)
 .help = "statistics. reslimit: diffraction_resolution_limit reported by CrystFEL; ioversigma: <I/sigma(I)> in each pattern; resnatsnr1: resolution where <I/sigma(I)>_bin drops below 1; pr: profile_radius reported by CrystFEL; wilsonb: ML-estimate of Wilson B value; abdist: Andrews-Bernstein distance (unit cell dissimilarity)"
"""

def make_G6(uc):
    # taken from cctbx/uctbx/determine_unit_cell/target_uc.py
    uc = uctbx.unit_cell(uc).niggli_cell().parameters()
    """ Take a reduced Niggli Cell, and turn it into the G6 representation """
    a = uc[0] ** 2
    b = uc[1] ** 2
    c = uc[2] ** 2
    d = 2 * uc[1] * uc[2] * math.cos(uc[3])
    e = 2 * uc[0] * uc[2] * math.cos(uc[4])
    f = 2 * uc[0] * uc[1] * math.cos(uc[5])
    return [a, b, c, d, e, f]


def set_chunk_stats(chunk, stats, stat_choice, n_residues=None, ref_cell=None, space_group=None):
    if "reslimit" in stat_choice: stats["reslimit"].append(chunk.res_lim)
    else: stats["reslimit"].append(float("nan"))

    if "pr" in stat_choice: stats["pr"].append(chunk.profile_radius)
    else: stats["pr"].append(float("nan"))

    if "ioversigma" in stat_choice or "resnatsnr1" in stat_choice:
        iobs = chunk.data_array(space_group, False)
        iobs = iobs.select(iobs.sigmas()>0).merge_equivalents(use_internal_variance=False).array()
        binner = iobs.setup_binner(auto_binning=True)

        if "ioversigma" in stat_choice:
            stats["ioversigma"].append(flex.mean(iobs.data()/iobs.sigmas()))
        else:
            stats["ioversigma"].append(float("nan"))

        if "resnatsnr1" in stat_choice:
            res = float("nan")
            for i_bin in binner.range_used():
                sel = binner.selection(i_bin)
                tmp = iobs.select(sel)
                if tmp.size() == 0: continue
                sn = flex.mean(tmp.data()/tmp.sigmas())
                if sn <= 1:
                    res = binner.bin_d_range(i_bin)[1]
                    break

            stats["resnatsnr1"].append(res)
        else:
            stats["resnatsnr1"].append(float("nan"))
    else:
        stats["ioversigma"].append(float("nan"))
        stats["resnatsnr1"].append(float("nan"))

    if "abdist" in stat_choice:
        from cctbx.uctbx.determine_unit_cell import NCDist
        G6a, G6b = make_G6(ref_cell), make_G6(chunk.cell)
        abdist = NCDist(G6a, G6b)
        stats["abdist"].append(abdist)
    else:
        stats["abdist"].append(float("nan"))
    
    if "wilsonb" in stat_choice:
        iso_scale_and_b = ml_iso_absolute_scaling(iobs, n_residues, 0)
        stats["wilsonb"].append(iso_scale_and_b.b_wilson)
    else:
        stats["wilsonb"].append(float("nan"))
# set_chunk_stats()

def run(params):
    import sys

    if params.pklout is None: params.pklout = os.path.basename(params.streamin)+".pkl"
    if params.datout is None: params.datout = os.path.basename(params.streamin)+".dat"

    ifs = open(params.streamin)
    ofs_dat = open(params.datout, "w")
    ofs_dat.write("file event indexed_by reslimit ioversigma resnatsnr1 pr wilsonb abdist a b c al be ga\n")

    tell_p, tell_c = 0, 0
    read_flag = False
    chunk = None

    chunk_ranges = []
    stats = dict(reslimit=[],
                 ioversigma=[],
                 resnatsnr1=[],
                 abdist=[],
                 pr=[],
                 wilsonb=[],
                 )

    while True:
        l = ifs.readline()
        if l == "": break
        tell_p, tell_c = tell_c, ifs.tell()

        if "----- Begin chunk -----" in l:
            if read_flag: del chunk_ranges[-1] # in case chunk is not properly closed
            read_flag = True
            chunk = crystfel.stream.Chunk()
            chunk_ranges.append([tell_p+1,0])
        elif read_flag and "----- End chunk -----" in l:
            read_flag = False
            if chunk.indexed_by is not None:
                chunk_ranges[-1][1] = tell_c
                sys.stdout.write("%.6d processed\r" % len(chunk_ranges))
                set_chunk_stats(chunk, stats, params.stats,
                                n_residues=params.n_residues,
                                ref_cell=params.ref_cell,
                                space_group=params.space_group)
                ofs_dat.write("%s %s %s %.3f %.3f %.3f %.3e %.3f %.3f "%(chunk.filename, chunk.event, chunk.indexed_by, chunk.res_lim, 
                                                                         stats["ioversigma"][-1], stats["resnatsnr1"][-1], stats["pr"][-1], stats["wilsonb"][-1], stats["abdist"][-1]))
                ofs_dat.write("%.3f %.3f %.3f %.2f %.2f %.2f\n" % chunk.cell)
                #ofs_dat.flush()
            else:
                del chunk_ranges[-1]

        elif read_flag:
            try: chunk.parse_line(l)
            except:
                print "\nError in reading line: '%s'" % l
                read_flag = False
                
    if read_flag: #  Unfinished chunk
        print "\nWarning: unclosed chunk."
        del chunk_ranges[-1]
    
    stats["chunk_ranges"] = chunk_ranges
    pickle.dump(stats, open(params.pklout,"w"), -1)

    print
    print
    print "Use sort_stream.py %s %s ioversigma-"%(params.streamin, params.pklout)
    print "The suffix -: sort by descending order +: sort by ascending order"
# run()

if __name__ == "__main__":
    import sys

    if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
        print "All parameters:\n"
        iotbx.phil.parse(master_params_str).show(prefix="  ", attributes_level=1)
        quit()

    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args

    for arg in args:
        if os.path.isfile(arg) and params.streamin is None:
            params.streamin = arg

    if params.streamin is None:
        print "Give stream file"
        quit()

    run(params)
