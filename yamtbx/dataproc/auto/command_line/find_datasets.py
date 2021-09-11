"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals
import iotbx.phil
from libtbx.utils import multi_out
from yamtbx.dataproc.bl_logfiles import BssJobLog
from yamtbx.dataproc import dataset
from yamtbx.dataproc.command_line import spot_finder_gui
from yamtbx.dataproc.command_line import spot_finder_generic
import os
import pickle
import traceback

master_params_str = """
logtype = *sp8 pf
 .type = choice(multi=True)
pklout = None
 .type = path
find_spots = False
 .type = bool
use_cuda = False
 .type = bool
help = False
 .type = bool
"""

def check_files_exist(template, n_images, datadir, check_compressed=True):
    filenames = template_to_filenames(os.path.basename(template), 1, n_images)
    filenames = [os.path.join(datadir, s) for s in filenames]
    found = 0
    for f in filenames:
        if os.path.isfile(f):
            found += 1
        else:
            for ext in (".bz2", ".gz", ".xz"):
                if os.path.isfile(f+ext):
                    found += 1
                    break
    return found
# check_files_exist()

def test_gonio_coords_equal(centers):
    if len(centers) <= 1:
        return True

    cen0 = centers[0]
    for cen in centers[1:]:
        if cen != cen0:
            return False
    return True
# test_gonio_coords_equal()

def search_spots(img_template, n_images, config_mgr):
    filenames = dataset.template_to_filenames(img_template, 1, n_images)

    img_files, idxes = [], []
    for i, imgf in enumerate(filenames):
        if os.path.isfile(imgf):
            img_files.append(imgf)
            idxes.append(i+1)

    # Assuming the images are all on the same configuration.
    detkey = config_mgr.get_key_by_img(img_files[0])
    params = config_mgr.get_params_by_key(detkey)
    
    file_stats = spot_finder_generic.run(img_files, params)
    return [(idx, f, stat) for idx, (f, stat) in zip(idxes, file_stats)]
# search_spots()

def run(params, topdirs):
    import sys

    LogClass = dict(sp8=BssJobLog,
                    pf=None)

    out_root = os.path.dirname(params.pklout)
    out = multi_out()
    out.register("log", open(os.path.join(out_root, "find_datasets.log"), "w"), atexit_send_to=None)
    out.register("stdout", sys.stdout)

    if params.find_spots:
        shikalog = open(os.path.join(out_root, "shika.log"), "w")
        shikalog.write("prefix idx nspots\n")

    config_mgr = spot_finder_gui.ConfigManager(use_cuda=params.use_cuda)

    logobjects = []

    for topdir in topdirs:
        print("Looking into %s\n" % topdir, file=out)
        for root, dirnames, filenames in os.walk(topdir, followlinks=True):
            logs = [x for x in filenames if x.endswith(".log")]
            for log in logs:
                log = os.path.join(root, log)

                for logtype in params.logtype:
                    if 1:#try:
                        try:
                            logobj = LogClass[logtype](log)
                        except:
                            print(traceback.format_exc(), file=out)
                            print("\nException raised when parsing %s\n"%log, file=out)
                            raise

                        logobj.jobs = [x for x in logobj.jobs if x.job_mode != "XAFS"]
                        if len(logobj.jobs) > 0:
                            print("Found job log:", log, file=out)
                            spots = []
                            for job in logobj.jobs:
                                if params.find_spots:
                                    spots.append(search_spots(os.path.join(root, os.path.basename(job.filename)),
                                                              job.n_images, config_mgr))
                                    for idx, f, stat in spots[-1]:
                                        if stat is None: continue
                                        n_spots = stat.spots.get_n_spots("hi_pass_resolution_spots")
                                        shikalog.write("%s %6d %4d\n" % (os.path.join(root,
                                                                                      os.path.basename(job.filename)),
                                                                         idx, n_spots))
                                        shikalog.flush()
                                else:
                                    spots.append([])
                            logobjects.append([log, logobj, spots])
                    #except:
                    #    pass

    print(file=out)

    logobjects.sort(key=lambda x:x[0])

    for log, logobj, spots in logobjects:
        print(log)
        for job, spt in zip(logobj.jobs, spots):
            jobtype = "?"
            if job.advanced_centering == {} or test_gonio_coords_equal(job.advanced_centering.get("centers",[])):
                jobtype = "Single  " 
            elif job.advanced_centering.get("mode", "") == "vector_centering":
                jobtype = "Helical "
            elif job.advanced_centering.get("mode", "") == "multiple_centering":
                jobtype = "MultiCen"
            else:
                print("WARNING:: WHY REACH HERE?", job.advanced_centering, file=out)
                
            osc_range = job.osc_end - job.osc_start
            #nfound = check_files_exist(job.filename, job.n_images, os.path.dirname(log))
            nfound = len(dataset.find_existing_files_in_template(job.filename, 1, job.n_images,
                                                                 datadir=os.path.dirname(log),
                                                                 check_compressed=True))

            print(" %s osc_step=%6.3f osc_range=%5.1f found=%d/%d %s beam_size=%s" % (job.beamline, job.osc_step, osc_range,
                                                   nfound, job.n_images, jobtype,
                                                   "x".join(["%.1f"%x for x in job.beam_size])), file=out)
            if params.find_spots:
                n_spots = [stat.spots.get_n_spots("hi_pass_resolution_spots") for idx, f, stat in spt if stat is not None]
                gt20 = len([x for x in n_spots if x>=20])
                print(" spots < 5A: %d..%d (>=20 spots: %d frames)" % (min(n_spots), max(n_spots), gt20), file=out)

        print(file=out)

    pickle.dump(logobjects, open(params.pklout, "wb"), -1)
# run()

if __name__ == "__main__":
    import sys

    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args
    topdirs = []

    if "-h" in args: params.help = True

    if params.help:
        print("Parameter syntax:\n")
        iotbx.phil.parse(master_params_str).show(prefix="  ")
        print()
        print("Usage: find_datasets.py [data-root-dir/ | directory-list-file] [find_spots=true]")
        quit()

    for arg in args:
        if os.path.isdir(arg) and topdirs == []:
            topdirs.append(arg)
        if os.path.isfile(arg):
            topdirs.extend([x for x in [x.strip() for x in open(arg)] if os.path.isdir(x)])

    if topdirs == []:
        topdirs = [os.getcwd()]

    if params.pklout is None:
        params.pklout = "logobjects.pkl"

    run(params, topdirs)

    print()
    print("Want to prepare xds runs?")
    print("prep_xds_runs.py %s [anomalous=true]" % params.pklout)
