from yamtbx.util import get_temp_local_dir
from yamtbx.util import safe_copy
from yamtbx.dataproc import eiger
from yamtbx.dataproc.XIO import XIO
import h5py
import os
import shutil
import time
import traceback
import re
import pysqlite2.dbapi2 as sqlite3
import glob

def get_nspots(dbfile, prefix): # prefix must include _
    re_file = re.compile(re.escape(prefix)+"[0-9]+\..{1,3}$")
    print prefix

    con = sqlite3.connect(dbfile, timeout=10)
    cur = con.cursor()

    for itrial in xrange(60):
        try:
            c = con.execute("select imgf,nspot from stats")
            file_spots = c.fetchall()
            #print file_spots
            return filter(lambda x: re_file.search(str(x[0])), file_spots)
        except sqlite3.DatabaseError:
            print traceback.format_exc()
            print "DB failed. retrying (%d)" % itrial
            time.sleep(1)
            continue
    return []
# get_nspots()

def create_hitsonly_h5(master_h5, dbfile, hits_h5, spots_min):
    start_time = time.time()

    prefix = os.path.basename(master_h5)[:-len("master.h5")]

    n_images = XIO.Image(master_h5).header["Nimages"]
    print " %d images in %s" % (n_images, master_h5)

    flag_ok = False
    for i in xrange(100):
        file_spots = get_nspots(dbfile, prefix)
        print "%6d spot-finding results out of %6d images" % (len(file_spots), n_images)
        if len(file_spots) >= n_images:
            flag_ok = True
            break
        time.sleep(6)

    if not flag_ok:
        print "ERROR: Timeout!"
        return

    hits = filter(lambda x: x[1]>=spots_min, file_spots)

    print "There are %6d hits (>= %d spots) out of %6d images" % (len(hits), spots_min, n_images)

    shutil.copyfile(master_h5, hits_h5)

    h = h5py.File(hits_h5, "a")

    for k in h["/entry/data"].keys():
        del h["/entry/data"][k]

    for name, nsp in sorted(hits):
        print " hit: %s %d" %(name, nsp)
        frameno = int(os.path.splitext(name[len(prefix):])[0])
        data = eiger.extract_data(master_h5, frameno, return_raw=True)
        grpname = "%s%.6d"%(prefix, frameno)
        rootgrp = h["/entry/data"]
        rootgrp.create_group(grpname)
        if data is not None:
            eiger.compress_h5data(h, "/entry/data/%s/data"%grpname, data, chunks=None, compression="bslz4")
            rootgrp["%s/data"%grpname].attrs["n_spots"] = nsp
        else:
            print "  error: data not found (%s)" % name

    h.close()

    org_files = filter(lambda x: os.path.exists(x), eiger.get_masterh5_related_filenames(master_h5))
    org_size = sum(map(lambda x: os.path.getsize(x), org_files))/1024.0**2
    new_size = os.path.getsize(hits_h5)/1024.0**2

    eltime = time.time() - start_time
    print "SIZE: %.2f MB (%d) to %.2f MB (saved %.2f MB; %.1f%%) took %.2f sec" % (org_size, len(org_files),
                                                                                   new_size, org_size-new_size,
                                                                                   new_size/org_size*100., eltime)

    # TODO if original file moved to OVERWRITTEN_, then discard this onlyhits.h5 file!!
    #      Because hit-finding results may be contaminated, and onlyhits.h5 may be no use!
    #      To check this consistency, need master h5 creation date or something.
    
# run()

def run(master_h5, master_h5_ctime, dbfile, tmpdir=None, spots_min=3, remove_files=False):
    """
    Download must be finished when this script started!!
    Everything in tmpdir will be removed!
    """

    os.stat_float_times(False)

    if tmpdir: assert os.path.isdir(tmpdir)

    try:
        assert master_h5.endswith("_master.h5")
        prefix = os.path.basename(master_h5)[:-len("master.h5")]
    
        # If tmpdir given, copy of master.h5 and data h5 files should be there
        if tmpdir: master_h5_in_tmp = os.path.join(tmpdir, os.path.basename(master_h5))
        else: master_h5_in_tmp = master_h5 # If not, just use this
    
        # If tmpdir given, just use there; otherwise create new one.
        if not tmpdir: tmpdir = get_temp_local_dir("hitsonly", min_gb=1)
        hits_h5 = os.path.join(tmpdir, prefix+"onlyhits.h5")
        print "tmpdir is %s" %tmpdir
        create_hitsonly_h5(master_h5_in_tmp, dbfile, hits_h5, spots_min)

        if not os.path.isfile(hits_h5):
            raise Exception("Generation of %s failed" % hits_h5)

        print "ctime_master_h5 =", os.path.getctime(master_h5)
        if os.path.getctime(master_h5) != master_h5_ctime:
            raise Exception("Master h5 file (%s, %d) is changed! Discarding this hitsonly h5" % (master_h5, master_h5_ctime))

        #shutil.move(hits_h5, os.path.join(os.path.dirname(master_h5), os.path.basename(hits_h5)))
        safe_copy(hits_h5, os.path.normpath(os.path.dirname(master_h5)), move=True)

    finally:
        files_in_tmp = glob.glob(os.path.join(tmpdir, "*"))
        if files_in_tmp:
            print "Removing %d files in tmpdir:" % len(files_in_tmp)
            for f in files_in_tmp: print " %s" % f

        shutil.rmtree(tmpdir)

# run()

if __name__ == "__main__":
    import sys
    import optparse

    parser = optparse.OptionParser(usage="usage: %prog [options] master.h5 master.h5-in-local dbfile",
                                   description="")
    parser.add_option("--min-spots", action="store", dest="spots_min", type=int, default=3)
    parser.add_option("--ctime-master", action="store", dest="ctime_master", type=int)
    parser.add_option("--tmpdir", action="store", dest="tmpdir", type=str)

    opts, args = parser.parse_args(sys.argv[1:])

    master_h5, dbfile = args

    os.stat_float_times(False)

    master_h5_ctime = opts.ctime_master
    if master_h5_ctime is None or master_h5_ctime < 0: master_h5_ctime = os.path.getctime(master_h5)

    if opts.tmpdir=="None": opts.tmpdir = None

    print "Command args:"
    print "  %s %s --min-spots=%s --ctime-master=%s --tmpdir=%s" %(master_h5, dbfile, opts.spots_min, opts.ctime_master, opts.tmpdir)

    run(master_h5, master_h5_ctime, dbfile, tmpdir=opts.tmpdir, spots_min=opts.spots_min)
