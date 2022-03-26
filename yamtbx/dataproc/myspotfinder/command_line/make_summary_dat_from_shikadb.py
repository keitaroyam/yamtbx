from __future__ import print_function
from __future__ import unicode_literals
from yamtbx.dataproc.myspotfinder.command_line.spot_finder_gui import Stat
from yamtbx.dataproc import bl_logfiles
import sqlite3
import pickle
import os
import numpy

def read_db(scanlog, dbfile):
    con = sqlite3.connect(dbfile, timeout=10)
    try:
        c = con.execute("select filename,spots from spots")
    except sqlite3.OperationalError:
        print("# DB Error (%s)" % dbfile)
        return None

    results = dict([(str(x[0]), pickle.loads(str(x[1]))) for x in c.fetchall()])
    ret = []

    slog = bl_logfiles.BssDiffscanLog(scanlog)
    slog.remove_overwritten_scans()

    for scan in slog.scans:
        for imgf, (gonio, gc) in scan.filename_coords:
            #print imgf, (gonio, gc) 
            stat = Stat()
            if imgf not in results: continue
            snrlist = [x[2] for x in results[imgf]["spots"]]
            stat.stats = (len(snrlist), sum(snrlist), numpy.median(snrlist) if snrlist else 0)
            stat.spots = results[imgf]["spots"]
            stat.gonio = gonio
            stat.grid_coord = gc
            stat.scan_info = scan
            stat.thumb_posmag = results[imgf]["thumb_posmag"]
            stat.params = results[imgf]["params"]
            stat.img_file = imgf # os.path.join(self.ctrlFrame.current_target_dir, imgf)
            ret.append((scan.get_prefix(), stat.img_file, stat))
    return ret
# read_db()

def make_dat(datout, data):
    ofs = open(datout, "w")

    print("prefix x y kind data filename", file=ofs)
    for kind in ("total_integrated_signal","median_integrated_signal","n_spots"):
        for prefix, f, stat in data:
            gc = stat.grid_coord
            if gc is None:
                x, y = "na", "na"
                print("Warning: gc is None! %s"%f)
            else:
                x, y = gc

            d = stat.stats[("n_spots","total_integrated_signal","median_integrated_signal").index(kind)]
            print(prefix, x, y, kind, d, f, file=ofs)

# make_dat()

def run(scanlog, dbfile, datout):
    data = read_db(scanlog, dbfile)
    if data:
        make_dat(datout, data)

def run_from_args(args):
    scan_dir = args[0]
    scanlog = os.path.join(scan_dir, "diffscan.log")
    shika_db = os.path.join(scan_dir, "_spotfinder", "shika.db")
    datout = os.path.join(scan_dir, "_spotfinder", "summary.dat")

    print("scanlog=", scanlog)
    print("shika_db=", shika_db)
    print("datout=", datout)

    if os.path.exists(datout):
        print("%s already exists." % datout)
        quit()

    run(scanlog, shika_db, datout)


if __name__ == "__main__":
    import sys
    #scanlog, shika_db = sys.argv[1:]

    run_from_args(sys.argv[1:])
