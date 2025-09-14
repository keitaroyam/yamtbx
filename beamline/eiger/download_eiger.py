from __future__ import print_function
from __future__ import unicode_literals
import eigerclient
import h5py
import bitshuffle.h5
import tempfile
import time
import glob
import os
import shutil
import urllib.request, urllib.parse, urllib.error
import subprocess
import traceback
import re
import sys
from yamtbx.util import safe_copy
from yamtbx.util import get_temp_local_dir

#eiger_host = "192.168.163.204"
#default_tmpd = "/dev/shm"

now = lambda : time.strftime("%Y-%m-%d %H:%M:%S")

class ReduceMaster(object):
    def __init__(self, h5):
        self.h = h5
    # __init__()

    def remove_redundant(self, kinds):
        assert set(kinds).issubset(("flatfield","pixel_mask","trimbit"))

        detSp = self.h["/entry/instrument/detector/detectorSpecific"]
        for k in list(detSp.keys()):
            if k.startswith("detectorModule_"):
                for kk in kinds:
                    if kk not in detSp[k]: continue
                    print("  Removing", k+"/"+kk)
                    del detSp[k][kk]
    # remove_redundant()

    def find_large_dataset_visitor(self, path, obj):
        if type(obj) is h5py.Dataset and obj.size > 100:
            self.large_datasets.append(path)
    # find_large_dataset_visitor()

    def compress_large_datasets(self, compress):
        if not compress: return
        assert compress in ("bslz4", "bslz4_and_gzipshuf")

        self.large_datasets = []
        self.h.visititems(self.find_large_dataset_visitor)

        for path in self.large_datasets:
            print("  Compressing %s (size=%d)" % (path, self.h[path].size))
            data = self.h[path][:]
            del self.h[path]
            
            if compress=="bslz4":
                self.h.create_dataset(path, data.shape,
                                      compression=bitshuffle.h5.H5FILTER,
                                      compression_opts=(0, bitshuffle.h5.H5_COMPRESS_LZ4),
                                      chunks=None, dtype=data.dtype, data=data)
            elif compress == "bslz4_and_gzipshuf":
                if "pixel_mask" in path: # bslz4
                    self.h.create_dataset(path, data.shape,
                                          compression=bitshuffle.h5.H5FILTER,
                                          compression_opts=(0, bitshuffle.h5.H5_COMPRESS_LZ4),
                                          chunks=None, dtype=data.dtype, data=data)
                else: # gzip+shuf
                    self.h.create_dataset(path, data.shape,
                                          compression="gzip",shuffle=True,
                                          chunks=None, dtype=data.dtype, data=data)
            else:
                raise "Never reaches here."

    # compress_large_datasets()
# class ReduceMaster

def put_reversephi_info(h5):
    print("  Putting reverse-phi info")

    if "/entry/sample/depends_on" not in h5:
        h5["/entry/sample/depends_on"] = "/entry/sample/transformations/omega"

    grp = h5.require_group("/entry/sample/transformations")
    grp.attrs['NX_class'] = "NXtransformations"
    
    if "omega" not in grp:
        omega = h5["/entry/sample/goniometer/omega"]
        ds = grp.create_dataset("omega", omega.shape, dtype=omega.dtype, data=omega)
    else:
        ds = grp["omega"]

    fw_ver = h5["/entry/instrument/detector/detectorSpecific/eiger_fw_version"][()].decode()
    r = re.search(r"release-([0-9\.]+)", fw_ver)
    if r and [int(x) for x in r.group(1).split(".")] >= [2022,1,2]:
        print("   Newer firmware. Don't put omega vector")
    else:
        print("   Old firmware. Adding omega (-1,0,0)")
        ds.attrs['units'] = 'degree'
        ds.attrs['transformation_type'] = 'rotation'
        ds.attrs['vector'] = (-1., 0., 0.)
        ds.attrs['offset'] = (0., 0., 0.)
        ds.attrs['depends_on'] = '.'
# put_reversephi_info()

def fix_data_link(h5, bssid):
    keys = list(h5["/entry/data"].keys())
    filenames = {}
    for k in keys:
        filenames[k] = h5["/entry/data"].get(k, getlink=True).filename
        del h5["/entry/data"][k]

    for k in keys:
        assert filenames[k].startswith(bssid)
        f = filenames[k][len(bssid):]
        h5["/entry/data"][k] = h5py.ExternalLink(f, "/entry/data/data")

# fix_data_link()

def fix_omega(h5, omega_offset_by_trigger=None):
    ntrigger = h5["/entry/instrument/detector/detectorSpecific/ntrigger"][()]
    if ntrigger < 2: return

    nimages = h5["/entry/instrument/detector/detectorSpecific/nimages"][()]

    omega_increment = h5["/entry/sample/goniometer/omega_increment"][()]
    if omega_increment == 0: return

    print("  Need to fix omega values: omega_offset_by_trigger=%s" % omega_offset_by_trigger)

    omega = h5["/entry/sample/goniometer/omega"]
    omega_end = h5["/entry/sample/goniometer/omega_end"]
    omega_range_total = h5["/entry/sample/goniometer/omega_range_total"]

    if len(omega) != ntrigger * nimages:
        print("WARNING: len(omega) = %d != ntrigger*nimages = %d*%d" % (len(omega), ntrigger, nimages))

    print("   Original omega=", omega[:])

    if omega_offset_by_trigger is None:
        for i in range(1, ntrigger):
            omega[i*nimages:(i+1)*nimages] = omega[:nimages]
            omega_end[i*nimages:(i+1)*nimages] = omega_end[:nimages]
    else:
        for i in range(1, ntrigger):
            omega[i*nimages:(i+1)*nimages] += omega_offset_by_trigger*i
            omega_end[i*nimages:(i+1)*nimages] += omega_offset_by_trigger*i

    print("   Modified omega=", omega[:])

    # omega_range_total[()] =  # Leave it.
# fix_omega()

def move_original_files(prefix, wdir):
    files = glob.glob(os.path.join(wdir, "%s_data*.h5"%prefix))
    masterh5 = os.path.join(wdir, "%s_master.h5"%prefix)
    if os.path.isfile(masterh5): files.append(masterh5)
    hith5 = os.path.join(wdir, "%s_onlyhits.h5"%prefix)
    if os.path.isfile(hith5): files.append(hith5)

    if not files: return

    bdir = os.path.join(wdir, "OVERWRITTEN_%s"%time.strftime("%y%m%d-%H%M%S"))
    os.mkdir(bdir)
    print("Moving files with the same prefix to %s" % bdir)
    for f in files:
        os.rename(f, os.path.join(bdir, os.path.basename(f)))
# move_original_files()

def check_files(e, prefix, bssid):
    files = []

    try: files = e.fileWriterFiles()
    except: print("could not get file list")

    try:
        total_bytes = 0
        startt=time.time()
        for f in files:
            u = urllib.request.urlopen("http://%s/data/%s" % (e._host, f))
            total_bytes += int(u.info().getheaders("Content-Length")[0])
            u.close()

        print(now(), "On server: %d files %d bytes (%.3f sec)" % (len(files), total_bytes, time.time()-startt))
    except:
        pass
    pre = bssid + prefix if bssid else prefix

    files = [x for x in files if re.search(r"^"+re.escape(pre)+r"_master\.h5$", x) or re.search(r"^"+re.escape(pre)+r"_data_[0-9]+\.h5$", x)]
    return files
# check_files()

def modify_master(tmp, trg, bssid, omega_offset_by_trigger=None):
    print(" Modifying master..")
    h5 = h5py.File(tmp, "a")

    try:
        redmas = ReduceMaster(h5)

        # Remove unnecessary data in /entry/instrument/detector/detectorSpecific/detectorModule_*
        redmas.remove_redundant(("flatfield","pixel_mask","trimbit"))

        # Compress large data with bslz4 (for compatibility with Neggia plugin)
        #redmas.compress_large_datasets("bslz4")
        # Compress large data with bslz4 for pixel_mask and gzip+shuf for others (for compatibility with Neggia plugin and autoPROC)
        redmas.compress_large_datasets("bslz4_and_gzipshuf")
    except:
        print(traceback.format_exc())

    try:
        # Fix omega values if multi
        fix_omega(h5, omega_offset_by_trigger)
    except:
        print(traceback.format_exc())

    # Put reverse-phi info and fix the links to data.h5
    try:
        put_reversephi_info(h5)
    except:
        print(traceback.format_exc())

    if bssid: fix_data_link(h5, bssid)

    h5.close()

    # Run h5repack to clean up the removed space
    p = subprocess.Popen(["h5repack", tmp, trg], shell=False)
    p.wait()
# modify_master()

def download_files(e, files, wdir, bssid, tmpdir=None, omega_offset_by_trigger=None):
    """
    If bssid is not None, 'files' contains bssid+prefix_(master.h5|data_*.h5).

    When `tmpdir' is not None, download files to tmpdir once and then copy to destination.
    Files in tmpdir are kept, and will be used for hit-extraction.
    """

    failed_files = []

    for f in files:
        src = "http://%s/data/%s" % (e._host, f)
        #tmpfd, tmp = tempfile.mkstemp(prefix=f, dir=default_tmpd)
        #os.close(tmpfd)
        tmp = None
        if not bssid:
            trg = os.path.join(wdir, f)
        else:
            assert f.startswith(bssid)
            trg = os.path.join(wdir, f[len(bssid):])
        
        dl_failed = False
        for i in range(10):
            dl_failed = False
            print(now(), " donwloading %s (%dth try).." % (f, i+1))
            startt = time.time()
            try:
                u = urllib.request.urlopen(src)
                file_bytes = u.info().getheaders("Content-Length")[0]
                u.close()

                print(now(), "  file size on server %s = %s" % (f, file_bytes))
                tmpd = get_temp_local_dir("eigerdl", min_bytes=int(file_bytes), additional_tmpd=wdir)
                if tmpd is None:
                    print(now(), "  ERROR: no space available to download this file!")
                    dl_failed = True
                    break
                print(now(), "  to %s" % tmpd)
                tmp = os.path.join(tmpd, f)
                urllib.request.urlretrieve(src, tmp)
            except:
                print(traceback.format_exc())
                dl_failed = True # if success in last trial, dl_failed==False
                

            eltime = time.time() - startt

            if tmp and os.path.isfile(tmp):
                print(now(), "  done in %.2f sec. %.3f KB/s" % (eltime, os.path.getsize(tmp)/1024/eltime))
                break
            print() 

            # retry if downloading failed
            time.sleep(1)

        if dl_failed:
            print(now(), "  Download failed. Keeping on server: %s" % src)
            failed_files.append(f)

        if not tmp or not os.path.isfile(tmp): continue

        if f.endswith("_master.h5"):
            try:
                modify_master(tmp, tmp+"-fix.h5", bssid, omega_offset_by_trigger)
                if tmpdir: safe_copy(tmp+"-fix.h5", os.path.join(tmpdir, os.path.basename(trg)))
                startt = time.time()
                safe_copy(tmp+"-fix.h5", trg, move=True)
                eltime = time.time() - startt
                print(now(), "  local_copy done in %.2f sec. %.3f KB/s" % (eltime, os.path.getsize(trg)/1024/eltime))
            except:
                print(traceback.format_exc())
                if os.path.isfile(tmp):
                    if tmpdir: safe_copy(tmp, os.path.join(tmpdir, os.path.basename(trg)))
                    startt = time.time()
                    safe_copy(tmp, trg, move=True)
                    eltime = time.time() - startt
                    print(now(), "  local_copy done in %.2f sec. %.3f KB/s" % (eltime, os.path.getsize(trg)/1024/eltime))

            if os.path.isfile(tmp): os.remove(tmp)
        else:
            try:
                if tmpdir: safe_copy(tmp, os.path.join(tmpdir, os.path.basename(trg)))
                startt = time.time()
                safe_copy(tmp, trg, move=True)
                eltime = time.time() - startt
                print(now(), "  local_copy done in %.2f sec. %.3f KB/s" % (eltime, os.path.getsize(trg)/1024/eltime))
            except:
                dl_failed = True
                print(now(), "  local_copy failed. Keeping on server: %s" % src)
                failed_files.append(f)

        # delete from server
        if not dl_failed:
            e.fileWriterFiles(f, method="DELETE")

        os.rmdir(os.path.dirname(tmp))

    return failed_files
# download_files()

def check_all_files_done(prefix, wdir):
    missing_files = []

    master = os.path.join(wdir, "%s_master.h5" % prefix)
    if not os.path.isfile(master):
        print(" master file not exists")
        return [os.path.basename(master)]

    files = [master]

    h5 = h5py.File(master, "r")
    for k in list(h5["/entry/data"].keys()):
        dataf = os.path.join(wdir, "%s_%s.h5" % (prefix, str(k)))
        files.append(dataf)
        if not os.path.isfile(dataf):
            print(" data file %s not exists" % os.path.basename(dataf))
            missing_files.append(os.path.basename(dataf))

    if not missing_files:
        print(" all files exist:", " ".join([os.path.basename(x) for x in files]))

    return missing_files
# check_all_files_done()

def run(prefix, wdir, bssid, jobmode, opts):
    e = eigerclient.DEigerClient(host=opts.eiger_host)

    print("#####################################")
    print("Download %s_* to %s (hash %s jobmode %s)" % (prefix, wdir, bssid, jobmode))
    print("--hit-extract=%s --omega-offset-by-trigger=%s" % (opts.hit_extract, opts.omega_offset_by_trigger))
    print("#####################################")

    try: move_original_files(prefix, wdir) # FIXME!! (think if this is really ok)
    except: print(traceback.format_exc())
    
    last_dl_time = time.time()
    timeout = 60*60 if bssid else 60*10 # in seconds

    tmpdir = None
    #if opts.hit_extract and jobmode=="4": 
    #    from yamtbx.util import get_temp_local_dir
    #    tmpdir = get_temp_local_dir("forhitsonly", min_gb=1)

    i = 0
    failed_files = []

    while True:
        print(now(), "In %dth trial.." % (i+1))
        files = check_files(e, prefix, bssid)
        files = [x for x in files if x not in failed_files]

        if files:
            failed_files.extend(download_files(e, files, wdir, bssid, tmpdir, opts.omega_offset_by_trigger))
            last_dl_time = time.time()
        elif time.time() - last_dl_time > timeout:
            print(now(), "Download timeout!")
            sys.exit(1)

        sys.stdout.flush()
            
        missing_files = check_all_files_done(prefix, wdir)
        if bssid: missing_files = [bssid+x for x in missing_files]

        if not missing_files:
            print(now(), "Download %s to %s Success!" % (prefix, wdir))

            if opts.hit_extract and jobmode=="4":
                print("STARTING HIT EXTRACT")
                os.stat_float_times(False)
                master_h5 = os.path.join(wdir, "%s_master.h5" % prefix)
                ctime_master = os.path.getctime(master_h5)
                dbfile = os.path.join(os.path.dirname(master_h5), "_spotfinder", "shika.db")
                if not opts.no_sge:
                    args = ["qsub", "-q", opts.sge_q, "-cwd", "-N", "hitextract_%s"%prefix,
                            "-v", 'ctime_master=%d,master_h5="%s",tmpdir="%s",dbfile="%s"'%(ctime_master,
                                                                                      master_h5,
                                                                                      tmpdir,
                                                                                      dbfile),
                            "/blconfig/local_bss/yam/qsub_hit_extract_online.sh"]
                    print(" ".join(args))
                    p = subprocess.Popen(" ".join(args), shell=True, cwd=wdir, stdout=subprocess.PIPE)
                    p.wait()
                    print(p.stdout.read())
                    p = subprocess.Popen("qstat", shell=True, cwd=wdir, stdout=subprocess.PIPE)
                    p.wait()
                    print(p.stdout.read())
                else:
                    args = ["ssh", opts.ssh_host, """\
"cd '%s'; env ctime_master=%d master_h5='%s' tmpdir='%s' dbfile='%s' bash /oys/xtal/yamtbx/bl32xu/eiger/qsub_hit_extract_online.sh" > hitextract_%s.log 2>&1 & \
""" % (wdir, ctime_master, master_h5, tmpdir, dbfile, prefix)]
                    print(" ".join(args))
                    p = subprocess.Popen(" ".join(args), shell=True, cwd=wdir) # as background job

                """
                import hit_extract_online
                os.stat_float_times(False)
                master_h5 = os.path.join(wdir, "%s_master.h5" % prefix)
                hit_extract_online.run(master_h5=master_h5,
                                       master_h5_ctime=os.path.getctime(master_h5),
                                       master_h5_in_tmp=master_h5,
                                       dbfile=os.path.join(os.path.dirname(master_h5), "_spotfinder", "shika.db"),
                                       spots_min=0)
                """
            sys.exit()
        elif set(missing_files).issubset(set(failed_files)):
            print(now(), "Error occurred during downloading following files. Check the logfile!")
            for f in failed_files: print("  %s" % f)
            sys.exit(1)


        if not files:
            time.sleep(3)
        i += 1

    print(now(), "Download Failed!!!!!")
    sys.exit(1)
# run()

if __name__ == "__main__":
    import sys
    import optparse

    parser = optparse.OptionParser(usage="usage: %prog [options] prefix wdir [bssid] [jobmode]",
                                   description="")
    parser.add_option("--hit-extract", action="store_true", dest="hit_extract")
    parser.add_option("--sge-q", action="store", type=str, dest="sge_q", default="all.q@bss32xu")
    parser.add_option("--no-sge", action="store_true", dest="no_sge", help="when qsub is not available")
    parser.add_option("--ssh-host", action="store", type=str, dest="ssh_host", help="for hit-extract (when --no-sge)")
    parser.add_option("--omega-offset-by-trigger", action="store", type=float, dest="omega_offset_by_trigger", default=0)
    parser.add_option("--eiger-host", action="store", type=str, dest="eiger_host", default="192.168.163.204")

    opts, args = parser.parse_args(sys.argv[1:])


    pre, wdir = args[:2]
    bssid = None if len(args) < 3 else args[2]
    jobmode = None if len(args) < 4 else args[3]
    run(pre, wdir, bssid, jobmode, opts)
