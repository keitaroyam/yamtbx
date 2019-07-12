
"""
Reference:
 Python Multiprocessing with ZeroMQ
 http://taotetek.net/2011/02/02/python-multiprocessing-with-zeromq/
"""

import iotbx.phil
import libtbx.phil

import os
import stat
import time
import datetime
import getpass
import zmq
import re
import Queue
import collections
#import sqlite3
import pysqlite2.dbapi2 as sqlite3
import threading
import traceback
import numpy
import cPickle as pickle
import hashlib
from PIL import Image
from multiprocessing import Process
#import inotify.adapters # use yamtbx.python -mpip install inotify

from yamtbx.dataproc.myspotfinder import shikalog
from yamtbx.dataproc.myspotfinder import config_manager
from yamtbx.dataproc.myspotfinder import spot_finder_for_grid_scan
from yamtbx.dataproc import bl_logfiles
from yamtbx.dataproc import eiger
from yamtbx import util

#DEBUG for inotify
#import logging
#logger = logging.getLogger("inotify.adapters")
#logger.setLevel(logging.DEBUG)
#handlers = logging.StreamHandler()
#handlers.setLevel(logging.DEBUG)
#handlers.setFormatter(logging.Formatter("%(asctime)-15s %(levelname)s : %(message)s"))
#logger.addHandler(handlers)


master_params_str = """\
topdir = None
 .type = path
 .help = Root directory
bl = 32xu 41xu 26b2 44xu 45xu
 .type = choice(multi=False)
 .help = Choose beamline where you start SHIKA
date = "today"
 .type = str
 .help = Data collection date ("today" or %Y-%d-%m format)
blconfig = None
 .type = path
 .help = Override default blconfig path (/isilon/blconfig/bl$bl/)
nproc = 4
 .type = int
 .help = Number of processors used for spot finding
ports = 5557,5558,5559
 .type = ints(size=3,value_min=1024,value_max=49151)
 .help = Port numbers used by ZeroMQ.
dbdir = /isilon/cluster/log/shika/db
 .type = path
 .help = location to write sqlite3 db file.
logroot = /isilon/cluster/log/shika/
 .type = path
mode = *eiger_streaming bsslog zoo watch_ramdisk
 .type = choice(multi=False)
env = *oys ppu
 .type = choice(multi=False)
 .help = Excetution environment
eiger_host = "192.168.163.204"
 .type = str
 .help = "EIGER hostname or ip-address"

#incomplete_file_workaround = 0
# .type = float
# .help = wait given seconds after detecting new image
force_ssh_from = None
 .type = str
 .help = Users must not change this parameter.
only_check_in_last_hours = 1
 .type = float
 .help = "Only check diffscan.log modified during the last specified hours"
ramdisk_walk_interval = 2
 .type = float
"""

params = None

def retry_until_success(f, arg=None):
    args = (arg,) if arg is not None else ()
    return util.retry_until_noexc(f, args, ntry=30, outf=shikalog.warning)

class DiffScanManager:
    def __init__(self):
        self.clear()
    # __init__()

    def clear(self):
        self.scanlog = {} # directory: BssDiffscanLog object
        self.found_imgs = set()
        self.finished = {} # {filename:timestamp}; list of filenames of which analysis was completed
    # clear()

    def add_scanlog(self, slog):
        slog = os.path.abspath(slog)
        self.scanlog[os.path.dirname(slog)] = bl_logfiles.BssDiffscanLog(slog)
    # add_scanlog()

    def add_dir(self, slogdir):
        self.add_scanlog(os.path.join(slogdir, "diffscan.log"))
    # add_dir()

    def update_scanlogs(self):
        for logdir, slog in self.scanlog.items():
            if os.path.isfile(slog.scanlog):
                slog.parse()
            else:
                shikalog.error("diffraction scan log is not found!: %s" %slog.scanlog)
                continue

            # if no update since all images processed
            mtime = os.path.getmtime(slog.scanlog)
            if mtime == self.finished.get(slog.scanlog, -1): continue

            # Update 'processed files' using database
            dbfile = os.path.join(logdir, "_spotfinder", "shika.db")

            for _ in xrange(10):
                try:
                    if not os.path.exists(os.path.dirname(dbfile)): os.mkdir(os.path.dirname(dbfile))
                    con = sqlite3.connect(dbfile, timeout=30)
                    break
                except sqlite3.OperationalError:
                    shikalog.warning("Connecting to %s failed. Retrying" % dbfile)
            cur = con.cursor()

            # status TABLE
            #c = cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='status';")
            c = retry_until_success(cur.execute, "SELECT name FROM sqlite_master WHERE type='table' AND name='status';")
            if c.fetchone() is not None:
                #cur.execute("select filename from status")
                retry_until_success(cur.execute, "select filename from status")
                processed_files = map(lambda x:os.path.join(logdir, x[0]), cur.fetchall())
                print "debug::",processed_files 
                self.found_imgs.update(processed_files)
                # check if all images are processed
                if len(processed_files) == sum(map(lambda x: len(x.filename_idxes), slog.scans)):
                    shikalog.info("All %d images in %s are already processed. Will not check unless diffscan.log updated"%(len(processed_files), slog.scanlog))
                    self.finished[slog.scanlog] = mtime

            #con.commit()
            retry_until_success(con.commit)
            con.close()
        
    # update_scanlogs()

    def get_unprocessed_images(self, env=None):
        ret = []
        self.update_scanlogs()

        for slogdir, slog in self.scanlog.items():
            for scan in slog.scans:
                fcs = map(lambda x: (os.path.join(slogdir, x[0]), x[1]), scan.filename_idxes)
                #print "fix=", fcs
                #if env == "ppu": fcs_proxy = map(lambda x: (re.sub("^/isilon/users/", "/ramdisk/", x[0]), x[1]), fcs)
                if env == "ppu": f_mod = lambda x: re.sub("^/isilon/users/", "/ramdisk/", x)
                else: f_mod = lambda x: x
                unproc = filter(lambda x: x[0] not in self.found_imgs and os.path.isfile(f_mod(x[0])), fcs)
                ret.extend(map(lambda x:x+(scan,), unproc))

        self.found_imgs.update(map(lambda x: x[0], ret))
        return ret # (filename, idx, scan object)
    # get_unprocessed_images()

    def remove_found(self, files): # when user wants to recalculate..
        self.found_imgs.difference_update(files)
    # remove_found()

    def needs_to_be_processed(self, filename):
        """
        Check if the given file needs to be processed.
        No need to process if
        - not included in diffscan.log
        - first image in row (only if BSS creates such non-sense files)
        """

        scaninfo = self.get_scan_info(filename)
        if scaninfo is None:
            return False

        # return True here *if* BSS no longer creates such non-sense files.
        # this should be an option.

        if scaninfo.is_shutterless():
            r = scaninfo.get_file_number_based_on_template(filename)
            num = int(r.group(1))
            if scaninfo.hpoints > 1:
                return num%(scaninfo.hpoints+1) != 0 # if remainder is 0, discard the image.
            else:
                return num != 0 # discard 000.img
        else:
            return True
    # needs_to_be_processed()

    def get_grid_coord(self, filename):
        dirname = os.path.dirname(filename)
        if dirname not in self.scanlog:
            shikalog.warning("get_grid_coord(): directory is not found: %s" % dirname)
            return None

        return self.scanlog[dirname].get_grid_coord(os.path.basename(filename))
    # get_grid_coord()

    def get_scan_info(self, filename):
        dirname = os.path.dirname(filename)
        basename = os.path.basename(filename)

        if not dirname in self.scanlog:
            shikalog.warning("get_scan_info(): directory is not found: %s" % dirname)
            return None
        for scan in reversed(self.scanlog[dirname].scans):
            if scan.match_file_with_template(filename):
                return scan

        shikalog.warning("get_scan_info(): Not in scans: %s" % dirname)
        return None
    # get_scan_info()

    def get_gonio_xyz_phi(self, filename):
        dirname = os.path.dirname(filename)
        basename = os.path.basename(filename)

        if not dirname in self.scanlog:
            return None
        for scan in reversed(self.scanlog[dirname].scans):
            for f, c in scan.filename_coords:
                if basename == f:
                    if scan.is_shutterless():
                        return list(c[0]) + [scan.fixed_spindle]
                    else:
                        return list(c[0]) + [scan.osc_start]
        return None

    # get_gonio_xyz_phi()

# class DiffScanManager


class WatchScanlogThread:
    def __init__(self, queue, topdir, beamline=None, expdate="today"):
        self.queue = queue
        self.topdir = topdir
        self.interval = 5
        self.thread = None
        #self.latest_dir = None

        self.keep_going = True
        self.running = True

        self.beamline = beamline
        self.last_bsslog = None
        self.last_bsslog_line = 0

        self.expdate = None
        if expdate != "today": self.expdate = datetime.datetime.strptime(expdate, "%Y-%m-%d")

        self.thread = threading.Thread(None, self.run)
        self.thread.daemon = True

        self._cached_dirs = {}
        #self.thread.start()

    def start(self, interval=None):
        # Thread should be already started.
        # Just start to notify the latest directory.

        self.notify_latest_dir = True
        #wx.PostEvent(self.parent, EventLogWatcherStarted())

        if interval is not None:
            self.interval = interval

        # If accidentally stopped
        if not self.is_running():
            self.keep_going = True
            self.running = True
            self.thread = threading.Thread(None, self.run)
            self.thread.daemon = True
            self.thread.start()

    def stop(self):
        pass

    def is_running(self): return self.thread is not None and self.thread.is_alive()

    def find_in_directory(self, topdir):
        scanlogs = [] # (filename, date)

        for root, dirnames, filenames in os.walk(topdir):
            if "diffscan.log" in filenames:
                scanlog = os.path.join(root, "diffscan.log")
                scanlogs.append((scanlog, os.path.getmtime(scanlog)))

        return scanlogs
    # find_in_directory()

    def find_in_bsslog(self, topdir):
        def read_bsslog_line(l):
            if self.beamline in ("41xu","45xu") and "/ramdisk/" in l: l = l.replace("/ramdisk/","/isilon/users/", 1) # XXX should check if Pilatus or not!!
            if topdir not in l: return
            l = l[l.index(topdir):]
            if " " in l: l = l[:l.index(" ")]
            if ",.img" in l: l = l[:l.index(",.img")]
            return os.path.dirname(l)

        basedate = datetime.datetime.today() if self.expdate is None else self.expdate

        if self.last_bsslog is None:
            shikalog.debug("checking yesterday's bss log")
            self.last_bsslog = os.path.join(params.blconfig, "log",
                                            (basedate - datetime.timedelta(days=1)).strftime("bss_%Y%m%d.log"))

            if not os.path.isfile(self.last_bsslog):
                shikalog.info("Yesterday's log not found: %s"%self.last_bsslog)

        current_bsslog = os.path.join(params.blconfig, "log", basedate.strftime("bss_%Y%m%d.log"))

        if self.last_bsslog is not None and self.last_bsslog != current_bsslog and os.path.isfile(self.last_bsslog):
            shikalog.debug("reading last-log %s from %d" % (os.path.basename(self.last_bsslog), self.last_bsslog_line))
            for i, l in enumerate(open(self.last_bsslog)):
                if i <= self.last_bsslog_line: continue
                # read last log!
                found = read_bsslog_line(l)
                if found is not None: self._cached_dirs[found] = time.time()

            # reset for reading current log
            self.last_bsslog_line = 0

        if os.path.isfile(current_bsslog):
            shikalog.debug("reading curr-log %s from %d" % (os.path.basename(current_bsslog), self.last_bsslog_line))
            i = -1 # in case empty file
            for i, l in enumerate(open(current_bsslog)):
                if i <= self.last_bsslog_line: continue
                # read current log!
                found = read_bsslog_line(l)
                if found is not None: self._cached_dirs[found] = time.time()

            # set for next reading
            self.last_bsslog_line = i
        else:
            shikalog.info("bsslog not found: %s"%current_bsslog)
                
        self.last_bsslog = current_bsslog

        scanlogs = map(lambda x: os.path.join(x, "diffscan.log"), self._cached_dirs)
        uid = os.getuid()
        scanlogs = filter(lambda x: os.path.isfile(x) and os.stat(x).st_uid==uid, scanlogs)
        if params.only_check_in_last_hours is not None and params.only_check_in_last_hours > 0:
            now = time.time()
            last_seconds = params.only_check_in_last_hours*60*60
            scanlogs = filter(lambda x: (now-os.path.getmtime(x))<last_seconds, scanlogs)
        if scanlogs: shikalog.debug("found diffscan.log in bsslog: %s" % scanlogs)

        for k in self._cached_dirs.keys():
            # clear old cache
            if time.time() - self._cached_dirs[k] > 60*5: del self._cached_dirs[k]

        return map(lambda x: (x, os.path.getmtime(x)), scanlogs)
    # find_in_bsslog()

    def run_inner(self, method="bsslog"):
        assert method in ("bsslog", "os.walk")
        
        startt = time.time()

        if method == "bsslog":
            scanlogs = self.find_in_bsslog(self.topdir)
        else:
            scanlogs = self.find_in_directory(self.topdir)

        shikalog.debug("WatchScanlogThread.run_inner(method=%s) took %.3f sec for finding" % (method,
                                                                                            time.time()-startt))

        if len(scanlogs) > 0:
            scanlogs.sort(key=lambda x:x[1], reverse=True)
            for x in scanlogs: self.queue.put(x)
    # run_inner()

    def run(self):
        def mysleep():
            if self.interval < 1:
                time.sleep(self.interval)
            else:
                for i in xrange(int(self.interval/.5)):
                    if self.keep_going:
                        time.sleep(.5)
        # mysleep()

        shikalog.info("WatchScanlogThread loop STARTED")
        while self.keep_going:
            #shikalog.debug("in WatchScanlogThread timer")
            try:
                self.run_inner()
            except:
                shikalog.error("Error in WatchScanlogThread\n%s" % (traceback.format_exc()))
            mysleep()

        shikalog.info("WatchScanlogThread loop FINISHED")
        self.running = False
    # run()
# class WatchScanlogThread

class WatchDirThread:
    def __init__(self, queue, pushport):
        self.interval = 5
        self.thread = None
        self.queue = queue
        self.dirs = set()
        self.diffscan_manager = DiffScanManager()
        self.zmq_context = zmq.Context()
        self.ventilator_send = self.zmq_context.socket(zmq.PUSH)
        self.ventilator_send.bind("tcp://*:%d"%pushport)

    def start(self, interval=None):
        self.stop()

        self.keep_going = True
        self.running = True
        if interval is not None:
            self.interval = interval

        self.thread = threading.Thread(None, self.run)
        self.thread.daemon = True
        self.thread.start()

    def stop(self):
        if self.is_running():
            #shikalog.info("Stopping WatchDirThread.. Wait.")
            self.keep_going = False
            self.thread.join()
        else:
            pass
            #shikalog.info("WatchDirThread already stopped.")

    def is_running(self): return self.thread is not None and self.thread.is_alive()

    def run_inner(self):
        #shikalog.debug("in WatchDirThread timer")
        startt = time.time()

        while not self.queue.empty():
            scanlog, scanlogmtime = self.queue.get()
            self.diffscan_manager.add_scanlog(scanlog)


        if params.mode == "eiger_streaming":
            # Not sure this is needed, but clearly *not* needed when mode != "eiger_streaming"
            # because get_unprocessed_images() calls update_scanlogs()
            self.diffscan_manager.update_scanlogs()
        else:
            new_imgs = self.diffscan_manager.get_unprocessed_images(env=params.env)
            #print "new_imgs=", new_imgs
            for img, idx, scan in new_imgs:
                img_actual = img
                if params.env=="ppu": img_actual = re.sub("^/isilon/users/", "/ramdisk/", img)
                header = dict(file_prefix=scan.get_prefix()[:-1],
                              frame=idx-1,
                              raster_horizontal_number=scan.hpoints, raster_vertical_number=scan.vpoints,
                              raster_horizontal_step=scan.hstep, raster_vertical_step=scan.vstep,
                              raster_scan_direction=scan.scan_direction, raster_scan_path=scan.scan_path,
                              )
                shikalog.debug("Sending %s,%s" % (img, idx))
                msg = dict(imgfile=img, header=header,
                           cbf_data=open(img_actual, "rb").read())
                self.ventilator_send.send_pyobj(msg)

        shikalog.debug("WatchDirThread.run_inner took %.3f sec" % (time.time()-startt))

    def run(self):
        #shikalog.info("WatchDirThread loop STARTED")
        while self.keep_going:
            try:
                self.run_inner()
            except:
                shikalog.error("Error in WatchDirThread\n%s" % (traceback.format_exc()))

            if self.interval < 1:
                time.sleep(self.interval)
            else:
                for i in xrange(int(self.interval/.5)):
                    if self.keep_going:
                        time.sleep(.5)

        shikalog.info("WatchDirThread loop FINISHED")
        self.running = False
        #wx.PostEvent(self.parent, EventDirWatcherStopped()) # Ensure the checkbox unchecked when accidentally exited.

    # run()
# class WatchDirThread

def walk_nolink(top, topdown=True, onerror=None):
    # Original /misc/oys/xtal/cctbx/snapshots/dials-v1-8-3/base/lib/python2.7/os.py
    try:
        names = os.listdir(top)
    except os.error, err:
        if onerror is not None:
            onerror(err)
        return

    dirs, nondirs = [], []
    for name in names:
        try: st = os.lstat(os.path.join(top, name))
        except: continue # ignore deleted file

        if stat.S_ISDIR(st.st_mode):
            dirs.append((name,st))
        else:
            nondirs.append((name,st))

    if topdown:
        yield top, dirs, nondirs
    for name,_ in dirs:
        new_path = os.path.join(top, name)
        for x in walk_nolink(new_path, topdown, onerror):
            yield x
    if not topdown:
        yield top, dirs, nondirs
# walk_nolink()

class WatchRamdiskThread:
    def __init__(self, pushport, interval):
        self.interval = interval
        self.thread = None
        self.zmq_context = zmq.Context()
        self.ventilator_send = self.zmq_context.socket(zmq.PUSH)
        self.ventilator_send.bind("tcp://*:%d"%pushport)

    def start(self, interval=None):
        self.stop()

        self.keep_going = True
        self.running = True
        if interval is not None:
            self.interval = interval

        self.thread = threading.Thread(None, self.run)
        self.thread.daemon = True
        self.thread.start()

    def stop(self):
        if self.is_running():
            #shikalog.info("Stopping WatchRamdiskThread.. Wait.")
            self.keep_going = False
            self.thread.join()
        else:
            pass
            #shikalog.info("WatchRamdiskThread already stopped.")

    def is_running(self): return self.thread is not None and self.thread.is_alive()

    def run_inner_inotify(self):
        #shikalog.debug("in WatchRamdiskThread timer")
        startt = time.time()

        def is_processed(path, touch):
            if not os.path.exists(touch): return False
            mtime_touch = os.path.getmtime(touch)
            mtime_path = os.path.getmtime(path)
            return mtime_touch > mtime_path

        itree = inotify.adapters.InotifyTree('/ramdisk',
                                             mask=inotify.constants.IN_MOVED_TO|inotify.constants.IN_CLOSE_WRITE)
        for header, type_names, path, filename in itree.event_gen(yield_nones=False):
            print type_names, path, filename
            if "IN_ISDIR" in type_names: continue
            #if "IN_MOVED_TO" not in type_names: continue
            if not filename.endswith(".cbf"): continue
            imgf = os.path.join(path, filename)

            # Check if already processed
            touch = imgf+".touch"
            if is_processed(imgf, touch): continue

            # Send it!
            img_isilon = re.sub("^/ramdisk/", "/isilon/users/", imgf)
            header = dict(file_prefix=os.path.basename(img_isilon[:img_isilon.rindex("_")]),
                          frame=int(imgf[imgf.rindex("_")+1:imgf.rindex(".cbf")])-1)
            shikalog.debug("Sending %s,%s" % (header["file_prefix"], header["frame"]+1))
            msg = dict(imgfile=img_isilon, header=header,
                       cbf_data=open(imgf, "rb").read())
            self.ventilator_send.send_pyobj(msg)
            util.touch_file(touch) # mark as processed

        shikalog.debug("WatchRamdiskThread.run_inner stopped in %.3f sec" % (time.time()-startt))
    # run_inner_inotify()

    def run_inner_walk(self):
        def need_process(root, path_lst, filename_dict):
            # filenames is list of (filename, lstat); directory name not included
            # path_lst is (path, lstat(path)); directory name not included
            path = os.path.join(root, path_lst[0])
            lst = path_lst[1]
            # don't process if link
            if stat.S_ISLNK(lst.st_mode): return False
            touch = path_lst[0] + ".touch"
            # process if .touch not exists
            mtime_touch = filename_dict.get(touch, None)
            if mtime_touch is None: return lst.st_size >= 1000 and "Comment" in open(path).read(1000)
            
            mtime_path = lst.st_mtime
            if mtime_touch > mtime_path:
                return False
            else:
                shikalog.debug("Newer than touch (%s; %s <= %s)"%(path, mtime_touch, mtime_path))
                return lst.st_size >= 1000 and "Comment" in open(path).read(1000)
        # need_process()

        uid = os.getuid()
        start_time = time.time()
        n_dir = 0
        for root, dirnames, filenames in walk_nolink("/ramdisk"):
            n_dir += 1
            if os.stat(root).st_uid != uid: continue
            cbf_files = filter(lambda x: x[0].endswith(".cbf"), filenames)
            filename_dict = dict(filenames)
            new_files = filter(lambda x: need_process(root, x, filename_dict), cbf_files)
            new_files = map(lambda x: os.path.join(root, x[0]), new_files)

            for imgf in new_files:
                img_isilon = re.sub("^/ramdisk/", "/isilon/users/", imgf)
                header = dict(file_prefix=os.path.basename(img_isilon[:img_isilon.rindex("_")]),
                              frame=int(imgf[imgf.rindex("_")+1:imgf.rindex(".cbf")])-1)
                shikalog.debug("Sending %s" % img_isilon)
                msg = dict(imgfile=img_isilon, header=header,
                           cbf_data=open(imgf, "rb").read())
                self.ventilator_send.send_pyobj(msg)
                util.touch_file(imgf+".touch") # mark as processed
            
        shikalog.debug("WatchRamdiskThread.run_inner_walk finished in %.3f sec (%d dirs)" % (time.time()-start_time, n_dir))
    # run_inner_walk()

    def run(self):
        shikalog.info("WatchRamdiskThread loop STARTED")
        while self.keep_going:
            try:
                self.run_inner_walk()
            except:
                shikalog.error("Error in WatchRamdiskThread\n%s" % (traceback.format_exc()))

            if self.interval < 1:
                time.sleep(self.interval)
            else:
                for i in xrange(int(self.interval/.5)):
                    if self.keep_going:
                        time.sleep(.5)

        shikalog.info("WatchRamdiskThread loop FINISHED")
        self.running = False
    # run()
# class WatchRamdiskThread

class ResultsManager:
    def __init__(self, rqueue, dbdir):
        self.thread = threading.Thread(None, self.run)
        self.thread.daemon = True
        self.interval = 3

        self.dbdir = dbdir
        self.rqueue = rqueue
        self._diffscan_params = {}
    # __init__()

    def start(self):
        if not self.is_running():
            self.keep_going = True
            self.running = True
            self.thread.start()
    # start()

    def is_running(self): return self.thread is not None and self.thread.is_alive()

    def update_diffscan_params(self, msg):
        wdir = str(msg["data_directory"])
        prefix = str(msg["file_prefix"])

        key = os.path.join(wdir, "_spotfinder",  prefix)
        if key in self._diffscan_params: del self._diffscan_params[key]

        hstep = float(msg["raster_horizontal_step"])
        # raster_horizontal_number was raster_horizotal_number until bss_jul04_2017
        hpoint = int(msg.get("raster_horizotal_number", msg.get("raster_horizontal_number")))
        vstep = float(msg["raster_vertical_step"])
        vpoint = int(msg["raster_vertical_number"])


        if "raster_scan_direction" in msg and "raster_scan_path" in msg:
            scan_direction = str(msg["raster_scan_direction"])
            scan_path = str(msg["raster_scan_path"])
        else:
            scanlog = os.path.join(wdir, "diffscan.log")
            if not os.path.isfile(scanlog): return

            slog = bl_logfiles.BssDiffscanLog(scanlog)
            slog.remove_overwritten_scans()
            matched = filter(lambda x: x.get_prefix()==prefix+"_", slog.scans)
            if matched:
                scan_direction = matched[-1].scan_direction
                scan_path = matched[-1].scan_path
            else:
                return

        self._diffscan_params[key] = (vpoint, vstep, hpoint, hstep, scan_direction, scan_path)
        shikalog.info("_diffscan_params %s = %s" % (key, self._diffscan_params[key]))
    # update_diffscan_params()

    def get_raster_grid_coordinate(self, msg):
        try:
            header = msg["header"]
            hstep = float(header["raster_horizontal_step"])
            # raster_horizontal_number was raster_horizotal_number until bss_jul04_2017
            hpoint = int(header.get("raster_horizotal_number", header.get("raster_horizontal_number")))
            vstep = float(header["raster_vertical_step"])
            vpoint = int(header["raster_vertical_number"])
            scan_direction = str(header["raster_scan_direction"])
            scan_path = str(header["raster_scan_path"])
            return bl_logfiles.BssDiffscanLog.get_grid_coord_internal(vpoint, vstep, hpoint, hstep, msg["idx"], False, scan_direction, scan_path)
        except (KeyError, TypeError):
            shikalog.warning("not bringing info for grid coord: %s %s" % (msg.get("file_prefix"), msg.get("idx")))
            #pass

        wdir = str(msg["work_dir"])
        key = os.path.join(wdir, str(msg["file_prefix"]))
        if key in self._diffscan_params:
            vpoint, vstep, hpoint, hstep, scan_direction, scan_path = self._diffscan_params[key]
            return bl_logfiles.BssDiffscanLog.get_grid_coord_internal(vpoint, vstep, hpoint, hstep, msg["idx"], False, scan_direction, scan_path)
        else:
            shikalog.warning("_diffscan_params not available for %s" % key)
            scanlog = os.path.join(wdir, "..", "diffscan.log")
            gcxy = None
            if os.path.isfile(scanlog):
                slog = bl_logfiles.BssDiffscanLog(scanlog)
                slog.remove_overwritten_scans()
                gcxy = slog.calc_grid_coord(prefix=str(msg["file_prefix"])+"_", num=msg["idx"]) # may return None

            if gcxy is None: gcxy = [float("nan")]*2
            return gcxy
    # get_raster_grid_coordinate()


    def run(self):
        shikalog.info("ResultsManager loop STARTED")
        dbfile, summarydat, con, cur = None, None, None, None

        rcon = sqlite3.connect(os.path.join(self.dbdir, "%s.db"%getpass.getuser()), timeout=10)
        rcur = rcon.cursor()
        #c = rcur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='updates';")
        c = retry_until_success(rcur.execute, "SELECT name FROM sqlite_master WHERE type='table' AND name='updates';")
        if c.fetchone() is None:
            #rcur.execute("""create table updates (dirname text primary key,
            #                                      time real);""")
            retry_until_success(rcur.execute, """create table updates (dirname text primary key,
                                                  time real);""")

        while self.keep_going:
            try:
                messages = collections.OrderedDict()
                while not self.rqueue.empty():
                    msg = self.rqueue.get()
                    if "bss_job_mode" in msg:
                        shikalog.info("bssinfo= %s"%msg)
                        try:
                            self.update_diffscan_params(msg)
                        except:
                            shikalog.error("Error in update_diffscan_params%s" % (traceback.format_exc()))
                    else:
                        messages.setdefault(os.path.normpath(str(msg["work_dir"])), []).append(msg)

                for wdir in messages:
                    tmp = os.path.join(wdir, "shika.db")
                    if dbfile != tmp:
                        dbfile = tmp
                        con, cur = None, None
                        for _ in xrange(10):
                            try:
                                con = sqlite3.connect(dbfile, timeout=30)
                                break
                            except sqlite3.OperationalError:
                                shikalog.warning("Connecting to %s failed. Retrying" % dbfile)

                        if con is None:
                            shikalog.error("Could not connect to %s." % dbfile)
                        else:
                            cur = con.cursor()

                        summarydat = os.path.join(wdir, "summary.dat")
                        if not os.path.isfile(summarydat) or not os.path.getsize(summarydat):
                            open(summarydat, "w").write("prefix x y kind data filename\n")

                    canvas_data = {}
                            
                    for msg in messages[wdir]:
                        imgf = os.path.basename(str(msg["imgfile"]))
                        spots_is = map(lambda x: x[2], msg["spots"])

                        if cur is not None:
                            # 'status'
                            c = cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='status';")
                            if c.fetchone() is None:
                                cur.execute("""create table status (filename text primary key);""")
    
                            cur.execute("insert or replace into status values (?)", (imgf,)) # for backward 
    
                            # 'stats'
                            c = cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='stats';")
                            if c.fetchone() is None:
                                cur.execute("""create table stats (imgf text primary key, nspot real, total real, mean real);""")
    
                            cur.execute("insert or replace into stats values (?, ?,?,?)", 
                                        (imgf, len(msg["spots"]),
                                         sum(spots_is),
                                         sum(spots_is) / len(msg["spots"]) if len(msg["spots"])>0 else 0
                                         ))
    
                            # 'spots'
                            c = cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='spots';")
                            if c.fetchone() is None:
                                cur.execute("create table spots (filename text primary key, spots blob);")

                        # save jpg
                        if "jpgdata" in msg and msg["jpgdata"]:
                            jpgdir = os.path.join(wdir, 
                                                  "thumb_%s_%.3d" % (str(msg["file_prefix"]), msg["idx"]//1000))
                            try: os.mkdir(jpgdir)
                            except: pass
                            jpgout = os.path.join(jpgdir, imgf+".jpg")
                            open(jpgout, "wb").write(msg["jpgdata"])
                            del msg["jpgdata"]
                        elif "thumbdata" in msg and msg["thumbdata"]:
                            #assert len(msg["thumbdata"])==600*600*3
                            thumbw = int(numpy.sqrt(len(msg["thumbdata"])/3))
                            assert len(msg["thumbdata"]) == 3*thumbw*thumbw
                            prefix = str(msg["file_prefix"])
                            idx = (msg["idx"]-1)//100
                            tmpdir = os.path.join(os.path.expanduser("~"), ".shikatmp")
                            if not os.path.exists(tmpdir): os.mkdir(tmpdir)
                            tmpfile = os.path.join(tmpdir, "%s_%s_%.3d.pkl" % (hashlib.sha256(wdir).hexdigest(), prefix, idx))

                            if (prefix, idx) in canvas_data: 
                                canvas, ninc, _ = canvas_data[(prefix, idx)]
                            elif os.path.isfile(tmpfile):
                                shikalog.debug("loading thumbnail data from %s" % tmpfile)
                                canvas, ninc, _ = pickle.load(open(tmpfile))
                            else:
                                canvas, ninc = Image.new("RGB", (thumbw*10, thumbw*10), (0, 0, 0)), 0
                                
                            thumbimg = Image.frombytes("RGB", (thumbw, thumbw), msg["thumbdata"])
                            idx2 = (msg["idx"]-1)%100
                            x, y = idx2%10, idx2//10
                            canvas.paste(thumbimg, (x*thumbw, y*thumbw))
                            ninc += 1
                            #canvas.save(jpgout, "JPEG", quality=90, optimize=True) # IT MAY COST TOO MUCH TIME (SHOULD RUN THIS JUST ONCE?)

                            # calc max frame number
                            hpoint = int(msg["header"].get("raster_horizotal_number", msg["header"].get("raster_horizontal_number"))) # raster_horizontal_number was raster_horizotal_number until bss_jul04_2017
                            vpoint = int(msg["header"]["raster_vertical_number"])
                            n_max = hpoint * vpoint
                            idx_max = (n_max-1)//100
                            n_in_last = n_max - idx_max*100
                            
                            canvas_data[(prefix, idx)] = (canvas, ninc,
                                                          ninc == 100 or (idx_max==idx and ninc == n_in_last) # jpeg-completed flag
                                                          )
                            
                            del msg["thumbdata"]
                            

                        if cur is not None:
                            cur.execute("insert or replace into spots values (?, ?)",
                                        (imgf, sqlite3.Binary(pickle.dumps(msg, -1))))

                        # summary.dat
                        try:
                            #tmptmp = time.time()
                            gcxy = self.get_raster_grid_coordinate(msg)

                            with open(summarydat, "a") as ofs:
                                kinds = ("n_spots", "total_integrated_signal","median_integrated_signal")
                                data = (len(msg["spots"]), sum(spots_is), numpy.median(spots_is))
                                for k, d in zip(kinds, data):
                                    ofs.write("%s_ % .4f % .4f %s %s %s\n" % (str(msg["file_prefix"]),
                                                                              gcxy[0], gcxy[1], k, d, imgf))
                            #shikalog.info("##DEBUG time: %f" %( time.time()-tmptmp))
                        except:
                            shikalog.error("Error in summary.dat generation at %s\n%s" % (wdir, traceback.format_exc()))

                            
                    for prefix, idx in canvas_data:
                        tmpdir = os.path.join(os.path.expanduser("~"), ".shikatmp")
                        jpgdir = os.path.join(wdir, "thumb_%s" % prefix)
                        if not os.path.exists(tmpdir): os.mkdir(tmpdir)
                        if not os.path.exists(jpgdir): os.mkdir(jpgdir)
                        tmpfile = os.path.join(tmpdir, "%s_%s_%.3d.pkl" % (hashlib.sha256(wdir).hexdigest(), prefix, idx))                        
                        jpgout = os.path.join(jpgdir, "%s_%.6d-%.6d.jpg" % (prefix, idx*100+1, (idx+1)*100))
                        jpgtmp = os.path.join(jpgdir, ".tmp-%s_%.6d-%.6d.jpg" % (prefix, idx*100+1, (idx+1)*100))
                        
                        shikalog.info("saving thumbnail jpeg as %s" % jpgout)
                        
                        canvas_data[(prefix, idx)][0].save(jpgtmp, "JPEG", quality=50, optimize=True)
                        os.rename(jpgtmp, jpgout) # as it may take time
                        
                        if canvas_data[(prefix, idx)][2]:
                            if os.path.isfile(tmpfile): os.remove(tmpfile)
                        else:
                            shikalog.info("saving thumbnail data to %s" % tmpfile)
                            pickle.dump(canvas_data[(prefix, idx)], open(tmpfile, "w"), -1)
                            
                    while True:
                        try: con.commit()
                        except sqlite3.OperationalError:
                            shikalog.warning("sqlite3.OperationalError. Retrying.")
                            time.sleep(1)
                            continue
                        break

                    rcur.execute("insert or replace into updates values (?,?)", (wdir, time.time()))
                    rcon.commit()
                    shikalog.info("%4d results updated in %s" % (len(messages[wdir]), wdir))
            except:
                shikalog.error("Exception: %s" % traceback.format_exc())
            
            time.sleep(self.interval)

        self.running = False
        shikalog.info("ResultsManager loop FINISHED")
    # run()
# ResultsManager

def results_receiver(rqueue, pullport, results_manager):
    zmq_context = zmq.Context()
    receiver = zmq_context.socket(zmq.PULL)
    receiver.bind("tcp://*:%d"%pullport)

    last_times = []

    while True:
        msg = receiver.recv_pyobj()

        if "bss_job_mode" not in msg:
            last_times.append(time.time())
            if len(last_times) > 50: last_times = last_times[-50:]
            last_times = filter(lambda x: last_times[-1]-x < 50, last_times)
            hz = float(len(last_times))/(last_times[-1]-last_times[0]) if len(last_times)>1 else 0
            shikalog.info("%s %6.2f Hz (last %2d)" % (msg["imgfile"], hz, len(last_times)))

        rqueue.put(msg)
        #results_manager.start()
# results_receiver()

def worker(wrk_num, params):
    context = zmq.Context()
 
    # Set up a channel to receive work from the ventilator
    work_receiver = context.socket(zmq.PULL)
    work_receiver.connect("tcp://127.0.0.1:%d"%params.ports[0])

    eiger_receiver = context.socket(zmq.PULL)
    eiger_receiver.connect("tcp://%s:9999"%params.eiger_host)
 
    # Set up a channel to send result of work to the results reporter
    results_sender = context.socket(zmq.PUSH)
    results_sender.connect("tcp://127.0.0.1:%d"%params.ports[1])
 
    # Set up a channel to receive control messages over
    control_receiver = context.socket(zmq.SUB)
    control_receiver.connect("tcp://127.0.0.1:%d"%params.ports[2])
    control_receiver.setsockopt(zmq.SUBSCRIBE, "")
 
    # Set up a poller to multiplex the work receiver and control receiver channels
    poller = zmq.Poller()
    poller.register(work_receiver, zmq.POLLIN)
    poller.register(control_receiver, zmq.POLLIN)
    if params.mode == "eiger_streaming":
        poller.register(eiger_receiver, zmq.POLLIN)

    params_dict = {}
    for key in config_manager.sp_params_strs:
        params_str = config_manager.sp_params_strs[key] + config_manager.get_common_params_str()
        master_params = libtbx.phil.parse(spot_finder_for_grid_scan.master_params_str)
        working_params = master_params.fetch(sources=[libtbx.phil.parse(params_str)])
        params_dict[key] = working_params.extract()

    shikalog.info("worker %d ready" % wrk_num)

    # Loop and accept messages from both channels, acting accordingly
    while True:
        socks = dict(poller.poll())
 
        # the message came from work_receiver channel
        if socks.get(work_receiver) == zmq.POLLIN:
            msg = work_receiver.recv_json()
            imgfile = str(msg["imgfile"])
            pkey = config_manager.get_key_by_img(imgfile)
            #params_str = config_manager.sp_params_strs[pkey] + config_manager.get_common_params_str()
            
            shikalog.info("wrker%.2d: %s"%(wrk_num, msg))

            #master_params = libtbx.phil.parse(spot_finder_for_grid_scan.master_params_str)
            #working_params = master_params.fetch(sources=[libtbx.phil.parse(params_str)])
            #working_params.show()
            #dparams = working_params.extract()
            dparams = params_dict[pkey]
            dparams.work_dir = os.path.join(os.path.dirname(imgfile), "_spotfinder") # PPU-case!!?
            if os.path.exists(dparams.work_dir): assert os.path.isdir(dparams.work_dir)
            else:
                try: os.mkdir(dparams.work_dir)
                except: pass

            result = spot_finder_for_grid_scan.run(imgfile, dparams)
            result.update(msg)
            result["work_dir"] = dparams.work_dir
            result["params"] = dparams
            results_sender.send_pyobj(result)

        # the message from EIGER
        if socks.get(eiger_receiver) == zmq.POLLIN:
            frames = eiger_receiver.recv_multipart(copy = False)

            header, data = eiger.read_stream_data(frames)
            if util.None_in(header, data): continue

            #params_str = config_manager.sp_params_strs[("BL32XU", "EIGER9M", None, None)] + config_manager.get_common_params_str()
            #master_params = libtbx.phil.parse(spot_finder_for_grid_scan.master_params_str)
            #working_params = master_params.fetch(sources=[libtbx.phil.parse(params_str)])
            #working_params.show()
            dparams = params_dict[("BL32XU", "EIGER9M", None, None)] #working_params.extract()
            dparams.work_dir = os.path.join(str(header["data_directory"]), "_spotfinder")
            if os.path.exists(dparams.work_dir): assert os.path.isdir(dparams.work_dir)
            else:
                try: os.mkdir(dparams.work_dir)
                except: pass

            shikalog.info("Got data: %s"%header)

            imgfile = os.path.join(header["data_directory"],
                                   "%s_%.6d.img"%(str(header["file_prefix"]), header["frame"]+1))

            result = spot_finder_for_grid_scan.run(imgfile, dparams, data_and_header=(data, header))
            result["work_dir"] = dparams.work_dir
            result["params"] = dparams
            result["imgfile"] = imgfile
            result["template"] = "%s_%s.img"%(str(header["file_prefix"]), "?"*6)
            result["idx"] = header["frame"]+1
            results_sender.send_pyobj(result)

            #os.remove(imgfile)
            

        # the message for control
        if socks.get(control_receiver) == zmq.POLLIN:
            msg = control_receiver.recv_pyobj()
            if "params" in msg:
                params_dict = msg["params"]
                shikalog.info("worker %d: Parameters updated" % wrk_num)

    shikalog.info("Worker %d finished." % wrk_num)
# worker()

def run_from_args(argv):
    if "-h" in argv or "--help" in argv:
        print "All parameters:\n"
        iotbx.phil.parse(master_params_str).show(prefix="  ", attributes_level=1)
        return

    cmdline = iotbx.phil.process_command_line(args=argv,
                                              master_string=master_params_str)
    global params
    params = cmdline.work.extract()
    args = cmdline.remaining_args

    if params.topdir is None: params.topdir = os.getcwd()

    if params.mode == "zoo":
        if params.bl == "32xu": params.blconfig = "/isilon/BL32XU/BLsoft/PPPP/10.Zoo/ZooConfig/"
        elif params.bl == "26b2": params.blconfig = "/isilon/users/rikenbl/rikenbl/Zoo/ZooConfig/"

    if params.blconfig is None:
        params.blconfig = "/isilon/blconfig/bl%s" % params.bl

    shikalog.config(params.bl, "backend", params.logroot)

    if params.force_ssh_from is not None:
        shikalog.info("SSH_CONNECTION= %s" % os.environ.get("SSH_CONNECTION", ""))
        if "SSH_CONNECTION" not in os.environ:
            print
            print "ERROR!! Cannot get host information! SHIKA cannot start."
            print "Please contact staff. Need to use ssh or change parameter."
            print
            return

        ssh_from = os.environ["SSH_CONNECTION"].split()[0]
        if ssh_from != params.force_ssh_from:
            print
            print "ERROR!! Your host is not allowed! SHIKA cannot start here."
            if ssh_from[:ssh_from.rindex(".")] == "10.10.126" and 162 <= int(ssh_from[ssh_from.rindex(".")+1:]) <= 170:
                shikalog.debug("Access from Oyster computer.")
                print "You appear to be in Oyster. Please *LOGOUT* from Oyster, and try to start SHIKA from your local computer!!"
            else:
                print "Please contact staff. Need to access from the allowed host."
            print
            return

    #import subprocess
    #import sys
    #pickle.dump(params, open("/tmp/params.pkl","w"),-1)
    #pp = []
    for i in xrange(params.nproc):
        Process(target=worker, args=(i,params)).start()
        #p = subprocess.Popen(["%s -"%sys.executable], shell=True, stdin=subprocess.PIPE)
        #p.stdin.write("from yamtbx.dataproc.myspotfinder.command_line.spot_finder_backend import worker\nimport pickle\nworker(%d, pickle.load(open('/tmp/params.pkl')))"%i)
        #p.stdin.close()
        #pp.append(p)

    rqueue = Queue.Queue()
    results_manager = ResultsManager(rqueue=rqueue, dbdir=params.dbdir)

    if params.mode == "watch_ramdisk":
        ramdisk_watcher = WatchRamdiskThread(pushport=params.ports[0],
                                             interval=params.ramdisk_walk_interval)
        ramdisk_watcher.start()
    elif params.mode != "eiger_streaming": 
        queue = Queue.Queue()
        scanlog_watcher = WatchScanlogThread(queue, topdir=params.topdir,
                                             beamline=params.bl, expdate=params.date)
        dir_watcher = WatchDirThread(queue, pushport=params.ports[0])
        scanlog_watcher.start()
        dir_watcher.start()

    results_manager.start() # this blocks!??!

    results_receiver(rqueue=rqueue, pullport=params.ports[1], results_manager=results_manager)

    #for p in pp: p.wait()

    if params.nproc==0:
        while True:
            time.sleep(1)

# run_from_args()

if __name__ == "__main__":
    import sys
    run_from_args(sys.argv[1:])


