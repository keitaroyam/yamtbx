"""
(c) RIKEN 2017. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals

import wx
import wx.lib.newevent
import wx.lib.agw.pybusyinfo
import os
from yamtbx import util
from yamtbx.dataproc import bl_logfiles
from libtbx import adopt_init_args
from libtbx.utils import multi_out
import threading
import time
import traceback
import shutil
import fnmatch
import sys
import re
import glob

re_data_h5 = re.compile(r"^(.+)_data_([0-9]+).h5$")
re_eiger_h5 = re.compile(r"^(.+)_(master|data_([0-9]+)).h5$")
logfile_name = None #"bl32xu_datasync.log"

EventSyncError, EVT_SYNC_ERROR = wx.lib.newevent.NewEvent()
EventSyncDone, EVT_SYNC_DONE = wx.lib.newevent.NewEvent()

def get_mounted_dirs():
    # XXX Linux specific
    ret, out, err = util.call("df")
    dirs = [x.split()[-1] for x in out.splitlines()[1:]]
    dirs = [x for x in dirs if x.startswith("/media/")]
    return dirs
# get_mounted_dirs()

class SyncOptions(object):
    def __init__(self, sdir=None, ddir=None, exclude_dirs=[], exclude_files=[],
                 copy_shika_results=True, copy_scan_data=True, copy_hits=True):
        adopt_init_args(self, locals())
    # __init__()

    def show(self, out):
        out.write("""\
# Options:
#  source_dir = %(sdir)s
#  target_dir = %(ddir)s
#  copy_shika_results = %(copy_shika_results)s
#  copy_scan_files = %(copy_scan_data)s
#  copy_scan_hit_files = %(copy_hits)s
#  exclude_dirs = %(exclude_dirs)s
#  exclude_files = %(exclude_files)s
""" % self.__dict__)
    # show()

# class SyncOptions

class DataSyncThread(object):
    def __init__(self, parent):
        self.parent = parent
        self.interval = 10
        self.thread = None

        self.opts = SyncOptions()
        self.log_out = sys.stdout
        self._transferred = {}
        self._bsslog_cache = {}
        self.skipped = set()
    # __init__()

    def start(self, interval=None, opts=None):
        self.stop()

        if self.opts.sdir!=opts.sdir or self.opts.ddir!=opts.ddir:
            self._transferred = {}
            self._bsslog_cache = {}
            self.skipped = set()

        self.opts = opts

        self.log_out = multi_out()
        self.log_out.register("stdout", sys.stdout)
        self.log_out.register("file", open(os.path.join(self.opts.sdir, logfile_name), "a"))

        self.log_out.write(time.strftime("# DATA SYNC GUI - Starting at %Y-%m-%d %H:%M:%S\n"))
        self.opts.show(self.log_out)
        self.log_out.write("\n")
        self.log_out.flush()

        self.keep_going = True
        self.stop_now = False # for emergency stop
        self.running = True
        if interval is not None:
            self.interval = interval

        self.thread = threading.Thread(None, self.run)
        self.thread.daemon = True
        self.thread.start()
    # start()

    def stop(self, stop_now=False):
        if self.is_running():
            self.stop_now = stop_now
            self.keep_going = False
            self.thread.join()
    # stop()

    def is_running(self): return self.thread is not None and self.thread.is_alive()

    def run(self):
        counter = 0
        while self.keep_going:
            counter += 1
            if self.interval == 0: # Run only once
                self.keep_going = False
                continue

            try:
                self.log_out.write(time.strftime("# %Y-%m-%d %H:%M:%S syncing..\n"))
                self.log_out.flush()
                st = time.time()
                bytes_total, files_copied, skipped = self.do_sync()
                eltime = time.time() - st
                self.log_out.write(time.strftime("# %Y-%m-%d %H:%M:%S sync done. "))
                self.log_out.write("Total %.1f %s copied " % util.human_readable_bytes(bytes_total))
                self.log_out.write("in %.1f sec. " % eltime)
                self.log_out.write("(%.1f %s/sec.)\n\n" % util.human_readable_bytes(float(bytes_total)/eltime))
                self.log_out.flush()

                ev = EventSyncDone(bytes_total=bytes_total,
                                   files_copied=files_copied,
                                   skipped=skipped)
                wx.PostEvent(self.parent, ev)
            except:
                self.log_out.write("Error occurred in syncing.\n")
                self.log_out.write(traceback.format_exc()+"\n")
                self.log_out.flush()
                self.keep_going = False
                ev = EventSyncError(msg=traceback.format_exc())
                wx.PostEvent(self.parent, ev)

            if self.interval < 1:
                time.sleep(self.interval)
            else:
                for i in range(int(self.interval/.5)):
                    if self.keep_going:
                        time.sleep(.5)

        self.running = False
        self.log_out.write(time.strftime("# %Y-%m-%d %H:%M:%S Stopped"))
        self.log_out.write(" urgently.\n" if self.stop_now else ".\n")
        self.log_out.flush()
    # run()

    def get_total_bytes_sofar(self):
        ret = 0.
        for k in self._transferred: ret += self._transferred[k][1]
        return ret
    # get_total_bytes_sofar()

    def analyze_bsslogs_in_dir(self, wdir, files):
        """
        Analyze all BSS log files in `wdir', and return {prefix:"data"or"scan", ...}
        """
        if not files:
            return {}

        #self.log_out.write("debug:: checking logs in %s\n" % wdir)
        #self.log_out.flush()

        logfiles = [x for x in [os.path.join(wdir, x) for x in files] if x.endswith(".log")]
        logfiles_id = [(os.path.basename(x), os.path.getmtime(x), os.path.getsize(x)) for x in sorted(logfiles)]

        if wdir in self._bsslog_cache and self._bsslog_cache[wdir][0] == logfiles_id:
            #print "Using cached", self._bsslog_cache[wdir][1]
            return self._bsslog_cache[wdir][1]

        ret = {}
        tmp = {} # {prefix:[("data"|"scan", datetime.datetime object)]

        # Summarize valid data collections, and resolve overwritten issues
        for logf in logfiles:
            if os.path.basename(logf) == "diffscan.log":
                try:
                    slog = bl_logfiles.BssDiffscanLog(logf)
                    slog.remove_overwritten_scans()
                    for scan in slog.scans:
                        prefix = scan.get_prefix()[:-1] # ends with "_"
                        tmp.setdefault(prefix, []).append(("scan", scan.date.toordinal()))
                except:
                    print(traceback.format_exc())
                    continue
            else:
                try:
                    jlog = bl_logfiles.BssJobLog(logf)
                    jlog.annotate_overwritten_images(remove=True)
                    jobs = [x for x in jlog.jobs if x.job_mode.lower() not in ("xafs", "raster scan")]
                    for job in jlog.jobs:
                        t = job.scans[0].datetime.toordinal() if job.scans else -1
                        # job.prefix does not include "_"
                        tmp.setdefault(job.prefix, []).append(("data", t))
                except:
                    print(traceback.format_exc())
                    continue

        for p in tmp:
            newest = max(tmp[p], key=lambda x:x[1]) # take latest one
            ret[p] = newest[0]

        self._bsslog_cache[wdir] = (logfiles_id, ret)
        return ret
    # analyze_bsslogs_in_dir()

    def skip_directory(self, dname_rel): # XXX Not tested!
        """
        `dname' must be a relative path from source directory.
        assuming the path is already normalized.. true??
        """
        
        if not self.opts.exclude_dirs and self.opts.copy_shika_results:
            return False

        dirs = dname_rel.split(os.sep) # what if there is directory named "hoge/fuga"??

        if self.opts.exclude_dirs:
            for dname in dirs:
                if any([fnmatch.fnmatch(dname, x) for x in self.opts.exclude_dirs]):
                    return True

        if not self.opts.copy_shika_results and "_spotfinder" in dirs:
            return True

        return False
    # skip_directory()

    def skip_file(self, f, prefix_dict): # XXX Not tested!
        """
        `f' must be abspath.
        """
        if self.opts.exclude_files:
            fbase = os.path.basename(f)
            if any([fnmatch.fnmatch(fbase, x) for x in self.opts.exclude_files]):
                return True

        # If non-h5 file, do not skip (unless matched with user-specified pattern)
        if not f.endswith(".h5"): return False

        if not self.opts.copy_scan_data and not f.endswith("_onlyhits.h5"):
            #if not prefix_dict: return False

            r_h5_prefix = re_eiger_h5.search(os.path.basename(f))
            if r_h5_prefix:
                prefix = r_h5_prefix.group(1)
                #self.log_out.write("debug: %s %s\n" % (prefix, prefix_dict[prefix]))
                if prefix_dict.get(prefix, "") == "scan": return True
            """
            # If there is no diffscan.log in the same directory, the file never scan.
            # ..but, OVERWRITTEN_* directory has no log files.
            scanlog = os.path.join(os.path.dirname(f), "diffscan.log")
            if not os.path.exists(scanlog): return False

            fbase = os.path.basename(f)
            dirname = os.path.dirname(f)

            # If threre is <sameprefix>_onlyhits.h5, the related files must be scan files.
            r_h5_prefix = re_eiger_h5.search(fbase)
            if r_h5_prefix:
                prefix = r_h5_prefix.group(1)
                if os.path.exists(os.path.join(dirname, prefix+"_onlyhits.h5")):
                    return False

            # How can we check if this is scan file??
            # What if no onlyhits.h5.. maybe users want data if only few frames?
            if 0: return True # XXX check this!
            """

        # What if hit-extract failed??
        if not self.opts.copy_hits and f.endswith("_onlyhits.h5"):
            return True
        
        return False
    # skip_file()

    def need_copy(self, f_src, f_dst):
        time_tol = 2. # for FAT

        if not os.path.exists(f_dst):
            return True

        mtime_size = os.path.getmtime(f_src), os.path.getsize(f_src)
        mtime_size_saved = os.path.getmtime(f_dst), os.path.getsize(f_dst)

        print("debug::", f_dst, mtime_size[1] - mtime_size_saved[1], mtime_size[0]-mtime_size_saved[0])
        if mtime_size[1] == mtime_size_saved[1] and abs(mtime_size[0]-mtime_size_saved[0]) < time_tol:
            return False

        return True
    # need_copy()

    def do_sync(self):
        # TODO exit loop when requested
        sdir = self.opts.sdir
        ddir = self.opts.ddir
        log_out = self.log_out

        bytes_total, files_copied, skipped = 0, [], []

        for root, dirnames, filenames in os.walk(sdir, followlinks=True):
            root_rel = os.path.relpath(root, sdir)
            root_d = os.path.join(ddir, root_rel)

            if self.stop_now: return bytes_total, files_copied, skipped

            if self.skip_directory(os.path.join(sdir, root_rel)):
                if not root_rel in self.skipped:
                    log_out.write("skipping %s/\n"%root_rel)
                    self.skipped.add(root_rel)
                continue

            if not os.path.exists(root_d):
                log_out.write("mkdir %s\n"%root_d)
                os.makedirs(root_d)

            prefix_dict = self.analyze_bsslogs_in_dir(root, filenames)

            for f in filenames:
                if self.stop_now: return bytes_total, files_copied, skipped

                f_abs = os.path.join(root, f)
                f_rel = os.path.join(root_rel, f)
                if not os.path.exists(f_abs): continue # seems sometimes to happen..
                mtime, size = os.path.getmtime(f_abs), os.path.getsize(f_abs)
                if self.skip_file(f_abs, prefix_dict):
                    if f_rel not in self.skipped:
                        log_out.write("skipping %s\n"%f_rel)
                        self.skipped.add(f_rel)
                    skipped.append(f_rel)
                    continue

                if self.need_copy(f_abs, os.path.join(root_d, f)): #(mtime, size)):
                    log_out.write("copying %s "%f_rel)
                    st = time.time()
                    try:
                        util.safe_copy(f_abs, root_d)
                    except:
                        log_out.write(" failed.\n")
                        log_out.flush()
                        continue
                    self._transferred[f_abs] = (mtime, size)
                    files_copied.append(f_rel)
                    bytes_total += size
                    log_out.write(" %.1f %s/s\n" % util.human_readable_bytes(float(size)/(time.time()-st)))

            try:
                # we experienced OSError here once.. I think this can be ignored.
                shutil.copystat(root, root_d)
            except:
                log_out.write("Error in copying stat of %s\n" % root)
                log_out.write(traceback.format_exc()+"\n")

                
        return bytes_total, files_copied, skipped
    # do_sync()

# class DataSyncThread

class MainFrame(wx.Frame):
    def __init__(self, opts, parent=None, id=wx.ID_ANY):
        wx.Frame.__init__(self, parent=parent, id=id, title="User data backup for %s" % opts.bl,
                          size=(700,500))

        self.data_sync_thread = DataSyncThread(self)
        self.update_timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.on_update_timer, self.update_timer)
        self.update_timer.Start(5000)

        self.log_update_timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.on_log_update_timer, self.log_update_timer)
        self.log_update_timer.Start(500)
        self._logfile_pos = {}

        panel = wx.Panel(self, wx.ID_ANY)

        vbox = wx.BoxSizer(wx.VERTICAL)

        self.txtSourceDir = wx.TextCtrl(panel, wx.ID_ANY, size=(300, 25))
        self.txtSourceDir.SetValue(os.getcwd())
        self.txtDestDir = wx.TextCtrl(panel, wx.ID_ANY, style=wx.TE_PROCESS_ENTER, size=(300, 25))

        mounted_dirs = get_mounted_dirs()
        if mounted_dirs: self.txtDestDir.SetValue(mounted_dirs[0])
        
        self.txtSourceDir.Bind(wx.EVT_TEXT_ENTER, self.txtSource_enter)
        self.txtDestDir.Bind(wx.EVT_TEXT_ENTER, self.txtDestDir_enter)

        self.btnSourceDir = wx.Button(panel, wx.ID_ANY, "...", size=(25,25))
        self.btnDestDir = wx.Button(panel, wx.ID_ANY, "...", size=(25,25))
        self.btnUnmount = wx.Button(panel, wx.ID_ANY, "umount", size=(70,25))
        self.btnSourceDir.Bind(wx.EVT_BUTTON, self.btnSourceDir_click)
        self.btnDestDir.Bind(wx.EVT_BUTTON, self.btnDestDir_click)
        self.btnUnmount.Bind(wx.EVT_BUTTON, self.btnUnmount_click)

        grid1 = wx.GridBagSizer(3, 5)
        grid1.Add(wx.StaticText(panel, wx.ID_ANY, "Data location: "),
                  pos=(0,0), flag=wx.LEFT|wx.ALIGN_CENTER_VERTICAL, border=4)
        grid1.Add(self.txtSourceDir, pos=(0,1), flag=wx.EXPAND)
        grid1.Add(self.btnSourceDir, pos=(0,2), flag=wx.EXPAND)
        self.lblCopiedSize = wx.StaticText(panel, wx.ID_ANY, " (Total transferred size: ??? GB)")
        grid1.Add(self.lblCopiedSize, pos=(0,3), span=(1,2), flag=wx.LEFT|wx.ALIGN_CENTER_VERTICAL, border=4)
        grid1.Add(wx.StaticText(panel, wx.ID_ANY, "Destination: "),
                  pos=(1,0), flag=wx.LEFT|wx.ALIGN_CENTER_VERTICAL, border=4)
        grid1.Add(self.txtDestDir, pos=(1,1), flag=wx.EXPAND)
        grid1.Add(self.btnDestDir, pos=(1,2))
        grid1.Add(self.btnUnmount, pos=(1,3))
        self.lblDestLeft = wx.StaticText(panel, wx.ID_ANY, " (Space left: ??? GB)")
        grid1.Add(self.lblDestLeft, pos=(1,4), flag=wx.LEFT|wx.ALIGN_CENTER_VERTICAL, border=4)
        grid1.AddGrowableCol(1)
        vbox.Add(grid1)
        vbox.AddSpacer(5)
        
        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        vbox.Add(hbox1)

        self.cbCopyShikaResults = wx.CheckBox(panel, label="Copy SHIKA results")
        self.cbCopyShikaResults.SetValue(1)
        self.cbCopyScanOrg = wx.CheckBox(panel, label="Copy scan data (all images)")
        self.cbCopyScanHit = wx.CheckBox(panel, label="Copy scan data (hits only)")
        self.cbCopyScanHit.SetValue(1)
        sb1 = wx.StaticBox(panel, label="Copy options")
        svb1 = wx.StaticBoxSizer(sb1, wx.VERTICAL)
        svb1.Add(self.cbCopyShikaResults, flag=wx.LEFT|wx.TOP, border=5)
        svb1.Add(self.cbCopyScanOrg, flag=wx.LEFT|wx.TOP, border=5)
        svb1.Add(self.cbCopyScanHit, flag=wx.LEFT|wx.TOP, border=5)

        self.txtExcDirs = wx.TextCtrl(panel, wx.ID_ANY, style= wx.TE_MULTILINE, size=(160,80))
        self.txtExcDirs.SetValue("coot-*\n")
        self.txtExcFiles = wx.TextCtrl(panel, wx.ID_ANY, style= wx.TE_MULTILINE, size=(160,80))
        self.txtExcFiles.SetValue("0-coot-*\n.*\n")
        sb2 = wx.StaticBox(panel, label="Exclude options")
        svb2 = wx.StaticBoxSizer(sb2, wx.VERTICAL)
        grid2 = wx.GridBagSizer(2, 2)
        grid2.Add(wx.StaticText(panel, wx.ID_ANY, "directories: "),
                  pos=(0,0), flag=wx.LEFT|wx.ALIGN_CENTER_VERTICAL, border=4)
        grid2.Add(wx.StaticText(panel, wx.ID_ANY, "files: "),
                  pos=(0,1), flag=wx.LEFT|wx.ALIGN_CENTER_VERTICAL, border=4)
        grid2.Add(self.txtExcDirs, pos=(1,0), flag=wx.EXPAND, border=4)
        grid2.Add(self.txtExcFiles, pos=(1,1), flag=wx.EXPAND, border=4)
        svb2.Add(grid2)
        hbox1.Add(svb1, 1, flag=wx.EXPAND, border=4)
        hbox1.Add(svb2, 1, flag=wx.EXPAND, border=4)

        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        vbox.Add(hbox2)

        self.btnStart = wx.Button(panel, wx.ID_ANY, "Start", size=(200,50))
        self.btnStop = wx.Button(panel, wx.ID_ANY, "Stop after this round", size=(200,50))
        self.btnStopNow = wx.Button(panel, wx.ID_ANY, "Stop immediately", size=(200,50))
        self.btnStart.Bind(wx.EVT_BUTTON, self.btnStart_click)
        self.btnStop.Bind(wx.EVT_BUTTON, lambda e: self.btnStop_click(e, False))
        self.btnStopNow.Bind(wx.EVT_BUTTON, lambda e: self.btnStop_click(e, True))
        

        hbox2.AddStretchSpacer()
        #hbox2.AddSpacer(100)
        hbox2.Add(self.btnStart)
        #hbox2.AddSpacer(100)
        hbox2.Add(self.btnStop)
        hbox2.Add(self.btnStopNow)
        #hbox2.AddSpacer(100)
        #hbox2.AddStretchSpacer(prop=1)

        vbox.AddSpacer(5)

        notebook = wx.Notebook(panel, wx.ID_ANY)
        page1 = wx.Panel(notebook, wx.ID_ANY)
        page2 = wx.Panel(notebook, wx.ID_ANY)

        vboxp1 = wx.BoxSizer(wx.VERTICAL)
        page1.SetSizer(vboxp1)
        self.txtLog = wx.TextCtrl(page1, wx.ID_ANY, style= wx.TE_MULTILINE)
        self.txtLog.SetEditable(False)
        vboxp1.Add(self.txtLog, 1, wx.EXPAND)

        vboxp2 = wx.BoxSizer(wx.VERTICAL)
        page2.SetSizer(vboxp2)
        self.lstSummary = wx.ListCtrl(page2, wx.ID_ANY, style=wx.LC_REPORT)
        #self.lstSummary.SetEditable(False)
        for i,(v,s) in enumerate((("Date", 150), ("Transferred", 100), ("Files updated", 100), ("Skipped", 100))):
            self.lstSummary.InsertColumn(i,v)
            self.lstSummary.SetColumnWidth(i, s)

        vboxp2.Add(self.lstSummary, 1, wx.EXPAND)

        notebook.InsertPage(0, page1, "LogView")
        notebook.InsertPage(1, page2, "Summary")
        vbox.Add(notebook, 1, flag=wx.EXPAND)

        panel.SetSizer(vbox)

        self.Bind(wx.EVT_CLOSE, self.onClose)
        self.Bind(EVT_SYNC_ERROR, self.on_sync_error)
        self.Bind(EVT_SYNC_DONE, self.on_sync_done)
        self.CreateStatusBar()
        self.update_gui()
        self.Show()

    # __init__()

    def onClose(self, ev):
        if self.data_sync_thread.is_running():
            if wx.MessageDialog(None, "Syncing. Are you sure to exit?",
                                "Warning", style=wx.OK|wx.CANCEL).ShowModal() != wx.ID_OK:
                return
            else:
                self.btnStop_click()

        self.Destroy()
    # onClose()

    def update_gui(self):
        if self.data_sync_thread.is_running():
            self.btnStop.Enable()
            self.btnStopNow.Enable()
            self.btnStart.Disable()
            self.btnSourceDir.Disable()
            self.btnDestDir.Disable()
            self.btnUnmount.Disable()
            self.cbCopyShikaResults.Disable()
            self.cbCopyScanOrg.Disable()
            self.cbCopyScanHit.Disable()
            self.txtSourceDir.SetEditable(False)
            self.txtDestDir.SetEditable(False)
            self.txtExcDirs.SetEditable(False)
            self.txtExcFiles.SetEditable(False)
            #self.log_update_timer.Start(500)
            self.SetStatusText("Syncing..")
        else: # stopped
            self.btnStop.Disable()
            self.btnStopNow.Disable()
            self.btnStart.Enable()
            self.btnSourceDir.Enable()
            self.btnDestDir.Enable()
            self.btnUnmount.Enable()
            self.cbCopyShikaResults.Enable()
            self.cbCopyScanOrg.Enable()
            self.cbCopyScanHit.Enable()
            self.txtSourceDir.SetEditable(True)
            self.txtDestDir.SetEditable(True)
            self.txtExcDirs.SetEditable(True)
            self.txtExcFiles.SetEditable(True)
            #self.log_update_timer.Stop()
            self.SetStatusText("Not syncing.")
    # update_gui()

    def on_update_timer(self, ev):
        if not self.data_sync_thread.is_running():
            self.update_gui()
            return

        ddir = self.data_sync_thread.opts.ddir
        space_left = util.check_disk_free_bytes(ddir)
        self.lblDestLeft.SetLabel(" (Space left: %.2f %s)" % util.human_readable_bytes(space_left))
        
        bytes_sofar = self.data_sync_thread.get_total_bytes_sofar()
        self.lblCopiedSize.SetLabel(" (Total transferred size: %.2f %s)" % util.human_readable_bytes(bytes_sofar))
    # on_update_timer()

    def on_log_update_timer(self, ev=None):
        #if not self.data_sync_thread.is_running(): return
        if not self.data_sync_thread.opts.sdir: return

        log_f = os.path.join(self.data_sync_thread.opts.sdir, logfile_name)
        if not os.path.isfile(log_f): return

        size = os.path.getsize(log_f)
        append_str = ""

        if log_f not in self._logfile_pos or self._logfile_pos[log_f] > size:
            append_str = open(log_f).read(size)
        elif self._logfile_pos[log_f] == size:
            return
        else:
            f = open(log_f)
            pos = self._logfile_pos[log_f]
            f.seek(pos)
            append_str = f.read(size - pos)

        self._logfile_pos[log_f] = size

        if append_str:
            # Remove first lines if exceeded the limit
            max_lines = 100000
            curr_lines = self.txtLog.GetNumberOfLines()
            if curr_lines > max_lines:
                lls = [self.txtLog.GetLineLength(x) for x in range(curr_lines-3*max_lines//2)]
                self.txtLog.Remove(0, sum(lls)+len(lls))
                self.txtLog.SetInsertionPointEnd()

            self.txtLog.AppendText(append_str)

            tmp = [x for x in append_str.splitlines(True) if x.endswith("\n") and x.strip()]
            if tmp: self.SetStatusText(tmp[-1].strip().replace("# ", ""))

    # on_log_update_timer()

    def btnSourceDir_click(self, ev):
        current_dir = self.txtSourceDir.GetValue()
        if not os.path.isdir(current_dir): current_dir = ""

        dlg = wx.DirDialog(None, message="Choose a source directory", defaultPath=current_dir)
        if dlg.ShowModal() == wx.ID_OK:
            dirsel = dlg.GetPath()
            self.txtSourceDir.SetValue(dirsel)
        dlg.Destroy()
    # btnSourceDir_click()

    def btnDestDir_click(self, ev):
        current_dir = self.txtDestDir.GetValue()
        dlg = wx.DirDialog(None, message="Choose a destination directory", defaultPath=current_dir)
        if dlg.ShowModal() == wx.ID_OK:
            dirsel = dlg.GetPath()
            self.txtDestDir.SetValue(dirsel)
            self.txtDestDir_enter()
        dlg.Destroy()
    # btnDestDir_click()

    def btnUnmount_click(self, ev):
        def df(d):
            ret, out, err = util.call("df", '"%s"'%d) # XXX Linux specific?
            out_lines = [x for x in [x.strip() for x in out.splitlines()] if x]
            if len(out_lines) < 2: return None
            sp = out_lines[1].split()
            if not sp: return None
            return sp[0]

        ddir = self.txtDestDir.GetValue()
        if not ddir or not os.path.isabs(ddir) or not os.path.isdir(ddir):
            wx.MessageDialog(None, "Invalid directory.", "Error", style=wx.OK).ShowModal()
            return

        mounted_dev = df(ddir)

        if not mounted_dev:
            wx.MessageDialog(None, "Invalid directory.", "Error", style=wx.OK).ShowModal()
            return

        if wx.MessageDialog(None, "Unmounting %s. Are you sure?" % mounted_dev,
                            "Confirmation", style=wx.OK|wx.CANCEL).ShowModal() != wx.ID_OK:
            return

        ret, out, err = util.call("udisks", "--unmount %s" % mounted_dev) # XXX Linux specific?

        if out:
            wx.MessageDialog(None, out, "Error", style=wx.OK).ShowModal()
        elif mounted_dev == df(ddir):
            wx.MessageDialog(None, "Unmount failed.", "Error", style=wx.OK).ShowModal()
        else:
            wx.MessageDialog(None, "Successfully unmounted.", "OK", style=wx.OK).ShowModal()
    # btnUnmount_click()

    def txtSource_enter(self, ev=None):
        sdir = str(self.txtSourceDir.GetValue())

        if not sdir:
            wx.MessageDialog(None, "Please specify source directory.", "Error", style=wx.OK).ShowModal()
            return False

        if not os.path.isabs(sdir):
            wx.MessageDialog(None, "Please specify absolute path for source directory.", "Error", style=wx.OK).ShowModal()
            return False

        if not os.path.isdir(sdir):
            wx.MessageDialog(None, "Please specify source directory, not file", "Error", style=wx.OK).ShowModal()
            return False

        if not os.access(sdir, os.R_OK|os.X_OK):
            wx.MessageDialog(None, "You can't read the source directory.", "Error", style=wx.OK).ShowModal()
            return False

        sdir_norm = os.path.normpath(sdir)
        if sdir != sdir_norm: self.txtSourceDir.SetValue(sdir_norm)

        return True
    # txtSource_enter()

    def txtDestDir_enter(self, ev=None):
        ddir = str(self.txtDestDir.GetValue())
        self.lblDestLeft.SetLabel(" (Space left: ??? GB)")

        if not ddir:
            wx.MessageDialog(None, "Please specify destination directory.", "Error", style=wx.OK).ShowModal()
            return False

        if not os.path.isabs(ddir):
            wx.MessageDialog(None, "Please specify absolute path for destination directory.", "Error", style=wx.OK).ShowModal()
            return False

        if not os.path.exists(ddir):
            try:
                os.makedirs(ddir)
            except:
                wx.MessageDialog(None, "Creation of %s failed"%ddir, "Error", style=wx.OK).ShowModal()
                return False

        if not os.path.isdir(ddir):
            wx.MessageDialog(None, "Please specify destination directory, not file", "Error", style=wx.OK).ShowModal()
            return False

        ddir_norm = os.path.normpath(ddir)
        if ddir != ddir_norm: self.txtDestDir.SetValue(ddir_norm)

        space_left = util.check_disk_free_bytes(ddir)
        self.lblDestLeft.SetLabel(" (Space left: %.2f %s)" % util.human_readable_bytes(space_left))

        if not os.access(ddir, os.W_OK):
            wx.MessageDialog(None, "You don't have write permission in the destination directory.", "Error", style=wx.OK).ShowModal()
            return False

        return True
    # txtDestDir_enter()

    def prepare_opts(self):
        sdir = os.path.normpath(str(self.txtSourceDir.GetValue()))
        ddir = os.path.normpath(str(self.txtDestDir.GetValue()))
        if os.path.basename(ddir) != os.path.basename(sdir):
            # this would be user-friendly..
            ddir = os.path.join(ddir, os.path.basename(sdir))

        if not os.path.exists(ddir): os.mkdir(ddir)

        return SyncOptions(sdir=sdir,
                           ddir=ddir,
                           exclude_dirs=[x for x in self.txtExcDirs.GetValue().splitlines() if x],
                           exclude_files=[x for x in self.txtExcFiles.GetValue().splitlines() if x],
                           copy_shika_results=self.cbCopyShikaResults.GetValue(),
                           copy_scan_data=self.cbCopyScanOrg.GetValue(),
                           copy_hits=self.cbCopyScanHit.GetValue())
    # prepare_opts()

    def btnStart_click(self, ev):
        if not self.txtSource_enter(): return
        if not self.txtDestDir_enter(): return

        sdir = str(self.txtSourceDir.GetValue())
        ddir = str(self.txtDestDir.GetValue())
        
        if os.path.realpath(sdir) == os.path.realpath(ddir):
            wx.MessageDialog(None, "You cannot specify the same directory", "Error", style=wx.OK).ShowModal()            
            return


        if (os.path.realpath(ddir)+os.sep).startswith(os.path.realpath(sdir)+os.sep):
            wx.MessageDialog(None, "Destination path is a subdirectory of source path!", "Error", style=wx.OK).ShowModal()            
            return

        if os.stat(sdir).st_dev == os.stat(ddir).st_dev:
            if wx.MessageDialog(None, "Two directories may be in the same device. Proceed?",
                                "Warning", style=wx.OK|wx.CANCEL).ShowModal() != wx.ID_OK:
                return

        opts = self.prepare_opts()
        self.data_sync_thread.start(opts=opts)
        self.update_gui()
    # btnStart_click()

    def btnStop_click(self, ev=None, stop_now=True):
        if not self.data_sync_thread.is_running():
            wx.MessageDialog(None, "Not running now", "Info", style=wx.OK).ShowModal()
            return False

        busyinfo = wx.lib.agw.pybusyinfo.PyBusyInfo("Stopping..", title="Wait")
        try: wx.SafeYield()
        except: pass

        try:
            self.data_sync_thread.stop(stop_now)
        finally:
            busyinfo = None

        self.on_log_update_timer() # update log
        self.update_gui()
    # btnStop_click()

    def on_sync_error(self, ev):
        self.update_gui()
        wx.MessageDialog(None, ev.msg, "Error", style=wx.OK).ShowModal()
    # on_sync_error()

    def on_sync_done(self, ev):
        max_rows = 10000
        if self.lstSummary.GetItemCount() > max_rows: # just in case
            for i in range(self.lstSummary.GetItemCount()-max_rows):
                self.lstSummary.DeleteItem(0)

        n_copied = len(ev.files_copied)
        if "./%s" % logfile_name in ev.files_copied: n_copied -= 1
        
        pos = self.lstSummary.InsertStringItem(self.lstSummary.GetItemCount(),
                                               time.strftime("%Y-%m-%d %H:%M:%S"))
        self.lstSummary.SetStringItem(pos, 1, "%.2f %s"%(util.human_readable_bytes(ev.bytes_total)))
        self.lstSummary.SetStringItem(pos, 2, "%d"%n_copied)
        self.lstSummary.SetStringItem(pos, 3, "%d"%len(ev.skipped))
        self.lstSummary.EnsureVisible(pos)
    # on_sync_done()
# class MainFrame

def run(argv):
    import optparse
    parser = optparse.OptionParser(usage="usage: %prog [options]")

    parser.add_option("--bl", action="store", dest="bl", default="BL32XU", help="For Window title and logfile name")

    opts, args = parser.parse_args(argv)

    global logfile_name
    logfile_name = "%s_datasync.log" % opts.bl.lower()

    app = wx.App()
    app.TopWindow = MainFrame(opts)
    app.MainLoop()
# run()

if __name__ == "__main__":
    import sys
    run(sys.argv[1:])
