"""
(c) RIKEN 2016-2018. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import wx
import sys
import os
import glob
import re
import socket
import subprocess
import time
import tempfile
import json
import traceback
import getpass
import h5py
import urllib.request, urllib.parse, urllib.error
import numpy
import wx.lib.agw.pybusyinfo
from yamtbx.dataproc import eiger
from yamtbx.dataproc import adxv
from yamtbx.dataproc.XIO.plugins import eiger_hdf5_interpreter

def parse_range(s):
    ret = []
    if not s: return ret
    for x in s.split(","):
        if x.count("-") > 1: return None
        tmp = list(map(int, x.split("-")))
        if len(tmp) == 1: ret.append(tmp[0])
        else: ret.extend(list(range(tmp[0], tmp[1]+1)))
    return sorted(set(ret))
# parse_range()

def make_range_str(frames):
    # frames must be unique and sorted
    if not frames: return ""
    s = str(frames[0])
    if len(frames) == 1: return s
    diffs = [frames[i]-frames[i-1] for i in range(1, len(frames))]
    for i, d in enumerate(diffs):
        if d > 1:
            if i>0 and diffs[i-1] == 1: s += "-%d,%d"%(frames[i],frames[i+1])
            else: s += ",%d"%frames[i+1]
        
    if diffs[-1] == 1: s += "-%d" % frames[-1]
    return s
# make_range_str()

class MainFrame(wx.Frame):
    def __init__(self, parent=None, id=wx.ID_ANY, h5in=None, eiger_host="192.168.163.204", eiger_api_ver="1.6.1", bl="BL32XU", monitor_interval=5.):
        wx.Frame.__init__(self, parent=parent, id=id, title="Adxv launcher for Eiger",
                          size=(540,500))
        self.h5file = None
        self._lstTrigger_frame_request_workaround = None # dirty workaround..
        if h5in:
            if os.path.isdir(h5in):
                h5files = [x for x in sorted(glob.glob(os.path.join(h5in, "*.h5"))) if x.endswith(("_master.h5", "_onlyhits.h5"))]
                if h5files: self.set_h5file(h5files[0])
            else:
                self.set_h5file(h5in)
                
        self.onlyhits_keys = []
        self.n_images = -1
        self.n_images_each = -1
        self.adxv = adxv.Adxv()
        self.last_monitor_image = None
        self.bl = bl

        vbox = wx.BoxSizer(wx.VERTICAL)

        notebook = wx.Notebook(self, wx.ID_ANY)
        panel1 = wx.Panel(notebook, wx.ID_ANY)
        panel2 = wx.Panel(notebook, wx.ID_ANY)

        # panel1
        vbox0 = wx.BoxSizer(wx.VERTICAL)
        hbox00 = wx.BoxSizer(wx.HORIZONTAL)
        hbox00.Add(wx.StaticText(panel1, wx.ID_ANY, "File: "), flag=wx.LEFT|wx.ALIGN_CENTER_VERTICAL)
        self.cmbInfile = wx.ComboBox(panel1, wx.ID_ANY, size=(350,25))
        self.cmbInfile.Bind(wx.EVT_COMBOBOX, self.cmbInfile_onChange)
        self.cmbInfile.Bind(wx.EVT_TEXT_ENTER, self.cmbInfile_onChange)
        hbox00.Add(self.cmbInfile, 1, flag=wx.EXPAND|wx.RIGHT)
        self.btnInfile = wx.Button(panel1, wx.ID_ANY, "...", size=(25,25))
        self.btnInfile.Bind(wx.EVT_BUTTON, self.btnInfile_click)
        self.btnLatest = wx.Button(panel1, wx.ID_ANY, "Latest", size=(50,25))
        self.btnLatest.Bind(wx.EVT_BUTTON, self.btnLatest_click)
        hbox00.Add(self.btnInfile)
        hbox00.Add(self.btnLatest)
        vbox0.Add(hbox00, flag=wx.EXPAND|wx.TOP, border=4)

        hbox01 = wx.BoxSizer(wx.HORIZONTAL)
        vbox0.Add(hbox01, flag=wx.EXPAND|wx.TOP, border=4)

        self.lstTrigger = wx.ListCtrl(panel1, size=(180, 200), style=wx.LC_REPORT|wx.LC_SINGLE_SEL)
        self.lstTrigger.Bind(wx.EVT_LIST_ITEM_SELECTED, self.lstTrigger_onSelected)
        hbox01.Add(self.lstTrigger)

        vbox012 = wx.BoxSizer(wx.VERTICAL)
        hbox01.Add(vbox012, flag=wx.EXPAND|wx.LEFT, border=4)
        hbox0121 = wx.BoxSizer(wx.HORIZONTAL)
        vbox012.Add(hbox0121, flag=wx.EXPAND)
        self.spnFrame = wx.SpinCtrlDouble(panel1, wx.ID_ANY, size=(100,25), min=1)
        self.spnFrame.SetDigits(0) # this is integer actually. but I want arbitrary increments
        #self.sliFrame = wx.Slider(panel1)
        self.spnFrame.Bind(wx.EVT_SPINCTRLDOUBLE, self.onTextEnter)
        self.stFrameRange = wx.StaticText(panel1, wx.ID_ANY, "")
        hbox0121.Add(wx.StaticText(panel1, wx.ID_ANY, "Frame no: "), flag=wx.LEFT|wx.ALIGN_CENTER_VERTICAL)
        hbox0121.Add(self.spnFrame)#, flag=wx.EXPAND|wx.RIGHT)
        hbox0121.Add(self.stFrameRange, flag=wx.LEFT|wx.ALIGN_CENTER_VERTICAL)
        #hbox0121.Add(self.sliFrame, 1, flag=wx.EXPAND|wx.RIGHT)

        self.stShowing = wx.StaticText(panel1, wx.ID_ANY, "Nothing shown")
        vbox012.Add(self.stShowing, flag=wx.TOP, border=4)
        
        hbox0122 = wx.BoxSizer(wx.HORIZONTAL)
        vbox012.Add(hbox0122, flag=wx.TOP, border=4)
        self.chkSum = wx.CheckBox(panel1, label="Sum")
        self.chkSum.Bind(wx.EVT_CHECKBOX, self.chkSum_onChecked)
        self.spnSumFrame = wx.SpinCtrl(panel1, size=(80, 25))
        self.spnSumFrame.Bind(wx.EVT_SPINCTRL, self.onSpnSumFrame)
        self.spnSumDeg = wx.SpinCtrlDouble(panel1, size=(80, 25))
        self.spnSumDeg.Bind(wx.EVT_SPINCTRLDOUBLE, self.onSpnSumDeg)
        self.spnSumDeg.SetDigits(4)
        hbox0122.Add(self.chkSum)
        hbox0122.Add(self.spnSumFrame)
        hbox0122.Add(wx.StaticText(panel1, wx.ID_ANY, " frames or "), flag=wx.LEFT|wx.ALIGN_CENTER_VERTICAL)
        hbox0122.Add(self.spnSumDeg)
        hbox0122.Add(wx.StaticText(panel1, wx.ID_ANY, " degrees"), flag=wx.LEFT|wx.ALIGN_CENTER_VERTICAL)

        self.chkKeepFrame = wx.CheckBox(panel1, label="Keep frame number when trigger changed")
        vbox012.Add(self.chkKeepFrame, flag=wx.TOP, border=4)
        
        hbox0123 = wx.BoxSizer(wx.HORIZONTAL)
        vbox012.Add(hbox0123, flag=wx.TOP, border=4)
        self.llbtn = wx.Button(panel1, wx.ID_ANY, "<<")
        self.lbtn = wx.Button(panel1, wx.ID_ANY, "<")
        self.rbtn = wx.Button(panel1, wx.ID_ANY, ">")
        self.rrbtn = wx.Button(panel1, wx.ID_ANY, ">>")
        hbox0123.Add(self.llbtn, flag=wx.EXPAND|wx.RIGHT)
        hbox0123.Add(self.lbtn, flag=wx.EXPAND|wx.RIGHT)
        hbox0123.Add(self.rbtn, flag=wx.EXPAND|wx.RIGHT)
        hbox0123.Add(self.rrbtn, flag=wx.EXPAND|wx.RIGHT)
        self.lbtn.Bind(wx.EVT_BUTTON, lambda e: self.next_or_back(-1))
        self.rbtn.Bind(wx.EVT_BUTTON, lambda e: self.next_or_back(+1))
        self.llbtn.Bind(wx.EVT_BUTTON, lambda e: self.play(-1))
        self.rrbtn.Bind(wx.EVT_BUTTON, lambda e: self.play(+1))

        hbox0124 = wx.BoxSizer(wx.HORIZONTAL)
        vbox012.Add(hbox0124, flag=wx.TOP, border=4)
        hbox0124.Add(wx.StaticText(panel1, wx.ID_ANY, "Binning: "), flag=wx.LEFT|wx.ALIGN_CENTER_VERTICAL)
        self.txtBin = wx.TextCtrl(panel1, wx.ID_ANY, size=(100,25), style=wx.TE_PROCESS_ENTER)
        self.txtBin.SetValue("1")
        self.txtBin.Bind(wx.EVT_TEXT_ENTER, self.onTextEnter)
        hbox0124.Add(self.txtBin, flag=wx.EXPAND|wx.RIGHT)

        panel1.SetSizer(vbox0)
        notebook.InsertPage(0, panel1, "File")

        # panel2
        vbox1 = wx.BoxSizer(wx.VERTICAL)
        #hbox10 = wx.BoxSizer(wx.HORIZONTAL)
        #hbox10.Add(wx.StaticText(panel2, wx.ID_ANY, "Wavelength: "), flag=wx.LEFT|wx.ALIGN_CENTER_VERTICAL)
        #self.txtMonWavelen = wx.TextCtrl(panel2, wx.ID_ANY, size=(100,25))
        #self.txtMonWavelen.SetValue("1.0000")
        #hbox10.Add(self.txtMonWavelen, flag=wx.EXPAND|wx.RIGHT)
        #hbox10.Add(wx.StaticText(panel2, wx.ID_ANY, " Distance: "), flag=wx.LEFT|wx.ALIGN_CENTER_VERTICAL)
        #self.txtMonDistance = wx.TextCtrl(panel2, wx.ID_ANY, size=(100,25))
        #self.txtMonDistance.SetValue("300.0")
        #hbox10.Add(self.txtMonDistance, flag=wx.EXPAND|wx.RIGHT)
        #vbox1.Add(hbox10, flag=wx.EXPAND|wx.TOP, border=4)

        hbox11 = wx.BoxSizer(wx.HORIZONTAL)
        #hbox11.Add(wx.StaticText(panel2, wx.ID_ANY, "Beam center: "), flag=wx.LEFT|wx.ALIGN_CENTER_VERTICAL)
        #self.txtMonBeamxy = wx.TextCtrl(panel2, wx.ID_ANY, size=(150,25))
        #self.txtMonBeamxy.SetValue("1550, 1634.5")
        #hbox11.Add(self.txtMonBeamxy, flag=wx.EXPAND|wx.RIGHT)
        self.btnMonRefresh = wx.Button(panel2, wx.ID_ANY, "Start monitoring")
        self.btnMonRefresh.Bind(wx.EVT_BUTTON, self.on_btnMonRefresh)
        self.btnMonStop = wx.Button(panel2, wx.ID_ANY, "Stop")
        self.btnMonStop.Bind(wx.EVT_BUTTON, lambda e: self.monitor_timer.Stop())
        self.txtMonInterval = wx.TextCtrl(panel2, wx.ID_ANY, "%.1f"%monitor_interval, size=(50,25))
        hbox11.Add(self.btnMonRefresh, flag=wx.EXPAND|wx.RIGHT)
        hbox11.Add(self.btnMonStop, flag=wx.EXPAND|wx.RIGHT)
        hbox11.Add(wx.StaticText(panel2, wx.ID_ANY, " interval: "), flag=wx.LEFT|wx.ALIGN_CENTER_VERTICAL)
        hbox11.Add(self.txtMonInterval, flag=wx.EXPAND|wx.RIGHT)
        hbox11.Add(wx.StaticText(panel2, wx.ID_ANY, " sec."), flag=wx.LEFT|wx.ALIGN_CENTER_VERTICAL)
        vbox1.Add(hbox11, flag=wx.EXPAND|wx.TOP, border=4)

        hbox12 = wx.BoxSizer(wx.HORIZONTAL)
        self.chkMonRaiseW = wx.CheckBox(panel2, wx.ID_ANY, "Raise adxv windows")
        self.chkMonRaiseW.SetValue(1)
        hbox12.Add(self.chkMonRaiseW, flag=wx.EXPAND|wx.RIGHT)
        vbox1.Add(hbox12, flag=wx.EXPAND|wx.TOP, border=4)

        hbox13 = wx.BoxSizer(wx.HORIZONTAL)
        hbox13.Add(wx.StaticText(panel2, wx.ID_ANY, "Eiger DCU host: "), flag=wx.LEFT|wx.ALIGN_CENTER_VERTICAL)
        self.txtEigerHost = wx.TextCtrl(panel2, wx.ID_ANY, eiger_host, size=(150,25))
        hbox13.Add(self.txtEigerHost)
        hbox13.Add(wx.StaticText(panel2, wx.ID_ANY, "  API ver.: "), flag=wx.LEFT|wx.ALIGN_CENTER_VERTICAL)
        self.txtEigerAPIver = wx.TextCtrl(panel2, wx.ID_ANY, eiger_api_ver, size=(100,25))
        hbox13.Add(self.txtEigerAPIver)
        vbox1.Add(hbox13, flag=wx.EXPAND|wx.TOP, border=4)


        panel2.SetSizer(vbox1)
        notebook.InsertPage(1, panel2, "Monitor")
        notebook.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, self.on_notebook_page_changed)

        vbox.Add(notebook, flag=wx.EXPAND)
        self.txtInfo = wx.TextCtrl(self, wx.ID_ANY, style=wx.TE_MULTILINE)
        self.txtInfo.SetFont(wx.Font(12, wx.MODERN, wx.NORMAL, wx.NORMAL))
        self.txtInfo.SetEditable(False)

        vbox.Add(self.txtInfo, 1, flag=wx.EXPAND, border=4)

        self.SetSizer(vbox)

        self.play_timer = wx.Timer(self)
        self.play_relval = 0
        self.Bind(wx.EVT_TIMER, self.on_play_timer, self.play_timer)

        self.monitor_timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.on_monitor_timer, self.monitor_timer)

        self.Show()

        if self.h5file is not None:
            self.read_h5file()
    # __init__()

    def set_h5file(self, f):
        if not f.endswith(("_master.h5", "_onlyhits.h5")):
            r = re.search("(.*)_data_[0-9]+\.h5$", f)
            if not r:
                wx.MessageDialog(None, "Choose master h5 file (*_master.h5) or Hit-only h5 file (*_onlyhits.h5)",
                                 "Error", style=wx.OK).ShowModal()
                return False
            f = r.group(1)+"_master.h5"
            
        if not os.path.isfile(f):
            wx.MessageDialog(None, "Selected file does not exist",
                             "Error", style=wx.OK).ShowModal()
            return False

        self.h5file = f
        return True
    # set_h5file()
    
    def on_notebook_page_changed(self, ev):
        self.txtInfo.SetValue("")

        if ev.Selection == 0: # File mode
            self.monitor_timer.Stop()
    # on_notebook_page_changed()

    def on_btnMonRefresh(self, ev):
        self.monitor_timer.Stop()
        self.on_monitor_timer(None)
        try:
            msec = float(self.txtMonInterval.GetValue()) * 1000
            if msec <= 0 or msec != msec: raise ValueError
        except ValueError:
            wx.MessageDialog(None, "Invalid interval (%s)" % self.txtMonInterval.GetValue(),
                             "Error", style=wx.OK).ShowModal()
            return

        self.monitor_timer.Start(msec)
    # on_btnMonRefresh()

    def open_hdf5(self, frameno, raise_window=True):
        tmpdir = "/dev/shm" if os.path.isdir("/dev/shm") else tempfile.gettempdir()

        #busyinfo = wx.lib.agw.pybusyinfo.PyBusyInfo("Loading image..", title="Busy adxv_eiger")
        focused = self.FindFocus()
        try: wx.SafeYield()
        except: pass

        try:
            self.adxv.open_hdf5(self.h5file, frameno, tmpdir=tmpdir,
                                raise_window=raise_window,
                                binning=int(self.txtBin.GetValue()))
        except RuntimeError as e:
            self.txtInfo.SetValue("Error! %s\n\n%s" % (str(e), self.txtInfo.GetValue()))
        except:
            self.txtInfo.SetValue(traceback.format_exc() + "\n\n" + self.txtInfo.GetValue())
        finally:
            busyinfo = None

        if focused:
            focused.SetFocus()

    # open_hdf5()

    def change_frameno(self, frames):
        if self.onlyhits_keys:
            frames = "/entry/data/%s/data"%self.onlyhits_keys[frames[0]-1] # No support for summation
        else:
            if all([x < 1 for x in frames]): frames = [1]
            if all([x > self.n_images for x in frames]): frames = [self.n_images]
            frames = [x for x in frames if 0 < x <= self.n_images]
            if not frames: frames = [1]

        self.open_hdf5(frames, raise_window=False)
    # change_frameno()

    def chkSum_onChecked(self, ev):
        if self.chkSum.GetValue():
            nframe = int(self.spnSumFrame.GetValue())
            self.spnFrame.SetRange(1, self.n_images_each-nframe+1)
            self.spnFrame.SetIncrement(nframe)
        else:
            self.spnFrame.SetRange(1, self.n_images_each)
            self.spnFrame.SetIncrement(1)

        self.onTextEnter(None)
    # chkSum_onChecked()
        
    def onSpnSumFrame(self, ev):
        nframe = int(self.spnSumFrame.GetValue())
        inc = self.spnSumDeg.GetIncrement()
        self.spnSumDeg.SetValue(inc*nframe)
        if self.chkSum.GetValue():
            self.spnFrame.SetRange(1, self.n_images_each-nframe+1)
            self.spnFrame.SetIncrement(nframe)
            self.onTextEnter(None)
    # onSpnSumFrame()

    def onSpnSumDeg(self, ev):
        deg = float(self.spnSumDeg.GetValue())
        inc = self.spnSumDeg.GetIncrement()
        nframe = int(deg/inc+.5)
        self.spnSumFrame.SetValue(nframe)
        if self.chkSum.GetValue():
            self.spnFrame.SetRange(1, self.n_images_each-nframe+1)
            self.spnFrame.SetIncrement(nframe)
            self.onTextEnter(None)
    # onSpnSumDeg()

    def lstTrigger_onSelected(self, ev):
        if not self.chkKeepFrame.GetValue():
            self.spnFrame.SetValue(1)

        if self._lstTrigger_frame_request_workaround is not None:
            self.spnFrame.SetValue(self._lstTrigger_frame_request_workaround)
            self._lstTrigger_frame_request_workaround = None
            
        self.lstTrigger.EnsureVisible(self.lstTrigger.GetFirstSelected())
        self.onTextEnter(None)
    # lstTrigger_onSelected()
    
    def onTextEnter(self, ev):
        if self.onlyhits_keys:
            frames = [self.lstTrigger.GetFirstSelected()+1]
            self.stShowing.SetLabel("Showing %s"%self.onlyhits_keys[frames[0]-1])
        else:
            trigger = self.lstTrigger.GetFirstSelected() +1
            frame_sp = int(self.spnFrame.GetValue())
            frame = (trigger-1)*self.n_images_each + frame_sp
            sumframe = self.spnSumFrame.GetValue() if self.chkSum.GetValue() else 1
            max_frame = trigger*self.n_images_each
            frames = list(range(frame, min(frame+sumframe-1, max_frame)+1))
            frame_str = ("frame %d" % frame_sp) if len(frames)<2 else ("frames %d..%d" % (frame_sp, frame_sp+len(frames)-1))
            self.stShowing.SetLabel("Showing %s of trigger %d" % (frame_str, trigger))

        self.change_frameno(frames)
    # onTextEnter()

    def next_or_back(self, sign):
        tf = self.lstTrigger.GetFirstSelected() +1
        tmin, tmax = 1, self.lstTrigger.GetItemCount()
        
        if self.onlyhits_keys:
            if tmin <= tf+sign <= tmax:
                self.lstTrigger.Select(tf+sign-1)
        else:
            cf, incf = self.spnFrame.GetValue(), self.spnFrame.GetIncrement()
            fmin, fmax = self.spnFrame.GetMin(), self.spnFrame.GetMax()
            
            if cf+sign*incf < fmin:
                if tf > tmin:
                    #self.spnFrame.SetValue(fmax) # this does not work, because listctrl.Select() calls an event, which changes spnFrame value.
                    self._lstTrigger_frame_request_workaround = fmax
                    self.lstTrigger.Select(tf-1-1)
            elif fmax < cf+sign*incf:
                if tf < tmax:
                    #self.spnFrame.SetValue(fmin)
                    self._lstTrigger_frame_request_workaround = fmin
                    self.lstTrigger.Select(tf+1-1)
            else: #if fmin <= cf+sign*incf <= fmax:
                self.spnFrame.SetValue(cf+sign*incf)
                self.onTextEnter(None)
    # next_or_back()

    def play(self, rel):
        if self.play_timer.IsRunning():
            self.play_timer.Stop()
            if rel*self.play_relval > 0:
                self.play_stop()
                return

        if rel > 0: self.rrbtn.SetLabel("||")
        else: self.llbtn.SetLabel("||")

        self.play_relval = rel
        self.play_timer.Start(200)
    # play()

    def on_play_timer(self, ev):
        save = (self.spnFrame.GetValue(), self.lstTrigger.GetFirstSelected())
        self.next_or_back(self.play_relval)
        try: wx.Yield()
        except: pass

        if save == (self.spnFrame.GetValue(), self.lstTrigger.GetFirstSelected()):
            self.play_stop()
    # on_play_timer()

    def play_stop(self):
        self.play_timer.Stop()
        self.llbtn.SetLabel("<<")
        self.rrbtn.SetLabel(">>")
    # play_stop()

    def btnInfile_click(self, ev):
        dfile, ddir = "", ""
        if self.h5file:
            if os.path.exists(self.h5file): dfile = self.h5file
            else:
                ddir = os.path.dirname(self.h5file)
        else:
            ddir = os.getcwd()

        if not os.path.exists(ddir): ddir = ""
        
        dlg = wx.FileDialog(None, message="Choose a master.h5",
                            defaultFile=dfile, defaultDir=ddir,
                            wildcard="Mater HDF5 (*_master.h5 or *_onlyhits.h5)|*.h5")

        if dlg.ShowModal() == wx.ID_OK:
            if self.set_h5file(dlg.GetPath()):
                self.read_h5file()
            
        dlg.Destroy()
    # btnInfile_click()

    def btnLatest_click(self, ev):
        latestlog = os.path.join(os.environ["HOME"], ".bss_latest_file_%s.log" % self.bl)
        if not os.path.isfile(latestlog):
            wx.MessageDialog(None, "%s not found." % latestlog,
                             "Error", style=wx.OK).ShowModal()
            return

        latestfile = open(latestlog).read().strip()
        latestfile = os.path.join(os.path.dirname(latestfile), re.sub("_[0-9]+\..*$", "_master.h5", os.path.basename(latestfile)))

        if not os.path.isfile(latestfile):
            wx.MessageDialog(None, "Latest file (%s) not found." % latestfile,
                             "Error", style=wx.OK).ShowModal()
            return

        if self.set_h5file(latestfile):
            self.read_h5file()
    # btnLatest_click ()

    def cmbInfile_onChange(self, ev):
        f = self.cmbInfile.GetValue()
        if self.set_h5file(f):
            self.read_h5file()
    # cmbInfile_onChange()


    def read_h5file(self):
        self.cmbInfile.Clear()
        glob_f = lambda x: glob.glob(os.path.join(os.path.dirname(self.h5file), x))
        for f in sorted(glob_f("*_master.h5")+glob_f("*_onlyhits.h5")): self.cmbInfile.Append(f)
        self.cmbInfile.SetStringSelection(self.h5file)
        self.lstTrigger.ClearAll()

        widgets = self.spnFrame, self.chkSum, self.spnSumFrame, self.spnSumDeg, self.chkKeepFrame
        
        is_onlyhits = self.h5file.endswith("_onlyhits.h5")
        try:
            h5 = h5py.File(self.h5file, "r")
            h = eiger_hdf5_interpreter.Interpreter().getRawHeadDict(self.h5file)
            self.n_images = h["Nimages"]
            self.n_images_each = h["Nimages_each"]
        except:
            self.n_images = self.n_images_each = 0
            self.txtInfo.SetValue(traceback.format_exc())
            return

        try:
            safe_val = lambda x, k, a: x[k][()] if k in x else a
            yes_or_no = lambda x, k: ("no", "yes", "?")[safe_val(x, k, 2)]
            self.txtInfo.SetValue("""\
         Filename: %s
         Detector: %s (%s)
 Number of frames: %d (%d * %d trigger)
             Date: %s
       Omega step: %.4f deg
  Camera distance: %.2f mm
       Wavelength: %.5f A
       Frame time: %.6f sec 
      Beam center: (%d, %d) px
      Corrections: countrate=%s 
                   flatfield=%s
                   efficiency=%s
         ROI mode: %s
           Sensor: %s %.1f um
     Eiger fw ver: %s
    """ % (os.path.basename(self.h5file),
           safe_val(h5["/entry/instrument/detector"], "description", "?"), safe_val(h5["/entry/instrument/detector"], "detector_number", "?"),
           self.n_images, self.n_images_each, h["Ntrigger"],
           h["DateStr"].replace(b"T", b" "),
           h["PhiWidth"],
           h["Distance"]*1.e3,
           h["Wavelength"],
           safe_val(h5["/entry/instrument/detector"], "frame_time", float("nan")),
           h["BeamX"], h["BeamY"],
           yes_or_no(h5["/entry/instrument/detector"], "countrate_correction_applied"),
           yes_or_no(h5["/entry/instrument/detector"], "flatfield_correction_applied"),
           yes_or_no(h5["/entry/instrument/detector"], "efficiency_correction_applied"),
           safe_val(h5["/entry/instrument/detector/detectorSpecific"], "roi_mode", "?"),
           h.get("SensorMaterial", "?"), h.get("SensorThickness", float("nan"))*1.e6,
           safe_val(h5["/entry/instrument/detector/detectorSpecific"], "eiger_fw_version", "?")
    ))
        except:
            self.txtInfo.SetValue(traceback.format_exc())

        if is_onlyhits:
            for w in widgets: w.Disable()
            self.stFrameRange.SetLabel("")
            self.lstTrigger.InsertColumn(0, "Frame", width=55)
            self.lstTrigger.InsertColumn(1, "Spots", width=55)
            self.onlyhits_keys = sorted(h5["/entry/data"].keys())
            if not self.onlyhits_keys: return
            for i, v in enumerate(self.onlyhits_keys):
                nspots = h5["/entry/data"][v]["data"].attrs.get("n_spots")
                nspots = "%3d" % nspots if nspots is not None else ""
                frame = int(v[v.rindex("_")+1:])
                self.lstTrigger.InsertStringItem(i, str(frame))
                self.lstTrigger.SetStringItem(i, 1, nspots)
            self.lstTrigger.Select(0) # this calls an event to load image

        else:
            for w in widgets: w.Enable()
            self.onlyhits_keys = []
            self.stFrameRange.SetLabel(" (from 1 to %d)"%self.n_images_each)
            self.spnFrame.SetRange(1, self.n_images_each)
            self.spnFrame.SetIncrement(1)
            self.spnSumFrame.SetRange(1, self.n_images_each)
            self.spnSumDeg.SetIncrement(h["PhiWidth"])
            self.spnSumDeg.SetRange(h["PhiWidth"], h["PhiWidth"]*self.n_images_each)
            if h["PhiWidth"] > 1e-5:
                self.spnSumDeg.Enable()
                self.spnSumFrame.SetValue(int(1./h["PhiWidth"]+.5))
                self.onSpnSumFrame(None)
            else:
                self.spnSumFrame.SetValue(1)
                self.spnSumDeg.Disable()
                
            self.lstTrigger.InsertColumn(0, "Trigger")
            self.lstTrigger.InsertColumn(1, "Omega")
            omega = h5["/entry/sample/goniometer/omega"]
            for i in range(h["Ntrigger"]):
                self.lstTrigger.InsertStringItem(i, "%d"%(i+1))
                try:
                    oms = omega[i*self.n_images_each]
                    self.lstTrigger.SetStringItem(i, 1, "%.3f"%oms)
                except: pass

            self.lstTrigger.Select(0) # this calls an event to load image
    # read_h5file()

    def on_monitor_timer(self, ev):
        print("monitoring")
        response = urllib.request.urlopen("http://%s/monitor/api/%s/images/monitor" % (self.txtEigerHost.GetValue(), self.txtEigerAPIver.GetValue()))
        data = response.read()
        data_size = len(data)
        curlines = self.txtInfo.GetValue().splitlines()
        if len(curlines) > 20:
            self.txtInfo.SetValue("\n".join(curlines))

        if data_size < 10000:
            self.txtInfo.SetValue("%s: %s\n%s" % (time.ctime(), data, self.txtInfo.GetValue()))
            return
            
        #import cPickle as pickle
        #pickle.dump(data, open("monitor.pkl","wb"), -1)
        #data = pickle.load(open("monitor.pkl"))

        def request(query, base="detector"):
            r = urllib.request.urlopen("http://%s/%s/api/%s/%s"%(self.txtEigerHost.GetValue(), base, self.txtEigerAPIver.GetValue(), query)).read()
            return json.loads(r)

        nx = request("config/x_pixels_in_detector")["value"]
        ny = request("config/y_pixels_in_detector")["value"]

        if len(data) < 1000:
            data = self.last_monitor_image
            if data is None: return
        else:
            byte = data_size // (nx*ny)
            assert byte == 4 or byte == 2
            data = numpy.fromstring(data[8:-102], dtype=numpy.int32 if byte==4 else numpy.int16).reshape(ny,nx)

            if self.last_monitor_image is not None and (data==self.last_monitor_image).all():
                return

        wavelen = request("config/wavelength")["value"]
        beamx = request("config/beam_center_x")["value"]
        beamy = request("config/beam_center_y")["value"]

        distance = request("config/detector_distance")

        image_number = request("status/monitor_image_number", base="monitor")["value"]

        self.txtInfo.SetValue("""\
%s: Downloaded (frame: %.6d)
    Wavelength: %.5f A
    Detector_distance: %.2f %s
    Beam_xy: %.2f, %.2f
%s"""% (time.ctime(), image_number[1]+1, wavelen, distance["value"], str(distance["unit"]), beamx, beamy, self.txtInfo.GetValue()))

        from yamtbx.dataproc import cbf
        tmpdir = "/dev/shm" if os.path.isdir("/dev/shm") else tempfile.gettempdir()
        imgfile = os.path.join("/dev/shm", "adxvtmp-%s-%s.cbf"%(getpass.getuser(), os.getpid()))

        xp = request("config/x_pixel_size")
        yp = request("config/y_pixel_size")

        cbf.save_numpy_data_as_cbf(data.flatten(), size1=data.shape[1], size2=data.shape[0], title="",
                                   cbfout=imgfile,
                                   pilatus_header="""\
# Detector: Eiger
# Pixel_size %(px).5e %(pxu)s x %(py).5e %(pyu)s
# Wavelength %(Wavelength).6f A
# Detector_distance %(Distance).4e %(DistanceU)s
# Beam_xy (%(BeamX).1f, %(BeamY).1f) pixels
""" % dict(px=xp["value"], pxu=str(xp["unit"]), py=yp["value"], pyu=str(yp["unit"]), Wavelength=wavelen, Distance=distance["value"], DistanceU=str(distance["unit"]), BeamX=beamx, BeamY=beamy))

        self.adxv.open_image(imgfile, raise_window=self.chkMonRaiseW.GetValue())
        self.last_monitor_image = data
    # on_monitor_timer()

def run(argv):
    import optparse

    parser = optparse.OptionParser(usage="usage: %prog [options] foo_master.h5")

    parser.add_option("--bl", action="store", dest="bl", default="BL32XU", help="To look for ~/.bss_latest_file_<bl>.log")
    parser.add_option("--dcu-host", action="store", dest="dcu_host", default="192.168.163.204", help="Eiger DCU host (ip addr) for monitor mode")
    parser.add_option("--api-ver", action="store", dest="api_ver", default="1.6.1", help="Eiger API ver.")
    parser.add_option("--monitor-interval", action="store", dest="monitor_interval", default=5, type=float, help="Monitoring interval time in seconds")

    opts, args = parser.parse_args(argv)

    h5in = args[0] if args else None

    app = wx.App()
    app.TopWindow = MainFrame(parent=None, id=wx.ID_ANY, 
                              h5in=h5in,
                              eiger_host=opts.dcu_host, eiger_api_ver=opts.api_ver,
                              bl=opts.bl, monitor_interval=opts.monitor_interval)
    app.MainLoop()

if __name__ == "__main__":
    run(os.path.abspath(sys.argv[1:]))
