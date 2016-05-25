"""
(c) RIKEN 2016. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

import wx
import sys
import os
import re
import socket
import subprocess
import time
import json
import traceback
import getpass
import h5py
import urllib
import numpy
from yamtbx.dataproc import eiger
from yamtbx.dataproc import adxv
from yamtbx.dataproc.XIO.plugins import eiger_hdf5_interpreter

def parse_range(s):
    ret = []
    if not s: return ret
    for x in s.split(","):
        if x.count("-") > 1: return None
        tmp = map(int, x.split("-"))
        if len(tmp) == 1: ret.append(tmp[0])
        else: ret.extend(range(tmp[0], tmp[1]+1))
    return sorted(set(ret))
# parse_range()

def make_range_str(frames):
    # frames must be unique and sorted
    if not frames: return ""
    s = str(frames[0])
    if len(frames) == 1: return s
    diffs = map(lambda i: frames[i]-frames[i-1], xrange(1, len(frames)))
    for i, d in enumerate(diffs):
        if d > 1:
            if i>0 and diffs[i-1] == 1: s += "-%d,%d"%(frames[i],frames[i+1])
            else: s += ",%d"%frames[i+1]
        
    if diffs[-1] == 1: s += "-%d" % frames[-1]
    return s
# make_range_str()

class MainFrame(wx.Frame):
    def __init__(self, parent=None, id=wx.ID_ANY, h5in=None, eiger_host="192.168.163.204"):
        wx.Frame.__init__(self, parent=parent, id=id, title="Adxv launcher for Eiger",
                          size=(500,500))
        self.adxv_proc = None # subprocess object
        self.h5file = h5in
        self.n_images = -1
        self.eiger_host = eiger_host
        self.adxv = adxv.Adxv()
        self.last_monitor_image = None

        vbox = wx.BoxSizer(wx.VERTICAL)

        notebook = wx.Notebook(self, wx.ID_ANY)
        panel1 = wx.Panel(notebook, wx.ID_ANY)
        panel2 = wx.Panel(notebook, wx.ID_ANY)

        # panel1
        vbox0 = wx.BoxSizer(wx.VERTICAL)
        hbox00 = wx.BoxSizer(wx.HORIZONTAL)
        hbox00.Add(wx.StaticText(panel1, wx.ID_ANY, "File: "), flag=wx.LEFT|wx.ALIGN_CENTER_VERTICAL)
        self.txtInfile = wx.TextCtrl(panel1, wx.ID_ANY, size=(350,25))
        self.txtInfile.SetEditable(False)
        hbox00.Add(self.txtInfile, 1, flag=wx.EXPAND|wx.RIGHT)
        self.btnInfile = wx.Button(panel1, wx.ID_ANY, "...", size=(25,25))
        self.btnInfile.Bind(wx.EVT_BUTTON, self.btnInfile_click)
        self.btnLatest = wx.Button(panel1, wx.ID_ANY, "Latest", size=(50,25))
        self.btnLatest.Bind(wx.EVT_BUTTON, self.btnLatest_click)
        hbox00.Add(self.btnInfile)
        hbox00.Add(self.btnLatest)
        vbox0.Add(hbox00, flag=wx.EXPAND|wx.TOP, border=4)

        hbox01 = wx.BoxSizer(wx.HORIZONTAL)
        hbox01.Add(wx.StaticText(panel1, wx.ID_ANY, "Frame no: "), flag=wx.LEFT|wx.ALIGN_CENTER_VERTICAL)
        self.txtFrame = wx.TextCtrl(panel1, wx.ID_ANY, size=(100,25), style=wx.TE_PROCESS_ENTER)
        self.txtFrame.SetEditable(True)
        self.txtFrame.SetValue("1")
        self.txtFrame.Bind(wx.EVT_TEXT_ENTER, self.onTextEnter)
        hbox01.Add(self.txtFrame, flag=wx.EXPAND|wx.RIGHT)
        self.llbtn = wx.Button(panel1, wx.ID_ANY, "<<")
        hbox01.Add(self.llbtn, flag=wx.EXPAND|wx.RIGHT)
        self.lbtn = wx.Button(panel1, wx.ID_ANY, "<")
        hbox01.Add(self.lbtn, flag=wx.EXPAND|wx.RIGHT)
        self.rbtn = wx.Button(panel1, wx.ID_ANY, ">")
        hbox01.Add(self.rbtn, flag=wx.EXPAND|wx.RIGHT)
        self.rrbtn = wx.Button(panel1, wx.ID_ANY, ">>")
        hbox01.Add(self.rrbtn, flag=wx.EXPAND|wx.RIGHT)
        vbox0.Add(hbox01, flag=wx.EXPAND|wx.TOP, border=4)
        #self.lbtn.Bind(wx.EVT_BUTTON, lambda e: self.go_rel(-1))
        #self.rbtn.Bind(wx.EVT_BUTTON, lambda e: self.go_rel(+1))
        #self.llbtn.Bind(wx.EVT_BUTTON, lambda e: self.play(-1))
        #self.rrbtn.Bind(wx.EVT_BUTTON, lambda e: self.play(+1))
        self.lbtn.Bind(wx.EVT_BUTTON, lambda e: self.next_or_back(-1))
        self.rbtn.Bind(wx.EVT_BUTTON, lambda e: self.next_or_back(+1))
        self.llbtn.Bind(wx.EVT_BUTTON, lambda e: self.play(-1))
        self.rrbtn.Bind(wx.EVT_BUTTON, lambda e: self.play(+1))

        hbox02 = wx.BoxSizer(wx.HORIZONTAL)
        hbox02.Add(wx.StaticText(panel1, wx.ID_ANY, "Binning: "), flag=wx.LEFT|wx.ALIGN_CENTER_VERTICAL)
        self.txtBin = wx.TextCtrl(panel1, wx.ID_ANY, size=(100,25), style=wx.TE_PROCESS_ENTER)
        self.txtBin.SetValue("1")
        self.txtBin.Bind(wx.EVT_TEXT_ENTER, self.onTextEnter)
        hbox02.Add(self.txtBin, flag=wx.EXPAND|wx.RIGHT)
        vbox0.Add(hbox02, flag=wx.EXPAND|wx.TOP, border=4)

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
        hbox11.Add(self.btnMonRefresh, flag=wx.EXPAND|wx.RIGHT)
        hbox11.Add(self.btnMonStop, flag=wx.EXPAND|wx.RIGHT)
        vbox1.Add(hbox11, flag=wx.EXPAND|wx.TOP, border=4)

        hbox12 = wx.BoxSizer(wx.HORIZONTAL)
        self.chkMonRaiseW = wx.CheckBox(panel2, wx.ID_ANY, "Raise adxv windows")
        self.chkMonRaiseW.SetValue(1)
        hbox12.Add(self.chkMonRaiseW, flag=wx.EXPAND|wx.RIGHT)
        vbox1.Add(hbox12, flag=wx.EXPAND|wx.TOP, border=4)

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
            self.open_hdf5(1)
    # __init__()

    def on_notebook_page_changed(self, ev):
        self.txtInfo.SetValue("")

        if ev.Selection == 0: # File mode
            self.monitor_timer.Stop()
    # on_notebook_page_changed()

    def on_btnMonRefresh(self, ev):
        self.monitor_timer.Stop()
        self.on_monitor_timer(None)
        self.monitor_timer.Start(5000)
    # on_btnMonRefresh()

    def open_hdf5(self, frameno, raise_window=True):
        self.adxv.open_hdf5(self.h5file, frameno,
                            raise_window=raise_window,
                            binning=int(self.txtBin.GetValue()))
    # open_hdf5()

    def change_frameno(self, frames):
        if all(map(lambda x: x < 1, frames)): frames = [1]
        if all(map(lambda x: x > self.n_images, frames)): frames = [self.n_images]
        frames = filter(lambda x: 0 < x <= self.n_images, frames)
        if not frames: frames = [1]
        self.txtFrame.SetValue(make_range_str(frames))
        self.open_hdf5(frames, raise_window=False)
    # change_frameno()

    def onTextEnter(self, ev):
        frames = parse_range(self.txtFrame.GetValue())
        self.change_frameno(frames)
    # onTextEnter()

    def next_or_back(self, sign):
        s = make_range_str(parse_range(self.txtFrame.GetValue()))
        if "-" not in s: return self.go_rel(sign)
        val = 1
        for x in s.split(","):
            if "-" not in x: continue
            x = map(int, x.split("-"))
            val = max(val, x[1]-x[0]+1)
        return self.go_rel(sign*val)
    # next_or_back()

    def go_rel(self, rel):
        frames = parse_range(self.txtFrame.GetValue())
        frames = map(lambda x: x+rel, frames)
        self.change_frameno(frames)
    # go_rel()

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
        save = self.txtFrame.GetValue()
        self.next_or_back(self.play_relval)
        try: wx.Yield()
        except: pass

        if save == self.txtFrame.GetValue(): self.play_stop()
    # on_play_timer()

    def play_stop(self):
        self.play_timer.Stop()
        self.llbtn.SetLabel("<<")
        self.rrbtn.SetLabel(">>")
    # play_stop()

    def btnInfile_click(self, ev):
        dlg = wx.FileDialog(None, message="Choose a master.h5", defaultFile=self.h5file if self.h5file else "",
                            wildcard="Mater HDF5 (*_master.h5)|*_master.h5")

        if dlg.ShowModal() == wx.ID_OK:
            self.h5file = dlg.GetPath()
            self.read_h5file()
            self.open_hdf5(1)
            
        dlg.Destroy()
    # btnInfile_click()

    def btnLatest_click(self, ev):
        latestlog = os.path.join(os.environ["HOME"], ".bss_latest_file_BL32XU.log")
        if not os.path.isfile(latestlog):
            wx.MessageDialog(None, ".bss_latest_file_BL32XU.log not found in $HOME.",
                             "Error", style=wx.OK).ShowModal()
            return

        latestfile = open(latestlog).read().strip()
        latestfile = os.path.join(os.path.dirname(latestfile), re.sub("_[0-9]+\..*$", "_master.h5", os.path.basename(latestfile)))

        if not os.path.isfile(latestfile):
            wx.MessageDialog(None, "Latest file (%s) not found." % latestfile,
                             "Error", style=wx.OK).ShowModal()
            return

        self.h5file = latestfile
        self.read_h5file()
        self.open_hdf5(1)
    # btnLatest_click ()

    def read_h5file(self):
        self.txtInfile.SetValue(self.h5file)

        try:
            h5 = h5py.File(self.h5file, "r")
            h = eiger_hdf5_interpreter.Interpreter().getRawHeadDict(self.h5file)
            self.n_images = h["Nimages"]

            self.txtInfo.SetValue("""\
         Filename: %s
 Number of frames: %d
             Date: %s
 Phi start, width: %.6f, %.6f deg
  Camera distance: %.2f mm
       Wavelength: %.6f A
    Exposure time: %.6f sec
      Beam center: (%d, %d) px
    """ % (os.path.basename(self.h5file), self.n_images, h["DateStr"].replace("T", " "),
           h["PhiStart"], h["PhiWidth"], h["Distance"]*1.e3,
           h["Wavelength"], h["ExposureTime"], h["BeamX"], h["BeamY"]))
        except:
            self.n_images = 0
            self.txtInfo.SetValue(traceback.format_exc())

        self.change_frameno([1])
    # read_h5file()

    def on_monitor_timer(self, ev):
        print "monitoring"
        response = urllib.urlopen("http://%s/monitor/api/1.6.1/images/monitor" % (self.eiger_host))
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

        if len(data) < 1000:
            data = self.last_monitor_image
            if data is None: return
        else:
            byte = data_size // (3269*3110)
            assert byte == 4 or byte == 2
            data = numpy.fromstring(data[8:-102], dtype=numpy.int32 if byte==4 else numpy.int16).reshape(3269,3110)

            if self.last_monitor_image is not None and (data==self.last_monitor_image).all():
                return

        def request(query):
            r = urllib.urlopen("http://%s/detector/api/1.6.1/%s"%(self.eiger_host, query)).read()
            return json.loads(r)
            
        wavelen = request("config/wavelength")["value"]
        distance = request("config/detector_distance")["value"] * 1000.
        beamx = request("config/beam_center_x")["value"]
        beamy = request("config/beam_center_y")["value"]

        self.txtInfo.SetValue("""\
%s: Downloaded (%d)
    Wavelength: %.5f A
    Detector_distance: %.2f mm
    Beam_xy: %.2f, %.2f
%s"""% (time.ctime(), data_size, wavelen, distance, beamx, beamy, self.txtInfo.GetValue()))

        from yamtbx.dataproc import cbf
        imgfile = os.path.join("/dev/shm", "adxvtmp-%s-%s.cbf"%(getpass.getuser(), os.getpid()))
        cbf.save_numpy_data_as_cbf(data.flatten(), size1=data.shape[1], size2=data.shape[0], title="",
                                   cbfout=imgfile,
                                   pilatus_header="""\
# Detector: Eiger
# Pixel_size 75e-6 m x 75e-6 m
# Wavelength %(Wavelength).6f A
# Detector_distance %(Distance).3e m
# Beam_xy (%(BeamX).1f, %(BeamY).1f) pixels
""" % dict(Wavelength=wavelen, Distance=distance/1.e3, BeamX=beamx, BeamY=beamy))

        self.adxv.open_image(imgfile, raise_window=self.chkMonRaiseW.GetValue())
        self.last_monitor_image = data
    # on_monitor_timer()

def run(argv):
    h5in = argv[0] if argv else None

    app = wx.App()
    app.TopWindow = MainFrame(parent=None, id=wx.ID_ANY, 
                              h5in=h5in)
    app.MainLoop()

if __name__ == "__main__":
    run(os.path.abspath(sys.argv[1:]))
