#!/usr/bin/env yamtbx.python

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import sys
import os
import math
import tempfile
import getpass
import sqlite3
import matplotlib
matplotlib.interactive( True )
matplotlib.use( 'WXAgg' )
import matplotlib.figure
import matplotlib.backends.backend_agg
import matplotlib.backends.backend_wxagg
from matplotlib.ticker import FuncFormatter
from matplotlib.patches import Rectangle, Ellipse
import wx
import wx.lib.newevent
import wx.lib.agw.pybusyinfo
import wx.html
from wx.lib.mixins.listctrl import CheckListCtrlMixin, ListCtrlAutoWidthMixin
import datetime
import time
import glob
import pickle
import collections
import threading
import subprocess
import socket
import xmlrpc.client
import re
import copy
import zmq
from PIL import Image

import numpy
import scipy.spatial

import iotbx.phil
import libtbx.phil
from cctbx.array_family import flex

#from yamtbx.dataproc.myspotfinder import spot_finder_for_grid_scan
from yamtbx.dataproc import bl_logfiles
from yamtbx.dataproc.dataset import re_pref_num_ext
from yamtbx.util import get_number_of_processors, rotate_file
from yamtbx.dataproc.XIO import XIO
from yamtbx.dataproc.myspotfinder import shikalog
from yamtbx.dataproc.myspotfinder import config_manager
from yamtbx.dataproc.XIO.plugins import eiger_hdf5_interpreter

EventResultsUpdated, EVT_RESULTS_UPDATED = wx.lib.newevent.NewEvent()
EventTargetDirChanged, EVT_TARGET_DIR_CHANGED = wx.lib.newevent.NewEvent()
EventScanlogsUpdated, EVT_SCANLOGS_UPDATED = wx.lib.newevent.NewEvent()

gui_phil_str = """\
kuma_addr = None
 .type = str
 .help = "kuma address and port; like 192.168.163.5:1920"

imgview_host = None
 .type = str
 .help = "imgview address; like 192.168.163.5"

ask_directory = False
 .type = bool
 .help = Ask which directory to watch when program started

adxv = None
 .type = path
 .help = adxv command

bl = 32xu 41xu 26b2 44xu 45xu
 .type = choice(multi=False)
 .help = Choose beamline where you start SHIKA

readonly = False
 .type = bool
 .help = If readonly, any files will not be written by a program.

make_html = True
 .type = bool
 .help = make html report file

auto_mode = False
 .type = bool

ring_removal = False
 .type = bool
 .help = Automatically detect rings and set exclude resolution range

mode = zoo
 .type = choice(multi=False)
 .help = When ZOO, use mode=zoo.

dbdir = /isilon/cluster/log/shika/db
 .type = path
 .help = location to write sqlite3 db file.

subport = 5559
 .type = int
 .help = port for ZMQ-SUB to communicate with workers.
pushport = 5556
 .type = int
 .help = port for ZMQ-PUSH to communicate with workers.
"""

class Stat(object):
    def __init__(self):
        self.img_file = None
        self.stats = [] # n_spots, total, mean
        self.spots = []
        self.gonio = None
        self.grid_coord = None
        self.scan_info = None
        self.params = None
        self.thumb_posmag = None
        self.detector = ""

gui_params = None # replaced with phil params
current_stats = collections.OrderedDict()
zmq_context = zmq.Context()
control_send = zmq_context.socket(zmq.PUB)
ventilator_send = zmq_context.socket(zmq.PUSH)

def add_results(results):
    for f, stat in results: current_stats[f] = stat

def read_shika_auto_config(scandir):
    ret = {}
    cfgin = os.path.join(scandir, "shika_auto.config")
    if not os.path.isfile(cfgin): return ret

    for l in open(cfgin):
        l = l.strip()
        if l == "": continue
        if l.startswith("min_score="):
            ret["min_score"] = l[l.index("=")+1:].strip()
        elif l.startswith("min_dist="):
            ret["min_dist"] = l[l.index("=")+1:].strip()
        elif l.startswith("max_hits="):
            ret["max_hits"] = l[l.index("=")+1:].strip()
        else:
            shikalog.warning("Unrecognized config in %s: %s" % (cfgin, l))

    shikalog.info("Read auto-config from %s: %s" % (cfgin, ret))
    return ret
# read_shika_auto_config()

class ReportHTMLMakerThread(object):
    def __init__(self, parent, dont_work=False, make_html=True):
        self.parent = parent
        self.interval = 1
        self.thread = None
        self.queue = []
        self.lock = threading.Lock()
        self.dont_work = dont_work # if True, this thread will actually do nothing.
        self.make_html = make_html

        self.plotFrame = parent.plotFrame
        self.plot_data = None

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
            self.keep_going = False
            self.thread.join()
        else:
            pass

    def is_running(self):
        return self.thread is not None and self.thread.is_alive()

    def run(self):
        while self.keep_going:
            if not self.dont_work and len(self.queue) > 0:
                wdir, rotate = None, None
                # take from queue
                with self.lock:
                    self.queue = list(set(self.queue))
                    wdir, rotate = self.queue.pop(0)

                if self.make_html: self.make(wdir, rotate)
                self.make_dat(wdir)

            if self.interval < 1:
                time.sleep(self.interval)
            else:
                for i in range(int(self.interval/.5)):
                    if self.keep_going:
                        time.sleep(.5)

        self.running = False
    # run()

    """
    Create HTML report.
    XXX Currently, width is fixed (always 600).
    XXX Currently, only hi_pass_resolution_spots are used. If not available (in case XDS, other?), all are used.
    TODO resolve code duplication in prepare_plot()
    TODO Results should be ordered as in diffscan.log
    TODO Don't make the same plot again. (save & check hash of data)
    """
    def prepare_plot(self, f, kind, wdir, rotate):
        def normalize_max(v, maximum=400.):
            max_v = max(v)
            f = maximum / max_v if max_v > 0 else 1.
            return [f*x + 1. for x in v] # add 1 to make zero-value pickable # XXX when max_v is Inf?
        # normalize_max()

        scan_prefix = f[:f.index(" ")] if " (phi=" in f else f
        pngout = os.path.join(wdir, "plot_%s%s.png" % (scan_prefix, kind))
        if rotate:
            rotate_file(pngout)

        xs, ys, ds, imgfs = [], [], [], []
        zero_xs, zero_ys = [], [] # For values of zero
        for imgf, stat in self.plot_data[f]:
            gc = stat.grid_coord
            if gc is None:
                continue
            x, y = gc
            x *= 1000.
            y *= 1000.
            d = stat.stats[("n_spots","total_integrated_signal","median_integrated_signal").index(kind)]
            xs.append(x)
            ys.append(y)
            ds.append(d)
            imgfs.append(imgf)

            if d == 0:
                zero_xs.append(x)
                zero_ys.append(y)

        if len(xs) == 0:
            return "", ""

        win = (max(xs)-min(xs)+1000)/1000*400/80*1.7 # ad-hoc scale
        hin = (max(ys)-min(ys)+1000)/1000*400/80

        fig = matplotlib.figure.Figure(figsize=(win,hin), dpi=80) # figsize in inches
        ax = fig.add_subplot(111)
        p = ax.scatter(xs, ys, s=normalize_max(ds), c=ds, alpha=0.5) # s in points^2
        if max(ds) - min(ds) > 1e-5:
            fig.colorbar(p)
        ax.scatter(zero_xs, zero_ys, s=50, marker="x", c=[0]*len(zero_xs), alpha=0.5)
        ax.set_xlabel("horizontal [um]")
        ax.set_ylabel("vertical [um]")

        scaninfo = self.plot_data[f][0][1].scan_info
        if scaninfo is not None:
            vp, hp = scaninfo.vpoints, scaninfo.hpoints
            vs, hs = scaninfo.vstep*1000., scaninfo.hstep*1000.

            if 1 in (vp, hp) or len(self.plot_data[f]) <= hp:
                ax.set_aspect("auto")
            else:
                ax.set_aspect("equal")

            if vp == hp == 1:
                ax.set_xlim(-10, 10)
                ax.set_ylim(-10, 10)
            elif vp == 1:
                ax.set_xlim(min(xs) - hs, max(xs) + hs)
                ax.set_ylim(-10, 10)
            elif hp == 1:
                ax.set_xlim(-10, 10)
                ax.set_ylim(min(ys) - vs, max(ys) + vs)
            else:
                ax.set_xlim(min(xs) - hs, max(xs) + hs)
                ax.set_ylim(min(ys) - vs, max(ys) + vs)
        else:
            # Should never reach here.. but should we set limit here?
            pass

        canvas = matplotlib.backends.backend_agg.FigureCanvasAgg(fig)
        canvas.print_figure(pngout+".tmp", dpi=80, format="png")
        img_width = fig.get_figwidth() * 80
        img_height = fig.get_figheight() * 80

        map_str = '<map name="%smap">\n' % scan_prefix
        for x, y, imgf in zip(xs, ys, imgfs):
            tx, ty = ax.transData.transform((x,y))
            map_str += '  <area shape="circle" coords="%.2f,%.2f,10" title="%s" onClick=\'plotClick("%s", "%s")\'>\n' % (tx, img_height-ty, os.path.basename(imgf), scan_prefix, os.path.basename(imgf))

        map_str += "</map>"
        return pngout, map_str
    # prepare_plot()

    def make_dat(self, wdir):
        self.plot_data = self.plotFrame.data

        datout = os.path.join(wdir, "summary.gui.dat")
        ofs = open(datout, "w")
        kinds = [rb.GetLabelText() for rb in self.plotFrame.rb_kind]

        print("prefix x y kind data filename", file=ofs)
        for f in self.plot_data:
            for i, kind in enumerate(kinds):
                for imgf, stat in sorted(self.plot_data[f]):
                    gc = stat.grid_coord
                    if gc is None:
                        x, y = "na", "na"
                        shikalog.warning("gc is None! %s"%imgf)
                    else:
                        x, y = gc
                    d = stat.stats[("n_spots","total_integrated_signal","median_integrated_signal").index(kind)]
                    print(f[:f.rindex("(")-1], x, y, kind, d, os.path.basename(imgf), file=ofs)
    # make_dat()

    def make(self, wdir, rotate=False):
        self.plot_data = self.plotFrame.data
        shikalog.info("Making HTML report for %s"%wdir)
        startt = time.time()

        if gui_params.mode == "zoo": htmlout = os.path.join(wdir, "report_zoo.html")
        else: htmlout = os.path.join(wdir, "report.html")

        if rotate: rotate_file(htmlout)

        if gui_params.mode == "zoo": assert len(self.plot_data) <= 1

        kinds = [rb.GetLabelText() for rb in self.plotFrame.rb_kind]
        plots=""
        pngs = []
        for f in self.plot_data:
            scan_prefix = f[:f.index(" ")] if " (phi=" in f else f
            info = self.plot_data[f][0][1].scan_info

            if gui_params.mode == "zoo" and len(self.plot_data[f]) != info.vpoints*info.hpoints:
                continue

            if info is None: info = bl_logfiles.ScanInfo() # Empty info
            plots += '<table border=0 style="margin-bottom:0px">\n  <tr><td>\n'

            if gui_params.mode == "zoo":
                try:
                    im = Image.open(os.path.join(wdir, "../../before.ppm"))
                    im.save(os.path.join(wdir, "loop_before.jpg"))
                except:
                    import traceback
                    print("Can't convert loop image")
                    print(traceback.format_exc())
                plots += '  Loop image</td><td><img src="loop_before.jpg" /></td></tr>\n'
                plots += '  <tr><td>\n'

            plots += '  <table class="info"><tr><th>scan</th><td>%s</td></tr>\n' % scan_prefix
            plots += '    <tr><th>date</th><td>%s</td></tr>\n' % (info.date.strftime("%Y/%m/%d %H:%M:%S") if info.date!=0 else "??")

            if info.is_shutterless():
                plots += '    <tr><th>fixed spindle</th><td>%.2f&deg;</td></tr>\n' % info.fixed_spindle
                plots += '    <tr><th>frame rate</th><td>%.2f [Hz]</td></tr>\n' % info.frame_rate
            else:
                plots += '    <tr><th>osc. start</th><td>%.2f&deg;</td></tr>\n' % info.osc_start
                plots += '    <tr><th>osc. step</th><td>%.2f&deg;</td></tr>\n' % info.osc_step
                plots += '    <tr><th>exp. time</th><td>%.2f [sec]</td></tr>\n' % info.exp_time

            plots += '    <tr><th>beam size</th><td>h= %.1f, v= %.1f [&mu;m]</td></tr>\n' % (info.beam_hsize, info.beam_vsize)
            plots += '    <tr><th>attenuator</th><td>%s %.1f [&mu;m]</td></tr>\n' % info.attenuator
            plots += '    <tr><th>distance</th><td>%.2f [mm]</td></tr>\n' % info.distance
            plots += '    <tr><th>wavelength</th><td>%.4f [&#x212b;]</td></tr>\n' % info.wavelength
            plots += '    <tr><th>scan points</th><td>v=%d, h=%d</td></tr>\n' % (info.vpoints, info.hpoints)
            plots += '    <tr><th>scan steps</th><td>v=%.2f, h=%.2f [&mu;m]</td></tr>\n' % (info.vstep*1000., info.hstep*1000.)
            plots += '  </table>\n'
            
            if gui_params.mode == "zoo":
                dpi = 80.
                win_org, hin_org = self.plotFrame.plotPanel.figure.get_size_inches()
                if win_org < 1: self.plotFrame.plotPanel.figure.set_size_inches(7.5, 6) # needed if plot frame not shown. TODO more appropriate number?
                self.plotFrame.plotPanel.figure.canvas.print_figure(os.path.join(wdir, "%sselected_map.png"%scan_prefix),
                                                                    dpi=int(dpi), format="png")
                plots += '  <td><img name="%s" src="%sselected_map.png" usemap="#%smap" /><br />\n' % (scan_prefix, scan_prefix, scan_prefix)
                plots += '<map name="%smap">\n' % scan_prefix

                win, hin = self.plotFrame.plotPanel.figure.get_size_inches()
                vs, hs = info.vstep*1000., info.hstep*1000.
                for (x, y), imgf in zip(self.plotFrame.plotPanel.plotted_xy, self.plotFrame.plotPanel.current_plotted_imgfs):
                    tx1, ty1 = self.plotFrame.plotPanel.subplot.transData.transform((x-hs/2.,y-vs/2.))
                    tx2, ty2 = self.plotFrame.plotPanel.subplot.transData.transform((x+hs/2.,y+vs/2.))
                    img_height = hin*dpi
                    plots += '  <area shape="rect" coords="%.2f,%.2f,%.2f,%.2f" title="%s" onClick=\'plotClick("%s", "%s")\'>\n' % (tx1, img_height-ty1, tx2, img_height-ty2, os.path.basename(imgf), scan_prefix, os.path.basename(imgf))

                plots += '</map></td></tr></table><br>\n\n'
            else:
                for i, kind in enumerate(kinds):
                    pngout, mapstr = self.prepare_plot(f, kind, wdir, rotate)
                    pngs.append(pngout) # rename later
                    adds = ""
                    if i == 0:
                        plots += '  <td><img name="%s" src="%s" usemap="#%smap" /><br />\n' % (scan_prefix, os.path.basename(pngout), scan_prefix)
                        plots += '<form>\n'
                        adds = ' checked="checked"'
                    plots += '<input type="radio" name="spot_mode" value="%s" onClick="changeplot(this, \'%s\')"%s />%s<br />\n' % (kind, scan_prefix, adds, kind)
                plots += '</form>%s</td></tr></table><br>\n\n' % mapstr # The last mapstr is used. This is dirty way, though.
            plots += '<table border=0 style="margin-bottom:20px">\n  <tr><td>\n'
            plots += '<td style="border:solid 1px #999"><canvas id="%scanvas" width=600 height=600></canvas>\n' % scan_prefix
            plots += '<td id="%sinfo" valign="top"></tr></table>\n\n' % scan_prefix

        result = list(current_stats.items())
        if len(result) == 0:
            shikalog.warning("No results found. Exiting. %s"% wdir)
            return

        dbfile = os.path.join(wdir, "shika.db")
        con = sqlite3.connect(dbfile, timeout=10, isolation_level=None)
        con.execute('pragma query_only = ON;')
        print("Reading data from DB for making report html.")
        c = con.execute("select filename,spots from spots")
        dbspots = dict([(str(x[0]), pickle.loads(str(x[1]))) for x in c.fetchall()])
        spot_data = "var spot_data = {"
        for i, (f, stat) in enumerate(result):
            if stat is None: continue
            bf = os.path.basename(f)
            spots = dbspots[bf]["spots"]
            thumb_posmag = dbspots[bf]["thumb_posmag"]
            r = re.search("^(.*)_([0-9]+)\.[^0-9]+$", bf)
            prefix, num = r.group(1), int(r.group(2))
            spot_data += '"%s":[[' % bf
            for y,x,snr,d in spots:
                #x, y = spot.max_pxl_y(), spot.max_pxl_x()
                pos = thumb_posmag[0:2]
                mag = thumb_posmag[2]
                x, y = (x - pos[0])*mag, (y - pos[1])*mag
                spot_data += "[%d,%d]," % (x, y)

            spot_data += "], %.1f, %.1f, %d, %d]," % (stat.stats[1], stat.stats[2], stat.stats[0], num)

        spot_data += "};"
        spot_data = spot_data.replace("inf,", "Infinity,").replace("nan,", "NaN,")
        
        con.close()

        # Determine img picture extension
        img_ext = ".png" if os.path.exists(os.path.join(wdir, os.path.basename(result[0][0])+".png")) else ".jpg"
        jpg_dirs = "var jpg_dirs = {"
        flag_tiled_jpg = False
        if glob.glob(os.path.join(wdir, "thumb_*")):
            for res in result:
                r = re.search("^(.*)_([0-9]+)\.[^0-9]+$", os.path.basename(res[0]))
                prefix, num = r.group(1), int(r.group(2))
                jd = os.path.join("thumb_%s_%.3d" % (prefix, num//1000))
                if not os.path.exists(jd): flag_tiled_jpg = True  # THIS MAY CAUSE A PROBLEM..
                jpg_dirs += '"%s":"%s",' % (os.path.basename(res[0]), jd)
        else:
            for res in result:
                jpg_dirs += '"%s":".",' % os.path.basename(res[0])

        jpg_dirs += "};"

        ofs = open(htmlout, "w")
        ofs.write("""\
<!DOCTYPE html>
<html>
<head>
  <meta charset="UTF-8" />
  <title>SHIKA report</title>
  <script type="text/javascript">
  <!--
    function changeplot(obj, name){
     document.images[name].src = "plot_"+name+obj.value+".png";
    }
    %(spot_data)s
    %(jpg_dirs)s
""" % dict(spot_data=spot_data,
           jpg_dirs=jpg_dirs if not flag_tiled_jpg else ""))

        if flag_tiled_jpg: # FOR TILED JPEG
            ofs.write("""\
    function plotClick(scanprefix, imgfile) {
        var f = imgfile;
        var data = spot_data[f];
        var img = new Image();
        var idx = Math.floor((data[4]-1)/100);
        var n1 = idx*100+1;
        var n2 = (idx+1)*100;
        img.src = "thumb_" + scanprefix.slice(0,-1) + "/" + scanprefix + ("00000"+n1).slice(-6) + "-" + ("00000"+n2).slice(-6) + ".jpg"; // prefix ends with _

        var idx2 = (data[4]-1)%%100;
        var sx = idx2%%10;
        var sy = Math.floor(idx2/10);

        img.onload = (function(fn){
          return function(){
            var td = document.getElementById(scanprefix+"info");
            td.innerHTML = "<table border=0><tr><td>File name: <td>" + imgfile + "<tr><td>total signal: <td>" + data[1] + "<tr><td>median signal: <td>" + data[2] + "<tr><td>N_spots: <td>" + data[3] + "</table>";

            var t = data[0];
            var canvas = document.getElementById(scanprefix+"canvas");
            var ctx = canvas.getContext('2d');
            ctx.clearRect(0,0,canvas.width,canvas.height);
            ctx.drawImage(this, sx*600, sy*600, 600, 600, 0, 0, 600, 600);
""" % dict(img_ext=img_ext))
            
        else: # FOR SINGLE JPEGs
            ofs.write("""\
    function plotClick(scanprefix, imgfile) {
        var f = imgfile;
        var data = spot_data[f];
        var img = new Image();
        img.src = jpg_dirs[f] + "/" + f + "%(img_ext)s";
        img.onload = (function(fn){
          return function(){
            var td = document.getElementById(scanprefix+"info");
            td.innerHTML = "<table border=0><tr><td>File name: <td>" + imgfile + "<tr><td>total signal: <td>" + data[1] + "<tr><td>median signal: <td>" + data[2] + "<tr><td>N_spots: <td>" + data[3] + "</table>";

            var t = data[0];
            var canvas = document.getElementById(scanprefix+"canvas");
            var ctx = canvas.getContext('2d');
            ctx.clearRect(0,0,canvas.width,canvas.height);
            ctx.drawImage(this, 0, 0);
""" % dict(img_ext=img_ext))

        # Common parts
        ofs.write("""\
            for (var i = 0; i < t.length; i++) {
              ctx.rect(t[i][0]-6, t[i][1]-6, 12, 12);
            }
            ctx.strokeStyle = "red";
            ctx.lineWidth = 1;
            ctx.stroke();

            var center = [300,300];
            ctx.beginPath();
            ctx.strokeStyle = "blue";
            ctx.moveTo(center[0]-10, center[1]);
            ctx.lineTo(center[0]+10, center[1]);
            ctx.moveTo(center[0], center[1]-10);
            ctx.lineTo(center[0], center[1]+10);
            ctx.stroke();
          }
        }(f));
    }
  //-->
  </script>
  <style type="text/css">
  <!--
    table.info {
      border-collapse: separate;
      border-spacing: 7px;
    }
    table.info th {
      text-align: left;
    }

    table.images {
      border-collapse: collapse;
      border: solid 1px #999;
    }
    table.images caption {
      margin-top: 1em;
      text-align: left;
    }
    table.images th,
    table.images td {
      border: solid 1px #999;
    }
    table.images th {
      background: #E6E6E6;
      text-align: center;
      white-space: nowrap;
    }
  -->
  </style>
</head>

<body>
<h1>SHIKA report</h1>
<div align="right">
Created on %(date)s<br>
Original directory: %(wdir)s
</div>
<hr style="height: 1px;border: none;border-top: 1px #000000 dotted;" />

%(plots)s

</body>
</html>

""" % dict(plots=plots,
           date=datetime.datetime.today().strftime("%Y/%m/%d %H:%M:%S"),
           wdir=wdir,
           ))

        shikalog.debug("Renaming png files in %s" % wdir)
        for png in pngs:
            os.rename(png+".tmp", png)

        delt = time.time() - startt
        shikalog.info("HTML making Done (took %f s). Open? firefox %s"% (delt, htmlout))
    # make()
# class ReportHTMLMakerThread

class ConfigManager(object):
    def __init__(self):
        self.items = copy.copy(config_manager.sp_params_strs)
        self.common_params_str = config_manager.get_common_params_str()
    # __init__()

    def get_common_params_str(self): return self.common_params_str
    def set_common_params_str(self, s): self.common_params_str = s
    def get_specific_params_str(self, key): return self.items[key]
    def set_specific_params_str(self, key, s): self.items[key] = s

    def set_item(self, beamline, detector, binning, extra, params_str):
        self.items[(beamline, detector, binning, extra)] = params_str
    # set_item()

    def check_phil_valid(self, phil_str):
        master_params = libtbx.phil.parse(spot_finder_for_grid_scan.master_params_str)
        try:
            working_params, alldef = master_params.fetch(sources=[libtbx.phil.parse(phil_str)],
                                                         track_unused_definitions=True)
            working_params.extract()
            if len(alldef) > 0:
                return "Unknown parameters: " + ", ".join([x.path for x in alldef])
        except RuntimeError as e:
            return e.message
        return ""
    # check_phil_valid()

    def get_names(self):
        ret = []
        for k in self.items:
            s = "%s %s" % (k[0], k[1])
            ex = []
            if k[2] is not None:
                ex.append("%s bin" % k[2])
            if k[3] is not None:
                ex.append(k[3])
            if len(ex) > 0:
                s += " " + ", ".join(ex)

            ret.append((s, k))
        return ret
    # get_names()

    def keys(self): return list(self.items.keys())

    def get_params_by_key(self, key):
        params_str = self.get_common_params_str() + self.get_specific_params_str(key)
        master_params = libtbx.phil.parse(spot_finder_for_grid_scan.master_params_str)
        working_params = master_params.fetch(sources=[libtbx.phil.parse(params_str)])
        return working_params.extract()
    # get_params_by_key()

# class ConfigManager

class CheckListCtrl(wx.ListCtrl, CheckListCtrlMixin, ListCtrlAutoWidthMixin):
    """
    http://zetcode.com/wxpython/advanced/
    """
    def __init__(self, parent, style=wx.LC_REPORT | wx.SUNKEN_BORDER):
        wx.ListCtrl.__init__(self, parent, -1, style=style)
        CheckListCtrlMixin.__init__(self)
        ListCtrlAutoWidthMixin.__init__(self)
    # __init__()
# class CheckListCtrl

class ConfigFrame(wx.Frame):
    class CommonPanel(wx.Panel):
        def __init__(self, parent, manager):
            wx.Panel.__init__(self, parent)
            self.manager = manager
            sizer = wx.GridBagSizer()
            self.SetSizer(sizer)
            lab = wx.StaticText(self, wx.ID_ANY, "Common Settings")
            lab.SetFont(wx.Font(12, wx.FONTFAMILY_DEFAULT,wx.FONTSTYLE_NORMAL,wx.FONTWEIGHT_NORMAL))
            self.txtctrl = wx.TextCtrl(self, style=wx.TE_MULTILINE)
            self.btnRevert = wx.Button(self, wx.ID_ANY, "Revert")
            self.btnApply = wx.Button(self, wx.ID_ANY, "Apply")
            sizer.Add(lab, pos=(0,0), span=(1,2))
            sizer.Add(self.txtctrl, pos=(1,0), span=(1,2), flag=wx.EXPAND|wx.ALL, border=4)
            sizer.Add(self.btnRevert, pos=(2,0), flag=wx.EXPAND)
            sizer.Add(self.btnApply, pos=(2,1), flag=wx.EXPAND)
            sizer.AddGrowableRow(1)
            sizer.AddGrowableCol(0)
            sizer.AddGrowableCol(1)

            self.btnRevert.Bind(wx.EVT_BUTTON, self.btnRevert_onClick)
            self.btnApply.Bind(wx.EVT_BUTTON, self.btnApply_onClick)

            self.btnRevert_onClick(None)
        # __init__()

        def btnRevert_onClick(self, ev): self.txtctrl.SetValue(self.manager.get_common_params_str())

        def btnApply_onClick(self, ev):
            phil_str = self.txtctrl.GetValue()
            err = self.manager.check_phil_valid(phil_str)
            if err == "":
                self.manager.set_common_params_str(phil_str)
                self.GetParent().GetParent().send_control()
            else:
                wx.MessageDialog(None, "Wrong settings! Please resolve following error:\n\n"+err,
                                 "Error", style=wx.OK).ShowModal()
        # btnApply_onClick()
    # class CommonPanel

    class SpecificPanel(wx.Panel):
        def __init__(self, parent, manager, beamline):
            wx.Panel.__init__(self, parent)
            self.manager = manager
            sizer = wx.GridBagSizer()
            self.SetSizer(sizer)
            lab = wx.StaticText(self, wx.ID_ANY, "Specific Settings: ")
            lab.SetFont(wx.Font(12, wx.FONTFAMILY_DEFAULT,wx.FONTSTYLE_NORMAL,wx.FONTWEIGHT_NORMAL))
            self.cmbDet = wx.ComboBox(self, wx.ID_ANY, style=wx.CB_READONLY)
            self.txtctrl = wx.TextCtrl(self, style=wx.TE_MULTILINE)
            self.btnRevert = wx.Button(self, wx.ID_ANY, "Revert")
            self.btnApply = wx.Button(self, wx.ID_ANY, "Apply")
            sizer.Add(lab, pos=(0,0))
            sizer.Add(self.cmbDet, pos=(0,1), flag=wx.EXPAND)
            sizer.Add(self.txtctrl, pos=(1,0), span=(1,2), flag=wx.EXPAND|wx.ALL, border=4)
            sizer.Add(self.btnRevert, pos=(2,0), flag=wx.EXPAND)
            sizer.Add(self.btnApply, pos=(2,1), flag=wx.EXPAND)
            sizer.AddGrowableRow(1)
            sizer.AddGrowableCol(0)
            sizer.AddGrowableCol(1)

            self.btnRevert.Bind(wx.EVT_BUTTON, self.btnRevert_onClick)
            self.btnApply.Bind(wx.EVT_BUTTON, self.btnApply_onClick)
            self.cmbDet.Bind(wx.EVT_COMBOBOX, self.btnRevert_onClick) # Just reverting works.

            self.set_names(beamline)
            self.btnRevert_onClick(None)
        # __init__()

        def set_names(self, beamline=None):
            self.keys = {}
            self.cmbDet.Clear()

            for name, key in self.manager.get_names():
                self.cmbDet.Append(name)
                self.keys[name] = key

            self.cmbDet.Select(0)

            if beamline == "32xu":
                fltr = [x for x in enumerate(self.cmbDet.GetItems()) if "BL32XU" in x[1]]
                if len(fltr) > 0:
                    self.cmbDet.Select(fltr[0][0])
            elif beamline == "41xu":
                fltr = [x for x in enumerate(self.cmbDet.GetItems()) if "BL41XU" in x[1]]
                if len(fltr) > 0:
                    self.cmbDet.Select(fltr[0][0])
            elif beamline == "26b2":
                fltr = [x for x in enumerate(self.cmbDet.GetItems()) if "BL26B2" in x[1]]
                if len(fltr) > 0:
                    self.cmbDet.Select(fltr[0][0])
            else:
                shikalog.warning("Unknown beamline: %s" %beamline)
        # set_names()

        def btnRevert_onClick(self, ev):
            key = self.keys[self.cmbDet.GetValue()]
            self.txtctrl.SetValue(self.manager.get_specific_params_str(key))
        # btnRevert_onClick()

        def btnApply_onClick(self, ev):
            key = self.keys[self.cmbDet.GetValue()]

            phil_str = self.txtctrl.GetValue()
            err = self.manager.check_phil_valid(phil_str)
            if err == "":
                self.manager.set_specific_params_str(key, phil_str)
                self.GetParent().GetParent().send_control()
            else:
                wx.MessageDialog(None, "Wrong settings! Please resolve following error:\n\n"+err,
                                 "Error", style=wx.OK).ShowModal()

        # btnApply_onClick()

    # class SpecificPanel

    def __init__(self, parent=None, beamline=None):
        wx.Frame.__init__(self, parent=parent, id=wx.ID_ANY, title="Settings",
                          size=(800,600))

        self.manager = ConfigManager()
        self.splitter = wx.SplitterWindow(self, id=wx.ID_ANY)
        self.splitter.SetSashGravity(0.5)
        self.panel1 = self.CommonPanel(self.splitter, self.manager)
        self.panel2 = self.SpecificPanel(self.splitter, self.manager, beamline)
        self.splitter.SplitVertically(self.panel1, self.panel2)
        self.splitter.SetSashPosition(400)
        self.panel2.SetFocus() # for key capture

        self.Bind(wx.EVT_CLOSE, lambda e: self.Hide()) # Don't destroy this frame when closed
        self.Bind(wx.EVT_KEY_UP, self.OnKeyUp)
    # __init__()

    def OnKeyUp(self, event):
        if event.GetKeyCode() == wx.WXK_ESCAPE:
            self.panel1.btnRevert_onClick(None)
            self.panel2.btnRevert_onClick(None)
            self.Hide()
    # OnKeyUp()

    def send_control(self):
        params = {}
        for key in list(self.manager.keys()): params[key] = self.manager.get_params_by_key(key)
        control_send.send_pyobj(dict(params=params))
    # send_control()

# class ConfigFrame

class ControlPanel(wx.Panel):
    def __init__(self, mainFrame, params, parent=None, id=wx.ID_ANY):
        wx.Panel.__init__(self, parent=parent, id=id)
        self.mainFrame = mainFrame

        self.current_target_dir = None
        self.current_target_fpref = None

        self.vbox = vbox = wx.BoxSizer(wx.VERTICAL)

        self.treectrl = wx.TreeCtrl(self, size=(500, 450))
        self.il_for_treectrl = wx.ImageList(16,16)
        self.il_for_treectrl.Add(wx.ArtProvider.GetBitmap(wx.ART_FOLDER, wx.ART_OTHER, (16,16)))
        self.il_for_treectrl.Add(wx.ArtProvider.GetBitmap(wx.ART_FILE_OPEN, wx.ART_OTHER, (16,16)))
        self.il_for_treectrl.Add(wx.ArtProvider.GetBitmap(wx.ART_NEW_DIR, wx.ART_OTHER, (16,16)))
        self.il_for_treectrl.Add(wx.ArtProvider.GetBitmap(wx.ART_NORMAL_FILE, wx.ART_OTHER, (16,16)))
        self.treectrl.AssignImageList(self.il_for_treectrl)
        self.dic_for_tree = {}
        vbox.Add(self.treectrl, flag=wx.EXPAND|wx.TOP, border=4)

        self.configFrame = ConfigFrame(self, beamline=params.bl)

        self.btnConfig = wx.Button(self, wx.ID_ANY, "Edit settings")
        self.btnConfig.Disable()
        self.btnUpdate = wx.Button(self, wx.ID_ANY, "Recalculate result (if you need)")
        self.btnUpdate.Disable()
        self.btnShowPlot = wx.Button(self, wx.ID_ANY, "Show plot")
        self.btnShowPlot.SetFont(wx.Font(10, wx.FONTFAMILY_DEFAULT,
                                         wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_BOLD))
        self.btnSetExResRange = wx.Button(self, wx.ID_ANY, "Set exclude_resolution ranges (special)")
        self.user_defined_exclude_resolution_ranges = [] # x100 resolution range
        self.detected_exclude_resolution_ranges = [] # x100 resolution range

        hbox0 = wx.BoxSizer(wx.HORIZONTAL)
        hbox0.Add(wx.StaticText(self, wx.ID_ANY, "TopDir: "), flag=wx.LEFT|wx.ALIGN_CENTER_VERTICAL)
        self.txtTopDir = wx.TextCtrl(self, wx.ID_ANY, size=(350,25))
        self.txtTopDir.SetEditable(False)
        self.btnUpdateDir = wx.Button(self, wx.ID_ANY, "Update tree")
        hbox0.Add(self.txtTopDir, flag=wx.EXPAND|wx.RIGHT)
        hbox0.Add(self.btnUpdateDir)

        self.grpTarget = wx.StaticBox(self, wx.ID_ANY, "Target")
        self.vbox_grpTarget = wx.StaticBoxSizer(self.grpTarget, wx.VERTICAL)
        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        hbox1.Add(wx.StaticText(self, wx.ID_ANY, "Dir: "), flag=wx.LEFT|wx.ALIGN_CENTER_VERTICAL)
        self.cmbTargetDir = wx.ComboBox(self, wx.ID_ANY, size=(350,25), style=wx.CB_DROPDOWN)
        #self.btnTargetDir = wx.Button(self, wx.ID_ANY, "...", size=(25,25))
        self.chkTargetDir = wx.CheckBox(self, wx.ID_ANY, "Autofind")
        hbox1.Add(self.cmbTargetDir, flag=wx.EXPAND|wx.LEFT)
        #hbox1.Add(self.btnTargetDir, flag=wx.EXPAND|wx.LEFT)
        hbox1.Add(self.chkTargetDir, flag=wx.EXPAND|wx.LEFT)
        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        self.chkTrackLatest = wx.CheckBox(self, wx.ID_ANY, "Track the latest result (auto-select scan)")
        self.btnReload = wx.Button(self, wx.ID_ANY, "Reload")
        hbox2.Add(self.btnReload, flag=wx.EXPAND|wx.LEFT)
        hbox2.Add(self.chkTrackLatest, flag=wx.EXPAND|wx.LEFT, border=10)
        self.vbox_grpTarget.Add(hbox1)
        self.vbox_grpTarget.Add(hbox2)

        vbox.Add(self.btnConfig, flag=wx.EXPAND|wx.TOP, border=4)
        vbox.Add(self.btnUpdate, flag=wx.EXPAND|wx.TOP, border=4)
        vbox.Add(self.btnSetExResRange, flag=wx.EXPAND|wx.TOP, border=4)
        vbox.Add(self.btnShowPlot, flag=wx.EXPAND|wx.TOP, border=4)
        vbox.Add(hbox0, flag=wx.EXPAND|wx.TOP, border=4)
        vbox.Add(self.vbox_grpTarget, flag=wx.EXPAND|wx.TOP, border=4)

        self.btnConfig.Bind(wx.EVT_BUTTON, self.btnConfig_onClick)
        self.btnUpdate.Bind(wx.EVT_BUTTON, self.btnUpdate_onClick)
        self.btnUpdateDir.Bind(wx.EVT_BUTTON, self.btnUpdateDir_onClick)
        self.btnShowPlot.Bind(wx.EVT_BUTTON, self.btnShowPlot_click)
        self.btnSetExResRange.Bind(wx.EVT_BUTTON, self.btnSetExResRange_onClick)
        self.btnReload.Bind(wx.EVT_BUTTON, lambda e: self.mainFrame.load_results())
        #self.btnTargetDir.Bind(wx.EVT_BUTTON, self.btnTargetDir_click)
        self.chkTrackLatest.Bind(wx.EVT_CHECKBOX, self.chkTrackLatest_onCheck)
        self.chkTargetDir.Bind(wx.EVT_CHECKBOX, self.chkTargetDir_onCheck)
        self.cmbTargetDir.Bind(wx.EVT_COMBOBOX, self.cmbTargetDir_onSelect)
        self.treectrl.Bind(wx.EVT_TREE_SEL_CHANGED, self.treectrl_onSelChanged)
        # Radio button to toggle displayed spots
        self.vbox_rb = vbox1 = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(vbox1, flag=wx.EXPAND, border=4)
        self.rbuttons = []

        self.SetSizer(vbox)

        self.result_update_timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.on_result_update_timer, self.result_update_timer)
        self.db_last_latest_time = None
        self.result_update_timer.Start(5000)

        self.Bind(EVT_RESULTS_UPDATED, self.onResultsUpdated)
        self.Bind(EVT_TARGET_DIR_CHANGED, self.onTargetDirChanged)
        self.Bind(EVT_SCANLOGS_UPDATED, self.onScanlogsUpdated)

        # Default behavior
        self.chkTargetDir.SetValue(True)
        self.chkTrackLatest.SetValue(True)
        ev = wx.CommandEvent(wx.wxEVT_COMMAND_CHECKBOX_CLICKED)
        ev.SetInt(1)
        wx.PostEvent(self.chkTargetDir, ev)
        wx.PostEvent(self.chkTrackLatest, ev)

        if gui_params.mode == "zoo":
            # User is not allowed to touch this checkbox when Zoo mode
            self.chkTrackLatest.Disable()

    # __init__()

    def btnTargetDir_click(self, event):
        current_dir = os.path.join(self.txtTopDir.GetValue(), self.cmbTargetDir.GetValue())
        dlg = wx.DirDialog(None, message="Choose a directory to watch", defaultPath=current_dir)

        if dlg.ShowModal() == wx.ID_OK:
            dirsel = dlg.GetPath()
            if os.path.isdir(dirsel):
                pass
            else: # never reaches, maybe. it seems to automatically create new empty directory.. bleugh.
                wx.MessageDialog(None, "Directory does not exist!", "Error", style=wx.OK).ShowModal()
        dlg.Destroy()
    # btnTargetDir_click()

    def cmbTargetDir_onSelect(self, event):
        target_dir = os.path.join(self.txtTopDir.GetValue(), self.cmbTargetDir.GetValue())
        wx.PostEvent(self, EventTargetDirChanged(target=target_dir,fpref=None))
    # cmbTargetDir_onSelect()

    def get_selected_dir_fpref(self):
        item = self.treectrl.GetSelection()
        sel = self.treectrl.GetPyData(item)
        seldir = os.path.dirname(sel) if " (phi=" in sel else sel
        fpref = os.path.basename(sel) if " (phi=" in sel else None
        return seldir, fpref
    # get_selected_dir_fpref()

    def treectrl_onSelChanged(self, event):
        seldir, fpref = self.get_selected_dir_fpref()
        target_dir = os.path.join(self.txtTopDir.GetValue(), seldir)
        print("DEBUG::target_dir, fpref=", target_dir, fpref)
        print("DEBUG::current           ", self.current_target_dir, self.current_target_fpref)
        if target_dir != self.current_target_dir or fpref != self.current_target_fpref:
            wx.PostEvent(self, EventTargetDirChanged(target=target_dir,
                                                     fpref=fpref))
    # treectrl_onSelChanged()

    def onTargetDirFileChanged(self, fpref):
        """
        Called when scan prefix was changed in tree ctrl.
        Calls plotFrame.rb_clicked() to change plot.
        """
        self.mainFrame.plotFrame.rb_clicked(None)
        self.current_target_fpref = fpref
    # onTargetDirFileChanged()

    def onTargetDirChanged(self, event, need_load=True):
        """
        Called via EVT_TARGET_DIR_CHANGED. This event is invoked by
         - cmbTargetDir_onSelect()
         - treectrl_onSelChanged() 
         - on_result_update_timer() if self.chkTargetDir is checked
         - MainFrame.__init__()

        event.target is the new target directory. It must be an absolute path.
        """

        new_target = os.path.abspath(event.target) # to remove /. or /../hoge
        if new_target == self.current_target_dir:
            if event.fpref is not None:
                self.onTargetDirFileChanged(event.fpref)
            return

        # Clear shown results. TODO this part need to move to other class.
        mainframe = self.mainFrame
        mainframe.data = collections.OrderedDict()

        mainframe.plotFrame.data = collections.OrderedDict()
        mainframe.plotFrame.plotPanel.reset()
        mainframe.plotFrame.splotFrame.reset()
        mainframe.grid.clear()

        print("DEBUG::self.current_target_dir, new_target=", self.current_target_dir, new_target)
        self.current_target_dir = new_target
        self.current_target_fpref = None

        scanlog = os.path.join(new_target, "diffscan.log")
        if not os.path.isfile(scanlog):
            shikalog.warning("NOT FOUND: %s"% scanlog)
            return

        # Read auto-config
        cfg = read_shika_auto_config(new_target)
        if "min_dist" in cfg: mainframe.plotFrame.peakPanel.txtMinDist.SetValue(cfg["min_dist"])
        if "min_score" in cfg: mainframe.plotFrame.peakPanel.txtMinScore.SetValue(cfg["min_score"])
        if "max_hits" in cfg: mainframe.plotFrame.peakPanel.txtMaxHits.SetValue(cfg["max_hits"])

        self.clear_detected_exclude_resolution_ranges()

        # Update combo box
        cmb_items = self.cmbTargetDir.GetItems()
        new_path = os.path.relpath(new_target, self.txtTopDir.GetValue())
        if new_path in cmb_items:
            self.cmbTargetDir.Select(cmb_items.index(new_path))
        else:
            self.cmbTargetDir.Append(new_path)
            self.cmbTargetDir.Select(self.cmbTargetDir.GetCount()-1)

        # Stop timer (and will restart again if running)
        self.result_update_timer.Stop()

        if event.fpref is not None: # If coming from other directory, uncheck it (not to go to latest).
            if mainframe.ctrlFrame.chkTrackLatest.IsEnabled():
                mainframe.ctrlFrame.chkTrackLatest.SetValue(False)

        # Select item in tree (need this - change selected directory - before loading data)
        if event.fpref is None: # if fpref exists, it means user clicked a file.
            k = tuple(os.path.relpath(new_target, self.mainFrame.topdir).split(os.sep))
            if k in self.dic_for_tree:
                self.treectrl.EnsureVisible(self.dic_for_tree[k])
                self.treectrl.SelectItem(self.dic_for_tree[k])
                self.treectrl.Expand(self.dic_for_tree[k])
        else:
            self.onTargetDirFileChanged(event.fpref)

        # Load .pkl data
        if need_load:
            mainframe.load_results()

            # After loading data, do this to select some child
            #if event.fpref is None:
            #    mainframe.track_latest_result() # Because folder-icon was clicked. #XXX OK? <- NO!!! But if no diffscan.log, we can do this..

        self.result_update_timer.Start()
    # onTargetDirChanged()

    def onScanlogsUpdated(self, event):
        # Update directory tree
        dirs = [os.path.relpath(os.path.dirname(x[0]), self.mainFrame.topdir) for x in event.scanlogs]
        dic = self.dic_for_tree
        #print "DEBUG:: dirs=", dirs
        #print "DEBUG:: dic=", dic
        for d in dirs:
            sp = d.split(os.sep)
            for i in range(len(sp)):
                key, keypar = tuple(sp[:i+1]), tuple(sp[:i])
                print("  DEBUG:: key, keypar=", key, keypar)
                if key not in dic:
                    dic[key] = self.treectrl.AppendItem(dic[keypar], sp[i], image=0)
                    self.treectrl.EnsureVisible(dic[key])
                    self.treectrl.Expand(dic[key])
                    self.treectrl.SetPyData(dic[key], os.sep.join(sp[:i+1]))

        if self.chkTargetDir.GetValue():
            pass
    # onScanlogsUpdated()

    def btnShowPlot_click(self, event):
        self.mainFrame.plotFrame.Show()
        self.mainFrame.plotFrame.Raise()
    # btnShowPlot_click()

    def btnConfig_onClick(self, event):
        self.configFrame.Show()
        self.configFrame.Raise()
    # btnConfig_onClick()

    def chkTargetDir_onCheck(self, event):
        if event.IsChecked():
            shikalog.info("Scanlog watch timer start")
    # chkTargetDir_onCheck()

    def chkTrackLatest_onCheck(self, event):
        if event.IsChecked():
            self.mainFrame.track_latest_result()
    # chkTrackLatest_onCheck()

    def onResultsUpdated(self, ev):
        result = ev.result

        # When target directory is changed before the spotfinder is finished..
        if len(ev.result) > 0 and self.current_target_dir != os.path.dirname(ev.result[0][0]):
            shikalog.error("Mismatch!! %s %s" % (self.current_target_dir, ev.result[0][0]))
            return

        if len(result) == 0:
            return

        startt = time.time()
        d = wx.lib.agw.pybusyinfo.PyBusyInfo("Updating %d results.." % len(result), title="Busy SHIKA")
        try:
            try: wx.SafeYield()
            except: pass

            for f, stat in sorted(result):
                #if gui_params.ring_removal:
                #    self.add_detected_exclude_resolution_ranges(stat.ring_res_ranges, updategui=False)

                if f not in self.mainFrame.data:
                    item = MainFrame.Item(f)
                    self.mainFrame.data[f] = item

            #self.diffscan_manager.add_results(result) # Must be called before this function!
            #self.diffscan_manager.update_scanlogs() # Must be called outside!!
            self.update_rbuttons()
            self.mainFrame.update_result(append=True)
        finally:
            d = None

        shikalog.info("Updating took %f s. len(data)= %d, len(result)= %d." % (time.time() - startt, len(self.mainFrame.data), len(result)))
    #  onResultsUpdated()

    def on_result_update_timer(self, ev):
        dbfile = os.path.join(gui_params.dbdir, "%s.db"%getpass.getuser())
        if not os.path.isfile(dbfile): return
        topdir = str(self.txtTopDir.GetValue())

        con = sqlite3.connect(dbfile, timeout=10, isolation_level=None)
        con.execute('pragma query_only = ON;')
        cur = con.cursor()

        if self.db_last_latest_time is None:
            c = cur.execute("select * from updates where dirname like ?", ("%s%%"%topdir,))
        else:
            c = cur.execute("select * from updates where dirname like ? and time>?",
                            ("%s%%"%topdir, self.db_last_latest_time-5)) # 5 sec buffer for possible delay of nfs

        self.db_last_latest_time = time.time()

        updates = [(os.path.dirname(x[0]),x[1]) for x in c.fetchall()] # take upper directory
        if not updates: return

        updates.sort(key=lambda x:x[1])
        #self.db_last_latest_time = updates[-1][1] # don't save this; save current time instead (above)

        scanlogs = []
        for root, t in updates: # root is */_spotfinder
            scanlog = os.path.join(root, "diffscan.log")
            if os.path.isfile(scanlog):
                scanlogs.append((scanlog, os.path.getmtime(scanlog)))
                
        wx.PostEvent(self, EventScanlogsUpdated(scanlogs=scanlogs))

        if self.chkTargetDir.GetValue():
            wx.PostEvent(self, EventTargetDirChanged(target=updates[-1][0],fpref=None))

        if os.path.normpath(updates[-1][0]) == os.path.normpath(self.current_target_dir):
            self.mainFrame.load_results()
    # on_result_update_timer()

    def get_spot_draw_mode(self):
        for rbtn in self.rbuttons:
            if rbtn.GetValue():
                return rbtn.GetLabelText()
        return None
    # get_spot_draw_mode()

    def set_spot_draw_mode(self, mode):
        """
        Try to set mode (if found)
        """
        for rbtn in self.rbuttons:
            if rbtn.GetLabelText() == mode:
                rbtn.SetValue(True)
                return True
        return False
    # set_spot_draw_mode()

    def rb_clicked(self, event, call_from_runbutton=False, append=False):
        mode = self.get_spot_draw_mode()
        if mode is None:
            return

        self.mainFrame.grid.refresh_image() # To update spot drawing
        self.mainFrame.plotFrame.rb_clicked(None)

        if self.mainFrame.grid.current_img_file is None:
            self.mainFrame.grid.load(list(self.mainFrame.data.keys())[0])
        else:
            self.mainFrame.grid.update()
    # rb_clicked()

    def update_rbuttons(self):
        result = current_stats
        if len(result) < 1:
            return

        self.rbuttons = []
        self.vbox_rb.DeleteWindows()
        #for i, k in enumerate(spotfinder_info.all_keys):
        #    if i == 0:
        #        self.rbuttons.append(wx.RadioButton(self, wx.ID_ANY, k, style=wx.RB_GROUP))
        #    else:
        #        self.rbuttons.append(wx.RadioButton(self, wx.ID_ANY, k))
        self.rbuttons.append(wx.RadioButton(self, wx.ID_ANY, "show all spots", style=wx.RB_GROUP))
        self.rbuttons.append(wx.RadioButton(self, wx.ID_ANY, "do not show spots"))

        for rb in self.rbuttons:
            self.vbox_rb.Add(rb)
            rb.Bind(wx.EVT_RADIOBUTTON, self.rb_clicked)

        self.Fit()
    # update_rbuttons()

    def btnUpdate_onClick(self, event):
        if len(self.mainFrame.data) == 0:
            shikalog.debug("Recalculation does not make sense because no data.")
            return

        shikalog.debug("Recalculation button pressed.")

        if wx.MessageDialog(None, "All results in this directory will be recalculated. Are you sure?",
                            "Confirm", style=wx.YES_NO).ShowModal() == wx.ID_NO:
            shikalog.debug("Recalculation canceled.")
            return


        dbfile = os.path.join(os.path.dirname(list(self.mainFrame.data.keys())[0]), "_spotfinder", "shika.db")
        if os.path.isfile(dbfile):
            con = sqlite3.connect(dbfile, timeout=10)
            con.execute("delete from spots")
            con.execute("delete from stats")
            con.execute("delete from status")
            con.commit()

        headers = {}

        h5files = []
        for imgfile in list(self.mainFrame.data.keys()):
            frameno = int(re.search(".*_([0-9]*)\..*$", imgfile).group(1))
            if os.path.isfile(imgfile):
                ventilator_send.send_json(dict(imgfile=imgfile, idx=frameno))
            else:
                h5master = re.sub("_[0-9]*\.img$", "_master.h5", imgfile) # XXX if binned?
                if not os.path.isfile(h5master): continue
                h5files.append(h5master)

        # All frames in h5 master!
        for h5master in set(h5files):
            if h5master not in headers:
                h = eiger_hdf5_interpreter.Interpreter().getRawHeadDict(h5master)
                headers[h5master] = dict(beam_center_x=float(h["BeamX"]),  # numpy data cannot be serialized..
                                         beam_center_y=float(h["BeamY"]),
                                         distance=float(h["Distance"]),
                                         wavelength=float(h["Wavelength"]),
                                         pixel_size_x=float(h["PixelX"]),
                                         file_prefix=re.sub("_master.h5$","", os.path.basename(h5master)))

            for i in range(h["Nimages"]):
                headers[h5master]["frame"] = i
                imgfile = h5master.replace("_master.h5", "_%.6d.img"%(i+1))
                print("Sending", h5master, i)
                ventilator_send.send_json(dict(imgfile=imgfile, h5master=h5master, idx=i+1, header=headers[h5master]))

            

    # btnUpdate_onClick()

    def btnUpdateDir_onClick(self, event):
        """Update directory tree by os.walk()"""

        d = wx.lib.agw.pybusyinfo.PyBusyInfo("Finding subdirectories..", title="Busy SHIKA")

        try: wx.SafeYield()
        except: pass

        scanlogs = []
        for root, dirnames, filenames in os.walk(str(self.txtTopDir.GetValue())):
            if "diffscan.log" in filenames:
                scanlog = os.path.join(root, "diffscan.log")
                scanlogs.append((scanlog, os.path.getmtime(scanlog)))

        wx.PostEvent(self, EventScanlogsUpdated(scanlogs=scanlogs))
        d = None
    # btnUpdateDir_onClick()

    def btnSetExResRange_onClick(self, event):
        class Dlg(wx.Dialog):
            def __init__(self, parent, ranges):
                wx.Dialog.__init__(self, parent, wx.ID_ANY, "Exclude Resolution Ranges", size=(250, 150))

                vbox = wx.BoxSizer(wx.VERTICAL)
                self.txtComment = wx.TextCtrl(self, wx.ID_ANY, "", (95, 155), style=wx.TE_MULTILINE)
                hbox = wx.BoxSizer(wx.HORIZONTAL)
                btnOK = wx.Button(self, wx.ID_ANY, 'OK', size=(70, 30))
                btnCancel = wx.Button(self, wx.ID_ANY, 'Cancel', size=(70, 30))
                hbox.Add(btnOK, 1)
                hbox.Add(btnCancel, 1, wx.LEFT, 5)

                vbox.Add(wx.StaticText(self, wx.ID_ANY, "Exclude resolution ranges:"))
                vbox.Add(wx.StaticText(self, wx.ID_ANY, " e.g. 12.0 14.0"))
                vbox.Add(self.txtComment, 1, wx.GROW|wx.LEFT)
                vbox.Add(hbox, 1, wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, 10)
                #vbox.Add(wx.StaticText(self, wx.ID_ANY, "Note that this configuration\n does not affect html report."))
                self.SetSizer(vbox)

                self.txtComment.SetFocus()
                self.ranges = ranges
                self.txtComment.SetValue("\n".join(["%.2f %.2f" % (r[0]/100.,r[1]/100.) for r in self.ranges])+"\n")

                btnOK.Bind(wx.EVT_BUTTON, self.btnOK_click)
                btnCancel.Bind(wx.EVT_BUTTON, lambda e: self.Destroy())
            # __init__()

            def btnOK_click(self, event):
                try:
                    newranges = []
                    for l in self.txtComment.GetValue().splitlines():
                        if l.strip() == "":
                            continue
                        sp = list(map(float, l.strip().split()))
                        if len(sp) != 2:
                            raise Exception("More than or less than two numbers in line.")
                        if abs(sp[0] - sp[1]) < 1e-4:
                            raise Exception("Idential values.")
                        newranges.append((sp[0], sp[1]))
                except:
                    wx.MessageDialog(None, "Check inputs.\nSpecify two different real numbers in line.",
                                     "Error", style=wx.OK).ShowModal()
                    return

                del self.ranges[:]
                self.ranges.extend(newranges)
                self.Destroy()
            # btnOK_click()
        # class Dlg

        #exc_range = self.get_exclude_resolution_ranges() #copy.copy(self.user_defined_exclude_resolution_ranges)
        exc_range = copy.copy(self.user_defined_exclude_resolution_ranges)
        dlg = Dlg(self, exc_range)
        dlg.ShowModal()

        if exc_range != self.get_exclude_resolution_ranges():
            self.update_user_defined_exclude_resolution_ranges(exc_range)
    # btnSetExResRange_onClick()

    def add_detected_exclude_resolution_ranges(self, rr, updategui=True):
        rrr = set([(int(x[0]*100.), int(x[1]*100.)) for x in rr])
        orgset = set(self.detected_exclude_resolution_ranges)
        if not rrr.issubset(orgset):
            toadd = list(rrr.difference(orgset))
            #self.detected_exclude_resolution_ranges.extend(

            for r in toadd:
                flag_add = True
                for j, r2 in enumerate(self.detected_exclude_resolution_ranges):
                    if min(r2) <= min(r) and max(r) <= max(r2):
                        flag_add = False
                        break # skip this (redundant)
                    elif min(r) <= min(r2) and max(r2) <= max(r):
                        self.detected_exclude_resolution_ranges[j] = r
                        flag_add = False
                        break
                    elif min(r) <= min(r2) and max(r) <= max(r2):
                        self.detected_exclude_resolution_ranges[j] = (max(r2), min(r))
                        flag_add = False
                        break
                    elif min(r) <= max(r2) and max(r2) <= max(r):
                        self.detected_exclude_resolution_ranges[j] = (max(r), min(r2))
                        flag_add = False
                        break

                if flag_add:
                    self.detected_exclude_resolution_ranges.append(r)

            if orgset != set(self.detected_exclude_resolution_ranges):
                shikalog.info("detected exclude_resolution_range = %s" % [(x[0]/100.,x[1]/100.) for x in self.detected_exclude_resolution_ranges])

            if updategui: # XXX do this if range was actually updated
                self.update_rbuttons()
                self.mainFrame.update_result()
    # add_detected_resolution_ranges()

    def update_user_defined_exclude_resolution_ranges(self, rr):
        rr = set([(int(x[0]*100.), int(x[1]*100.)) for x in rr])
        if set(self.user_defined_exclude_resolution_ranges + self.detected_exclude_resolution_ranges) == rr:
            shikalog.debug("user_defined_exclude_resolution_ranges not changed")
        else:
            self.user_defined_exclude_resolution_ranges = list(rr)
            shikalog.info("user_defined_exclude_resolution_ranges = %s" % [(x[0]/100., x[1]/100.) for x in rr])
            self.mainFrame.load_results()

        self.detected_exclude_resolution_ranges = [] # user's definition overrides auto-detected ranges
    # update_user_defined_exclude_resolution_ranges()

    def clear_detected_exclude_resolution_ranges(self):
        self.detected_exclude_resolution_ranges = []
    # clear_detected_exclude_resolution_ranges()

    def get_exclude_resolution_ranges(self, stat=None):
        #rr = self.user_defined_exclude_resolution_ranges + self.detected_exclude_resolution_ranges
        rr =  [(x[0]/100., x[1]/100.) for x in self.user_defined_exclude_resolution_ranges]
        if stat is not None and hasattr(stat, "ring_res_ranges"):
            rr += stat.ring_res_ranges
        return rr #map(lambda x:(x[0]/100., x[1]/100.), rr)
    # get_exclude_resolution_ranges()
# class ControlPanel

class ScatterPlotFrame(wx.Frame):
    class ScatterPlotPanel(wx.Panel):
        def __init__(self, parent):
            wx.Panel.__init__(self, parent)#, size=(600,400))
            self.sizer = wx.BoxSizer(wx.VERTICAL)
            self.SetSizer(self.sizer)

            self.figure = matplotlib.figure.Figure(dpi=80)
            self.subplot = self.figure.add_subplot(111)
            self.points = []

            self.canvas = matplotlib.backends.backend_wxagg.FigureCanvasWxAgg(self, wx.ID_ANY, self.figure)
            self.sizer.Add(self.canvas, 1, flag=wx.LEFT|wx.TOP|wx.GROW)

            self.reset()
        # __init__()

        def reset(self):
            for p in self.points:
                p.remove()

            self.points = []
            self.SetSize((self.Size[0],self.Size[1])) # Clear drawn plot
        # reset()

    # class ScatterPlotPanel

    def __init__(self, parent=None, id=wx.ID_ANY):
        wx.Frame.__init__(self, parent=parent, id=id, title="Plot",
                          size=(600,600))
        self.Bind(wx.EVT_CLOSE, lambda e: self.Hide()) # Don't destroy this frame when closed
        self.splotPanel = self.ScatterPlotPanel(self) # scatter plot (I vs d^-2)
        self.reset = self.splotPanel.reset
        self.statusbar = self.CreateStatusBar()
        self.splotPanel.canvas.mpl_connect('motion_notify_event', self.canvas_onMouseMove)
    # __init__()

    def plot(self, spots, mode, res_outer, filename=None):
        s2_formatter = lambda x,pos: "inf" if x == 0 else "%.2f" % (1./math.sqrt(x))
        log_formatter = lambda x,pos: "%.2e" % (10.**x)

        xs, ys = [], []
        for y,x,i,d in spots:
            s2 = 1./d**2 if d > 0 else -1
            xs.append(s2)
            ys.append(math.log10(i))

        self.splotPanel.subplot.xaxis.set_major_formatter(FuncFormatter(s2_formatter))
        self.splotPanel.subplot.yaxis.set_major_formatter(FuncFormatter(log_formatter)) # set_yscale("log") didn't work.. why?
        p = self.splotPanel.subplot.scatter(xs, ys)
        self.splotPanel.points = [p]
        self.splotPanel.subplot.set_xlabel("resolution (s^2)")
        self.splotPanel.subplot.set_ylabel("intensity")
        if res_outer is not None:
            self.splotPanel.subplot.set_xlim(0, 1./res_outer**2)

        if filename is not None:
            self.SetTitle(os.path.basename(filename))
    # plot()

    def canvas_onMouseMove(self, event):
        if None not in (event.xdata, event.ydata):
            d = math.sqrt(1./event.xdata) if event.xdata > 0 else float("inf")
            self.statusbar.SetStatusText("d= %.2f, I= %.2f" % (d, 10**event.ydata))
        else:
            self.statusbar.SetStatusText("")
    # canvas_onMouseMove()
# class ScatterPlotFrame

class PlotFrame(wx.Frame):
    class PlotPanel(wx.Panel):
        def __init__(self, parent):
            wx.Panel.__init__(self, parent)#, size=(600,400))
            self.sizer = wx.BoxSizer(wx.VERTICAL)
            self.SetSizer(self.sizer)

            #matplotlib figure
            self.figure = matplotlib.figure.Figure(dpi=80) # this dpi is also needed.. right?
            self.figure.set_facecolor((0.7,0.7,1.))
            self.subplot = self.figure.add_subplot(111)
            self.colorbar = None
            self.points = []
            self.position_patch = None
            self.checked_patches = []

            #canvas
            self.canvas = matplotlib.backends.backend_wxagg.FigureCanvasWxAgg(self, wx.ID_ANY, self.figure)
            self.canvas.SetBackgroundColour(wx.Colour(100,255,255))
            self.canvas.mpl_connect('motion_notify_event', self.canvas_onMouseMove)
            self.canvas.mpl_connect('button_press_event', self.canvas_onMouseDown)
            self.canvas.mpl_connect('key_press_event', self.canvas_onKeyDown)
            self.sizer.Add(self.canvas, 1, flag=wx.LEFT|wx.TOP|wx.GROW)

            # for efficient annotation update
            self._blit_cache = None
            self._resize_id = self.figure.canvas.mpl_connect('resize_event',
                                                             self._handle_resize)

            self.reset()
        # __init__()

        def _handle_resize(self, *args):
            self.figure.canvas.mpl_disconnect(self._resize_id)
            self._blit_cache = None
            ###self._init_draw() # => self._drawn_artists = self._init_func()
            if self.position_patch is not None: self.position_patch.set_visible(False)
            for p in self.checked_patches: p.set_visible(False)
            self._resize_id = self.figure.canvas.mpl_connect('draw_event',
                                                             self._end_redraw)
        # _handle_resize()

        def _end_redraw(self, evt):
            self._post_draw()
            self.figure.canvas.mpl_disconnect(self._resize_id)
            self._resize_id = self.figure.canvas.mpl_connect('resize_event',
                                                             self._handle_resize)
        # _end_redraw()

        def _post_draw(self):
            if self._blit_cache is None:
                self._blit_cache = self.figure.canvas.copy_from_bbox(self.figure.axes[0].bbox)

            if self.position_patch is not None and self.position_patch.axes is not None:
                self.position_patch.set_visible(True)
                self.position_patch.axes.draw_artist(self.position_patch)

            for p in self.checked_patches:
                if p.hide: continue
                p.set_visible(True)
                p.axes.draw_artist(p)

            self.figure.canvas.blit(self.figure.axes[0].bbox)
        # _post_draw()

        def _pre_draw(self):
            if self.position_patch is None and len(self.checked_patches)==0:
                return

            self.figure.canvas.restore_region(self._blit_cache)
        # _pre_draw()

        def reset(self):
            plotFrame = self.GetParent().GetParent().GetParent().GetParent()

            for p in self.points:
                p.remove()

            self.points = []
            self.plotted_xy = None
            self.plotted_data = []
            self.kdtree = None # for fast lookup of nearest neighbour
            self.current_plotted_imgfs = []
            self.current_idx_mouse_on = None
            self._blit_cache = None
            self.subplot.set_title("")
            self.remove_annotate()
            plotFrame.peakPanel.clear()
            self.SetSize((self.Size[0],self.Size[1])) # Clear drawn plot
            #self.canvas.draw()
        # reset()

        def remove_annotate(self, refresh=True):
            if self.position_patch is not None:
                self.position_patch.remove()
                self.position_patch = None

            for p in self.checked_patches: p.remove()
            self.checked_patches = []

            if refresh:
                self.SetSize((self.Size[0],self.Size[1])) # Refresh drwan plot
        # remove_annotate()

        def canvas_onMouseMove(self, event):
            plotFrame = self.GetParent().GetParent().GetParent().GetParent()
            if None not in (event.xdata, event.ydata) and self.plotted_xy is not None:
                dist, idx = self.kdtree.query((event.xdata, event.ydata), k=1, p=1)
                x, y = self.plotted_xy[idx]
                imgf = os.path.basename(self.current_plotted_imgfs[idx])
                data = self.plotted_data[idx]
                scaninfo = plotFrame.find_data_by_filename(self.current_plotted_imgfs[idx]).scan_info
                vp, vs = scaninfo.vpoints, scaninfo.vstep*1000.
                hp, hs = scaninfo.hpoints, scaninfo.hstep*1000.
                dx, dy = abs(x-event.xdata), abs(y-event.ydata)
                if (vp==hp==1 and dx<.5 and dy<.5) or (vp==1 and dx < hs/2) or (hp==1 and dy < vs/2) or dx < hs/2 and dy < vs/2:
                    plotFrame.statusbar.SetStatusText("x= %.1f, y= %.1f, data= %.1f, file= %s" % (x, y, data, imgf))
                    self.current_idx_mouse_on = idx
                else:
                    plotFrame.statusbar.SetStatusText("")
                    self.current_idx_mouse_on = None
        # canvas_onMouseMove()

        def canvas_onMouseDown(self, event):
            # Sometimes 'button_press_event' does not seem to release mouse..
            # Does this fix a bug?
            if self.canvas.HasCapture():
                self.canvas.ReleaseMouse()

            idx = self.current_idx_mouse_on
            if idx is None:
                return

            plotFrame = self.GetParent().GetParent().GetParent().GetParent()
            plotFrame.peakPanel.select_imgf(self.current_plotted_imgfs[idx]) # which calls self.select_imgf()

            if event.button == 3: # right click
                # want to show popupmenu.
                pass
        # canvas_onMouseDown()

        def select_imgf(self, imgf):
            plotFrame = self.GetParent().GetParent().GetParent().GetParent()
            mainFrame = plotFrame.GetParent()

            imgfb = os.path.basename(imgf)
            print("Selected:", imgfb)

            # Update main window
            mainFrame.grid.load(imgf)

            # Show scatter plot
            plotFrame.splotFrame.reset()

            f, kind = plotFrame.get_selected_f_kind()
            sel = [x for x in plotFrame.data[f] if os.path.basename(x[0])==imgfb]
            if len(sel) == 0:
                return

            res_outer = sel[0][1].params.distl.res.outer
            plotFrame.splotFrame.plot(spots=sel[0][1].spots,
                                      mode=mainFrame.ctrlFrame.get_spot_draw_mode(),
                                      res_outer=res_outer,
                                      filename=imgf)

        # select_imgf()

        def canvas_onKeyDown(self, event):
            plotFrame = self.GetParent().GetParent().GetParent().GetParent()
            mainframe = plotFrame.GetParent()
            lc = plotFrame.peakPanel.listctrl
            isel = lc.GetFirstSelected()
            if event.key == "up":
                if isel > 0:
                    lc.Select(isel-1)
                    lc.EnsureVisible(isel-1)
            elif event.key == "down":
                if isel < lc.GetItemCount():
                    lc.Select(isel+1)
                    lc.EnsureVisible(isel+1)
            elif event.key == " ":
                lc.ToggleItem(isel)
            elif event.key in "hl": # <- H, -> L
                idx = list(mainframe.data.keys()).index(mainframe.grid.current_img_file)
                inc = 1 if event.key=="h" else -1
                if 0<= idx + inc < len(mainframe.data):
                    self.select_imgf(list(mainframe.data.keys())[idx+inc])
            elif event.key in "jk": # j:down, k:up
                imgf = mainframe.grid.current_img_file
                found = plotFrame.find_data_by_filename(imgf)
                sc, gc = found.scan_info, found.grid_coord
                inc = 1 if event.key=="j" else -1
                if sc.vpoints==1: return
                if sc.hpoints==1:
                    idx = list(mainframe.data.keys()).index(imgf)
                    if 0<= idx + inc < len(mainframe.data):
                        self.select_imgf(list(mainframe.data.keys())[idx+inc])
                else:
                    newgc = (gc[0], gc[1] - inc*sc.vstep)
                    fnd = plotFrame.find_data_by_gc(mainframe.ctrlFrame.current_target_fpref, newgc)
                    if fnd is not None: self.select_imgf(fnd[0])
        # canvas_onKeyDown()
    # class PlotPanel

    class PeakPanel(wx.Panel):
        def __init__(self, parent):
            wx.Panel.__init__(self, parent)#, size=(600,400)):
            vbox = wx.BoxSizer(wx.VERTICAL)
            self.SetSizer(vbox)

            self.listctrl = CheckListCtrl(self, style=wx.LC_REPORT|wx.LC_SINGLE_SEL)
            self.listctrl.InsertColumn(0, "score", wx.LIST_FORMAT_RIGHT, width=80) # with checkbox
            self.listctrl.InsertColumn(1, "x", wx.LIST_FORMAT_RIGHT, width=50)
            self.listctrl.InsertColumn(2, "y", wx.LIST_FORMAT_RIGHT, width=50)
            vbox.Add(self.listctrl, 1, wx.EXPAND)
            vbox.Add((5,5))
            self.stxtSelNum = wx.StaticText(self, label="  0 positions checked")
            vbox.Add(self.stxtSelNum)#, flag=wx.EXPAND|wx.ALL, border=1)
            self.btn_uncheck_all = wx.Button(self, wx.ID_ANY, "Uncheck all")
            vbox.Add(self.btn_uncheck_all, flag=wx.EXPAND|wx.ALL, border=1)

            self.btn_tell_kuma_checked = wx.Button(self, wx.ID_ANY, "Checked positions to KUMA")
            vbox.Add(self.btn_tell_kuma_checked, flag=wx.EXPAND|wx.ALL, border=1)

            # Auto-select
            self.sb = wx.StaticBox(self, label="Automatic select")
            sbsizer = wx.StaticBoxSizer(self.sb, wx.VERTICAL)
            hbox1 = wx.BoxSizer(wx.HORIZONTAL)
            hbox1.Add(wx.StaticText(self, label="Min score: "), flag=wx.RIGHT|wx.ALIGN_CENTER_VERTICAL)
            self.txtMinScore = wx.TextCtrl(self, wx.ID_ANY, "9", (95, 105))
            hbox1.Add(self.txtMinScore)
            sbsizer.Add(hbox1)
            hbox2 = wx.BoxSizer(wx.HORIZONTAL)
            hbox2.Add(wx.StaticText(self, label="Min distance (um): "), flag=wx.RIGHT|wx.ALIGN_CENTER_VERTICAL)
            self.txtMinDist = wx.TextCtrl(self, wx.ID_ANY, "15", (95, 105))
            hbox2.Add(self.txtMinDist)
            sbsizer.Add(hbox2)
            hbox3 = wx.BoxSizer(wx.HORIZONTAL)
            hbox3.Add(wx.StaticText(self, label="Max hits: "), flag=wx.RIGHT|wx.ALIGN_CENTER_VERTICAL)
            self.txtMaxHits = wx.TextCtrl(self, wx.ID_ANY, "Inf", (95, 105))
            hbox3.Add(self.txtMaxHits)
            sbsizer.Add(hbox3)
            self.btn_auto_select = wx.Button(self, wx.ID_ANY, "Auto-select")
            sbsizer.Add(self.btn_auto_select, flag=wx.EXPAND|wx.ALL, border=1)
            vbox.Add(sbsizer, flag=wx.EXPAND|wx.ALL, border=1)

            # wx.EVT_LIST_ITEM_FOCUSED is invoked when simple click, but also invoked after destroyed.. not knowing why..
            self.listctrl.Bind(wx.EVT_LIST_ITEM_SELECTED, self.listctrl_item_selected)
            self.listctrl.Bind(wx.EVT_LIST_KEY_DOWN, self.listctrl_item_key_down)
            self.btn_tell_kuma_checked.Bind(wx.EVT_BUTTON, self.btn_tell_kuma_checked_clicked)
            self.btn_uncheck_all.Bind(wx.EVT_BUTTON, self.btn_unckeck_all_clicked)
            self.btn_auto_select.Bind(wx.EVT_BUTTON, self.btn_auto_select_clicked)
            self.listctrl.OnCheckItem = self.listctrl_item_checked

            self.clear = self.listctrl.DeleteAllItems
            self.perm = []
        # __init__()

        def update_list(self):
            lc = self.listctrl
            lc.DeleteAllItems()
            p = self.GetParent().GetParent().plotPanel

            perm = sorted(list(range(len(p.plotted_data))), key=lambda x: -p.plotted_data[x])

            for i, k in enumerate(perm):
                lc.InsertStringItem(i, "%d"%p.plotted_data[k])
                lc.SetStringItem(i, 1, "%.1f"%p.plotted_xy[k][0])
                lc.SetStringItem(i, 2, "%.1f"%p.plotted_xy[k][1])

            self.perm = perm
        # update_list()

        def listctrl_item_selected(self, event):
            p = self.GetParent().GetParent().plotPanel
            
            idx = event.GetIndex()
            if len(self.perm) <= idx:
                print("error.")
                return

            p.select_imgf(p.current_plotted_imgfs[self.perm[idx]])
            p.SetFocus()
        # listctrl_item_selected()

        def listctrl_item_key_down(self, event):
            if event.GetCode() == wx.WXK_SPACE:
                self.listctrl.ToggleItem(self.listctrl.GetFirstSelected())
        # listctrl_item_key_down()
            
        def btn_tell_kuma_checked_clicked(self, event):
            p = self.GetParent().GetParent().plotPanel
            mainFrame = p.GetParent().GetParent().GetParent().GetParent().GetParent() # toooo dirty!!

            if mainFrame.kuma_addr is None:
                shikalog.warning("KUMA address (host and port) is not set.")
                wx.MessageDialog(None, "KUMA address is not set!",
                                 "Error", style=wx.OK).ShowModal()
                return

            gonios = []

            for i in range(self.listctrl.GetItemCount()):
                if self.listctrl.IsChecked(i):
                    score = self.listctrl.GetItem(i, 0).GetText()
                    imgf = p.current_plotted_imgfs[self.perm[i]]
                    gonio_xyz_phi = mainFrame.get_gonio_xyz_phi_from_imgf(imgf)
                    comment = "%s: %s" % (score, os.path.splitext(os.path.basename(imgf))[0])
                    gonios.append((gonio_xyz_phi, comment))

            for gonio, comment in gonios:
                print("tranferring", gonio, comment)
                mainFrame.tell_kuma(gonio, comment, with_dialog=False)
        # btn_tell_kuma_checked_clicked()

        def save_selected_info(self, filename):
            # NEED refactoring  - duplicated code with btn_tell_kuma_checked_clicked!!
            p = self.GetParent().GetParent().plotPanel
            mainFrame = p.GetParent().GetParent().GetParent().GetParent().GetParent() # toooo dirty!!

            ofs = open(filename, "w")
            ofs.write("filename gx gy gz phi score\n")

            for i in range(self.listctrl.GetItemCount()):
                if self.listctrl.IsChecked(i):
                    score = self.listctrl.GetItem(i, 0).GetText()
                    imgf = p.current_plotted_imgfs[self.perm[i]]
                    gonio_xyz_phi = mainFrame.get_gonio_xyz_phi_from_imgf(imgf)

                    data = (os.path.basename(imgf),) + tuple(gonio_xyz_phi) + (score,)
                    ofs.write("%s %s %s %s %s %s\n"%data)
        # save_selected_info()
        
        def btn_unckeck_all_clicked(self, event):
            if sum([self.listctrl.IsChecked(x) for x in range(self.listctrl.GetItemCount())]) == 0:
                return

            if wx.MessageDialog(None, "All items will be unchecked and this *CANNOT* be undone. Are you sure?",
                                "Confirm", style=wx.YES_NO|wx.NO_DEFAULT).ShowModal() == wx.ID_YES:
                for i in range(self.listctrl.GetItemCount()):
                    if self.listctrl.GetItem(i).GetImage() == 1: self.listctrl.SetItemImage(i, 0)
                self.listctrl_item_checked(None, None)
        # btn_unckeck_all_clicked()

        def btn_auto_select_clicked(self, event):
            lc = self.listctrl
            min_score = float(self.txtMinScore.GetValue())
            min_dist_sqr = float(self.txtMinDist.GetValue())**2
            max_hits = float(self.txtMaxHits.GetValue()) # to treat Inf

            # Register already-checked items
            registered = []
            for i in range(lc.GetItemCount()):
                if lc.IsChecked(i):
                    x, y = [float(lc.GetItem(i, x).GetText()) for x in range(1, 3)]
                    registered.append((x, y))

            # Traverse listctrl
            checked = []
            for i in range(lc.GetItemCount()):
                if len(registered) >= max_hits:
                    break

                score, x, y = [float(lc.GetItem(i, x).GetText()) for x in range(3)]

                if score < min_score:
                    break
                
                # dumn method.. would be slow if many positions registered
                min_ever = None
                if len(registered) > 0:
                    min_ever = min([(a[0]-x)**2+(a[1]-y)**2 for a in registered])

                if min_ever is None or min_ever >= min_dist_sqr:
                    checked.append(i)
                    registered.append((x, y))
                    
            shikalog.info("Auto-selection found %d positions (min_score= %.1f, min_dist= %.1f, max_hits= %.0f)" % (len(checked), min_score, math.sqrt(min_dist_sqr), max_hits))

            for i in checked:
                if lc.GetItem(i).GetImage() == 0:
                    lc.SetItemImage(i, 1)

            if len(checked) > 0:
                self.listctrl_item_checked(None, None)
        # btn_auto_select_clicked()
        
        def select_imgf(self, imgf):
            p = self.GetParent().GetParent().plotPanel
            for i in range(self.listctrl.GetItemCount()):
                if imgf == p.current_plotted_imgfs[self.perm[i]]:
                    self.listctrl.Select(i)
                    self.listctrl.EnsureVisible(i)
                    break
        # select_imgf()

        def listctrl_item_checked(self, index, flag):
            p = self.GetParent().GetParent().plotPanel
            plotFrame = p.GetParent().GetParent().GetParent().GetParent()
            imgfs = [p.current_plotted_imgfs[self.perm[i]] for i in [x for x in range(self.listctrl.GetItemCount()) if self.listctrl.IsChecked(x)]]
            self.stxtSelNum.SetLabel("%3d positions checked" % len(imgfs))
            plotFrame.annotate_checked(imgfs)
        # listctrl_item_checked()
    # class PeakPanel

    def __init__(self, parent=None, id=wx.ID_ANY):
        wx.Frame.__init__(self, parent=parent, id=id, title="Plot",
                          size=(800,600))
        self.Bind(wx.EVT_CLOSE, lambda e: self.Hide()) # Don't destroy this frame when closed
        self.Bind(wx.EVT_KEY_UP, self.OnKeyUp)

        self.splitter1 = wx.SplitterWindow(self, id=wx.ID_ANY)
        self.leftPanel = wx.Panel(self.splitter1)
        self.peakPanel = self.PeakPanel(self.splitter1)

        self.splitter1.SetSashGravity(1.0)
        self.splitter1.SplitVertically(self.leftPanel, self.peakPanel)
        self.splitter1.SetSashPosition(600)

        # Left panel
        self.splitter = wx.SplitterWindow(self.leftPanel, id=wx.ID_ANY)
        self.splitter.SetSashGravity(1.0)
        self.plotPanel = self.PlotPanel(self.splitter)
        self.panel = panel = wx.Panel(self.splitter)
        self.splitter.SplitHorizontally(self.plotPanel, self.panel)
        self.splitter.SetSashPosition(500)

        self.leftPanel.SetSizer(wx.BoxSizer(wx.HORIZONTAL))
        self.leftPanel.GetSizer().Add(self.splitter, 1, wx.EXPAND) # expand left panel

        vbox = wx.BoxSizer(wx.VERTICAL) # includes hbox and splotPanel

        hbox = wx.BoxSizer(wx.HORIZONTAL) # includes vbox11 and vbox12
        vbox.Add(hbox, 1, wx.EXPAND)

        vbox11 = wx.BoxSizer(wx.VERTICAL)
        hbox.Add(vbox11, 1,flag=wx.EXPAND|wx.LEFT, border=4)
        self.rb_kind = []
        self.rb_kind.append(wx.RadioButton(panel, wx.ID_ANY, "total_integrated_signal", style=wx.RB_GROUP))
        self.rb_kind.append(wx.RadioButton(panel, wx.ID_ANY, "median_integrated_signal"))
        self.rb_kind.append(wx.RadioButton(panel, wx.ID_ANY, "n_spots"))
        self.rb_kind[-1].SetValue(True) # Set n_spot as default.
        for rb in self.rb_kind:
            vbox11.Add(rb)
            rb.Bind(wx.EVT_RADIOBUTTON, self.rb_clicked)

        # right of rb_kind
        vbox12 = wx.BoxSizer(wx.VERTICAL)
        hbox.Add(vbox12, 1,flag=wx.EXPAND|wx.LEFT, border=4)
        self.chkShowAnnotate = wx.CheckBox(panel, wx.ID_ANY, "Show selected position")
        self.chkShowAnnotate.SetValue(True)
        self.chkShowAnnotate.Bind(wx.EVT_CHECKBOX, self.chkShowAnnotate_onCheck)
        vbox12.Add(self.chkShowAnnotate)
        self.chkGridPlot = wx.CheckBox(panel, wx.ID_ANY, "Grid plot")
        self.chkGridPlot.SetValue(True)
        self.chkGridPlot.Bind(wx.EVT_CHECKBOX, self.chkGridPlot_onCheck)
        vbox12.Add(self.chkGridPlot)


        self.splotFrame = ScatterPlotFrame(self) # scatter plot (I vs d^-2)

        panel.SetSizer(vbox)

        self.data = collections.OrderedDict()

        self.statusbar = self.CreateStatusBar()
    # __init__()

    def find_data_by_filename(self, filename):
        # TODO probably inefficient way.
        # Better to find data in main frame?
        for fpref in self.data:
            fltr = [s for s in self.data[fpref] if s[0]==filename]
            if len(fltr) > 0:
                return fltr[0][1]
        return None
    # find_data_by_filename()

    def find_data_by_gc(self, fpref, gc):
        # TODO inefficient way.
        if fpref not in self.data: return None
        tocmp = lambda x: [int(y*1e4+.5) for y in x] # in .1 micron precision
        gc = tocmp(gc)
        fltr = [s for s in self.data[fpref] if tocmp(s[1].grid_coord)==gc]
        if len(fltr) > 0:
            return fltr[0] # return filename, too. (not fltr[0][1])
        return None
    # find_data_by_gc()

    def OnKeyUp(self,event):
        if event.ControlDown() and event.GetKeyCode() == ord("R"):
            self.splotFrame.Show()
            self.splotFrame.Raise()

    def rb_clicked(self, event):
        """
        Find selected radio button and make a plot.
        """

        f, kind = self.get_selected_f_kind()
        if None not in (f, kind):
            self.SetTitle("Plot - %s" % f)
            self.plot(f, kind)
        self.splitter.SizeWindows() # Fit plots
        self.plotPanel.figure.canvas.draw()
        self.plotPanel._post_draw() # cache the plot

        if gui_params.auto_mode and None not in (f, kind) and f in self.data:
            work_dir = os.path.join(os.path.dirname(self.data[f][0][0]), "_spotfinder")
            scaninfo = self.data[f][0][1].scan_info
            assert scaninfo is not None
            vp, hp = scaninfo.vpoints, scaninfo.hpoints
            self.peakPanel.btn_unckeck_all_clicked(None)
            self.peakPanel.btn_auto_select_clicked(None)
            tmp = f if " (" not in f else f.split()[0] # just in case (must always include " (phi=..")
            outf = os.path.join(work_dir, "%sselected.dat"%tmp)
            shikalog.info("Auto-saving for KUMA: %s"%outf)
            self.peakPanel.save_selected_info(outf)

            mainframe = self.GetParent()
            mainframe.html_maker_thread.make_dat(work_dir)

            if len(self.data[f]) == vp*hp:
                open(outf, "a").write("#scan_complete\n")

            ### START HTML MAKER (here annotation is finished)
                if gui_params.mode == "zoo":
                    mainframe.html_maker_thread.make(work_dir, False)

        mainframe = self.GetParent()

        # Select best result (in track-the-latest-result mode)
        #if 1: # I think this should always be done.. <= Very annoying if this was done during scan!!
        if mainframe.ctrlFrame.chkTrackLatest.GetValue():
            data, imgfs = self.plotPanel.plotted_data, self.plotPanel.current_plotted_imgfs
            if len(data) != len(imgfs) or len(data) == 0:
                shikalog.error("Mismatch or zero length; len(data)= %d, len(imgfs)= %d" % (len(data), len(imgfs)))
                return
            max_i, max_v = max(enumerate(data), key=lambda x:(x[1], x[0])) # Better score and later image
            shikalog.info("Selecting the best score image: %s (score= %.2f)" % (imgfs[max_i], max_v))
            mainframe.grid.load(imgfs[max_i])
    # rb_clicked()

    def get_selected_f_kind(self):
        seldir, file_sel = self.GetParent().ctrlFrame.get_selected_dir_fpref()

        kind_sel = [rb for rb in self.rb_kind if rb.GetValue()]
        if file_sel != "" and len(kind_sel) > 0:
            f = file_sel
            kind = kind_sel[0].GetLabelText()
            return f, kind
        else:
            return None, None
    # get_selected_f_kind()

    def decide_fpref(self, f, scaninfo):
        fpref = re_pref_num_ext.search(os.path.basename(f)).group(1)
        if scaninfo is not None:
            if scaninfo.is_shutterless():
                fpref += " (phi=%.2f)" % (scaninfo.fixed_spindle)
            else:
                fpref += " (phi=%.2f)" % (scaninfo.osc_start)
        return fpref
    # decide_fpref()

    def set_data(self, result, append=False):
        def find_changing_gonio_axis(gonios):
            if len(gonios) < 2:
                return [0]

            ret = [] # list of True/False
            for i in range(3):
                i_diff_max = max([g[i]-gonios[0][i] for g in gonios[1:]])
                if i_diff_max >= 1e-4:
                    ret.append(i)

            return ret
        # find_changing_gonio_axis()

        #changing_axis = find_changing_gonio_axis([stat.gonio for f,stat in result])
        #print "Changing=", changing_axis
        self.data = collections.OrderedDict()
        mainframe = self.GetParent()

        sorted_result = list(result.items())
        try:
            sorted_result.sort(key=lambda x:self.diffscan_manager.get_scan_info(x[0]).date)
        except:
            shikalog.warning("Can't sort result by date.")

        f_fpref = []
        for f, stat in sorted_result:
            if os.path.dirname(f) != mainframe.ctrlFrame.current_target_dir:
                continue

            if stat is None: continue

            fpref = self.decide_fpref(f, stat.scan_info)

            #self.data.setdefault(fpref, []).append((os.path.basename(f), stat))
            self.data.setdefault(fpref, []).append((f, stat))
            f_fpref.append((os.path.dirname(f), fpref))

        #print "DEBUG:: f_fpref=", f_fpref

        ## append item on treectrl
        seen = set()
        for f, fpref in (x for x in f_fpref if not (x in seen or seen.add(x))):
            dic4t = mainframe.ctrlFrame.dic_for_tree
            keypar = tuple(os.path.relpath(f, mainframe.topdir).split(os.sep))
            #if keypar==('.',): keypar = () # this doesn't fix a problem..
            key = keypar + (fpref,)
            if key not in dic4t:
                #print "keys=", dic4t.keys(), "key=", key, "keypar=", keypar
                dic4t[key] = mainframe.ctrlFrame.treectrl.AppendItem(dic4t[keypar], fpref, image=3)
                mainframe.ctrlFrame.treectrl.SetPyData(dic4t[key], os.sep.join(key))

    # set_data()

    def annotate(self, imgf, lifetime=0):
        if not self.chkShowAnnotate.GetValue():
            return

        if imgf not in self.plotPanel.current_plotted_imgfs:
            self.plotPanel.SetSize((self.plotPanel.Size[0], self.plotPanel.Size[1]))
            shikalog.error("%s is not in current_plotted_imgfs" % imgf)
            return

        gc = self.find_data_by_filename(imgf).grid_coord
        sc = self.find_data_by_filename(imgf).scan_info
        vp, vs = sc.vpoints, sc.vstep
        hp, hs = sc.hpoints, sc.hstep
        if vp==1: vs = 5e-3
        if hp==1: hs = 5e-3
        if self.plotPanel.position_patch is None:
            self.plotPanel.position_patch = Rectangle(((gc[0]-hs/2.)*1000., (gc[1]-vs/2.)*1000.),
                                                      hs*1000., vs*1000., fill=None, edgecolor="cyan",
                                                      alpha=1, linewidth=2)
            p = self.plotPanel.position_patch

            self.plotPanel.figure.canvas.draw()
            self.plotPanel._post_draw()
            self.plotPanel.subplot.add_patch(p)
            p.axes.draw_artist(p)
            p.axes.figure.canvas.blit(p.axes.bbox)
        else:
            p = self.plotPanel.position_patch
            self.plotPanel._pre_draw()
            p.xy = ((gc[0]-hs/2.)*1000., (gc[1]-vs/2.)*1000.)
            self.plotPanel._post_draw()

    # annotate()

    def annotate_checked(self, imgfs):
        """
        show checked (in peak list) positions.
        """

        for p in self.plotPanel.checked_patches:
            #p.set_visible(False)
            p.hide = True # add property

        self.plotPanel._pre_draw()

        for i, imgf in enumerate(imgfs):
            if imgf not in self.plotPanel.current_plotted_imgfs:
                continue

            gc = self.find_data_by_filename(imgf).grid_coord
            sc = self.find_data_by_filename(imgf).scan_info
            vp, vs = sc.vpoints, sc.vstep
            hp, hs = sc.hpoints, sc.hstep
            if vp==1: vs = 1e-3
            if hp==1: hs = 1e-3
            
            if len(self.plotPanel.checked_patches) <= i:
                p = Ellipse((gc[0]*1000., gc[1]*1000.),
                            hs*1000., vs*1000., fill=None, edgecolor="green",
                            alpha=1, linewidth=2)
                p.hide = False
                self.plotPanel.checked_patches.append(p)
                self.plotPanel.subplot.add_patch(p)
            else:
                self.plotPanel.checked_patches[i].center = (gc[0]*1000., gc[1]*1000.)
                self.plotPanel.checked_patches[i].hide=False#set_visible(True)

        self.plotPanel._post_draw()
    # annotate_checked()

    def plot_grid(self, xs, ys, ds, scaninfo):
        #import scipy.interpolate

        xlim = min(xs), max(xs)
        ylim = min(ys), max(ys)

        if scaninfo is not None:
            vs, hs = scaninfo.vstep*1000., scaninfo.hstep*1000.
            if scaninfo.vpoints == 1: vs = 5
            if scaninfo.hpoints == 1: hs = 5
        else:
            vs, hs = 5, 5

        zi = numpy.zeros((int((ylim[1]-ylim[0])/vs+1.5),
                          int((xlim[1]-xlim[0])/hs+1.5)))
        for x, y, d in zip(xs, ys, ds):
            i,j = int((y-ylim[0])/vs+0.5), int((x-xlim[0])/hs+0.5)
            zi[i,j] = d

        p1 = self.plotPanel.subplot.imshow(zi, origin='lower',
                                           extent=[min(xs)-hs/2, max(xs)+hs/2,
                                                   min(ys)-vs/2, max(ys)+vs/2],
                                           interpolation='none', cmap="YlOrRd")#PuRd

        if max(ds) - min(ds) > 1e-5: # If all values equal (maybe), colorbar() will cause segmentation fault.
            cax = self.plotPanel.colorbar.ax if self.plotPanel.colorbar is not None else None
            self.plotPanel.colorbar = self.plotPanel.figure.colorbar(p1, cax=cax)

        return p1,
    # plot_grid()

    def plot_circles(self, xs, ys, ds, zero_xs, zero_ys):
        def normalize(v, m=100., sd=60.):
            vm = float(sum(v))/float(len(v))
            vsd = math.sqrt(sum([(x-vm)**2 for x in v])/float(len(v)))
            if vsd < 1e-12:
                return [m for x in range(len(v))]
            return [sd*(x-vm)/vsd+m for x in v]
        # normalize()

        def normalize_max(v, maximum=400.):
            max_v = max(v)
            f = maximum / max_v if max_v > 0 else 1.
            return [f*x + 1. for x in v] # add 1 to make zero-value pickable
        # normalize_max()

        p1 = self.plotPanel.subplot.scatter(xs, ys, s=normalize_max(ds), c=ds, alpha=0.5)
        if max(ds) - min(ds) > 1e-5: # If all values equal (maybe), colorbar() will cause segmentation fault.
            cax = self.plotPanel.colorbar.ax if self.plotPanel.colorbar is not None else None
            self.plotPanel.colorbar = self.plotPanel.figure.colorbar(p1, cax=cax)

        p2 = self.plotPanel.subplot.scatter(zero_xs, zero_ys, s=50, marker="x", c=[0]*len(zero_xs), alpha=0.5)
        return p1, p2
    # plot_circles()

    def plot(self, f, kind):
        if len(self.data) == 0:
            return

        ctrlframe = self.GetParent().ctrlFrame
        mode = ctrlframe.get_spot_draw_mode()

        if mode == "do not show spots":
            return

        # Clear plot
        self.plotPanel.reset()


        xs, ys, ds, imgfs = [], [], [], []
        zero_xs, zero_ys = [], [] # For values of zero
        for imgf, stat in self.data[f]:
            gc = stat.grid_coord
            if gc is None:
                shikalog.warning("gc is None! %s"%imgf)
                continue
            x, y = gc
            x *= 1000.
            y *= 1000.
            d = stat.stats[("n_spots","total_integrated_signal","median_integrated_signal").index(kind)]

            xs.append(x)
            ys.append(y)
            ds.append(d)
            imgfs.append(imgf)

            if d == 0:
                zero_xs.append(x)
                zero_ys.append(y)

        if len(xs) == 0:
            return

        scaninfo = self.data[f][0][1].scan_info

        if self.chkGridPlot.GetValue():
            self.plotPanel.points = self.plot_grid(xs, ys, ds, scaninfo)
        else:
            self.plotPanel.points = self.plot_circles(xs, ys, ds, zero_xs, zero_ys)

        self.plotPanel.subplot.set_xlabel("horizontal [um]")
        self.plotPanel.subplot.set_ylabel("vertical [um]")

        if scaninfo is not None:
            vp, hp = scaninfo.vpoints, scaninfo.hpoints
            vs, hs = scaninfo.vstep*1000., scaninfo.hstep*1000.

            self.plotPanel.subplot.set_title("%d out of %d (h=%d,v=%d) processed" % (len(xs), vp*hp, hp, vp))

            if 1 in (vp, hp) or len(self.data[f]) <= hp:
                self.plotPanel.subplot.set_aspect("auto")
            else:
                self.plotPanel.subplot.set_aspect("equal")

            # Set limits
            if vp == hp == 1:
                self.plotPanel.subplot.set_xlim(-10, 10)
                self.plotPanel.subplot.set_ylim(-10, 10)
            elif vp == 1:
                self.plotPanel.subplot.set_xlim(-hs*hp/2 - hs, hs*hp/2 + hs)
                self.plotPanel.subplot.set_ylim(-10, 10)
            elif hp == 1:
                self.plotPanel.subplot.set_xlim(-10, 10)
                self.plotPanel.subplot.set_ylim(-vs*vp/2 - vs, vs*vp/2 + vs)
            else:
                self.plotPanel.subplot.set_xlim(-hs*hp/2 - hs, hs*hp/2 + hs)
                self.plotPanel.subplot.set_ylim(-vs*vp/2 - vs, vs*vp/2 + vs)
        else:
            # Should never reach here.. but should we set limit here?
            pass

        self.plotPanel.current_plotted_imgfs = imgfs
        self.plotPanel.plotted_xy = numpy.column_stack((xs, ys))
        self.plotPanel.kdtree = scipy.spatial.cKDTree(self.plotPanel.plotted_xy)
        self.plotPanel.plotted_data = ds

        self.peakPanel.update_list()
    # plot()

    def chkShowAnnotate_onCheck(self, event):
        if self.chkShowAnnotate.GetValue():
            imgf = self.GetParent().grid.current_img_file
            self.annotate(imgf)
        else:
            self.plotPanel.remove_annotate()
    # chkShowAnnotate_onCheck()

    def chkGridPlot_onCheck(self, event):
        self.rb_clicked(None)
    # chkGridPlot_onCheck()
# class PlotFrame

class ImageSpotPanel(wx.Panel):
    def __init__(self, parent, size):
        wx.Panel.__init__(self, parent, size=size)
        self.parent = parent
        self.img = None
        self._imgin = None
        self._sxsyshsw = (0, 0, None, None)
        self.stats = None
        self._pos = None
        self._mag = None
        self.current_dmin, self.current_width = None, None
        self.Bind(wx.EVT_PAINT, self.Draw)
    # __init__()

    def set_image(self, imgin, posmag, sx=0,sy=0,sh=None,sw=None):
        #self._bitmap = wx.Bitmap(imgin) # This is very slow if many many images loaded!
        self._imgin = imgin
        self._sxsyshsw = (sx, sy, sh, sw)
        self._pos = posmag[0:2]
        self._mag = posmag[2]
        self.Refresh()
    # set_image()

    def clear(self):
        self._imgin = None
        dc = wx.PaintDC(self)
        dc.Clear()
    # clear()

    def set_stats(self, stats):
        self.stats = stats
        self.Refresh()
    # set_stats()

    def Draw(self, ev):
        dc = wx.PaintDC(ev.GetEventObject())
        rect = self.GetClientRect()
        sx, sy, sh, sw = self._sxsyshsw
        if self._imgin is None: return
        if not os.path.isfile(self._imgin): return
        
        if (sh, sw).count(None) == 2:
            _image = wx.MemoryDC(wx.Bitmap(self._imgin))
        else:
            # Reading directly into wx.Bitmap is extremely slow!
            wx.Log_EnableLogging(False)
            try:
                im = wx.Image(self._imgin)
                if not im.IsOk(): raise
                im = im.GetSubImage(wx.Rect(sx,sy,sh,sw))
                if not im.IsOk(): raise
                _image = wx.MemoryDC(im.ConvertToBitmap())
            except:
                shikalog.warning("Thumbnail load failed: %s" % self._imgin)
                return
            finally:
                wx.Log_EnableLogging(True)

        width, height = _image.GetSize()
            
        if width > rect.width-2:
            width = rect.width-2
        if height > rect.height-2:
            height = rect.height-2

        draw_rect = wx.Rect(rect.x, rect.y, width, height)
        
        dc.Blit(draw_rect.x, draw_rect.y, draw_rect.width, draw_rect.height, _image, 0, 0, wx.COPY, True)
        
        self.draw_spots(dc, draw_rect)
        self.draw_beamcenter(dc, draw_rect)
    # Draw()

    def draw_spots(self, dc, draw_rect):
        """
        draw_rect is the region of the diffraction image in the dc
        """

        if self.stats is None:
            return

        spots = self.stats.spots
        ctrlframe = self.parent.GetParent().GetParent().ctrlFrame
        mode = ctrlframe.get_spot_draw_mode()
        if mode == "do not show spots":
            return

        dc.SetBrush(wx.Brush(wx.BLUE, wx.TRANSPARENT))
        dc.SetPen(wx.Pen("red"))
        w, h = 7, 7
        for y, x, snr, d in spots:
            x, y = draw_rect.x + (x - self._pos[0])*self._mag, draw_rect.y + (y - self._pos[1])*self._mag
            rect = (x-w, y-h, w*2+1, h*2+1)
            #if draw_rect.ContainsRect(rect):
            if draw_rect.Contains((x, y)):
                dc.DrawRectangleRect(rect)
    # draw_spots()

    def draw_beamcenter(self, dc, draw_rect):
        """
        Just add + mark on the center of image.
        NOTE that image is assumed to be centered on beam position!
        """
        l = 10
        w, h = draw_rect.width, draw_rect.height
        xc, yc = draw_rect.x + w/2, draw_rect.y + h/2
        dc.SetPen(wx.Pen("blue"))
        dc.DrawLine(xc - l, yc, xc + l, yc)
        dc.DrawLine(xc, yc - l, xc, yc + l)
    # draw_beamcenter()
# class ImageSpotPanel

class ImageResultPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        self.parent = parent
        self.r2d = None # to be a function

        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.vbox)

        self.view1 = wx.html.HtmlWindow(self, style=wx.NO_BORDER, size=(600,90))
        self.view1.SetStandardFonts()

        self.panel1 = wx.Panel(self, size=(600, 10))
        panel1_hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.panel1.SetSizer(panel1_hbox)
        self.dtxt = wx.StaticText(self.panel1, size=(200, 10))
        self.lbtn = wx.Button(self.panel1, wx.ID_ANY, "<")
        self.rbtn = wx.Button(self.panel1, wx.ID_ANY, ">")
        panel1_hbox.Add(self.dtxt, 1, flag=wx.ALIGN_CENTER_VERTICAL)
        panel1_hbox.Add(self.lbtn, 0, flag=wx.EXPAND)
        panel1_hbox.Add(self.rbtn, 0, flag=wx.EXPAND)

        self.panel2 = ImageSpotPanel(self, size=(600,600))

        self.view3 = wx.html.HtmlWindow(self, style=wx.NO_BORDER, size=(600,200))
        self.view3.SetStandardFonts()

        self.vbox.Add(self.view1, 1, flag=wx.EXPAND)
        self.vbox.Add(self.panel1, 0, flag=wx.EXPAND)
        self.vbox.Add(self.panel2, 0, flag=wx.EXPAND)
        self.vbox.Add(self.view3, 1, flag=wx.EXPAND)
        self.clear()

        self.lbtn.Bind(wx.EVT_BUTTON, lambda e: self.Scroll(+1))
        self.rbtn.Bind(wx.EVT_BUTTON, lambda e: self.Scroll(-1))
        self.Bind(wx.EVT_KEY_UP, self.OnKeyUp)
        self.panel2.Bind(wx.EVT_MOTION, self.onMouseMoveInImage)
        self.panel2.Bind(wx.EVT_LEAVE_WINDOW, lambda e: self.dtxt.SetLabel(""))
        self.panel2.Bind(wx.EVT_LEFT_DOWN, self.onMouseClick)
    # __init__()

    def is_the_same_fpref(self, filename, current_fpref):
        """
        Test the image has the same fpref as provided
        """
        mainframe = self.GetParent().GetParent()
        fpref_ = mainframe.plotFrame.decide_fpref(filename,
                                                  mainframe.data[filename].stat.scan_info)
        return current_fpref == fpref_
    # is_the_same_fpref()

    def Scroll(self, inc):
        mainframe = self.GetParent().GetParent()
        if self.current_img_file is None or len(mainframe.data) < 2:
            return

        data = list(mainframe.data.keys())

        if self.current_img_file not in data:
            return

        idx = data.index(self.current_img_file)

        if 0<= idx + inc < len(data):
            self.load(data[idx+inc])
    # Scroll()

    def OnKeyUp(self, event):
        keycode = event.GetKeyCode()
        if keycode in (wx.WXK_UP, wx.WXK_DOWN):
            inc = -1 if keycode==wx.WXK_UP else 1
            self.Scroll(inc)
    # OnKeyUp()

    def clear(self):
        self.view3.SetPage("<b>Scan info</b><br><br>No scan selected.")
        self.view1.SetPage("<b>Image info</b><br><br>No image selected.")
        self.panel2.clear()
        self.current_img_file = None
    # clear()

    def update(self):
        mainframe = self.GetParent().GetParent()

        img_file = self.current_img_file
        if img_file is None:
            return

        if img_file in mainframe.data:
            item = mainframe.data[img_file]
            spot_mode = mainframe.ctrlFrame.get_spot_draw_mode()
            if spot_mode == "do not show spots":
                return

            self.panel2.set_stats(item.stat)
            self.show_image_info(os.path.basename(item.img_file), item.stat.detector, item.stat.stats, spot_mode)
            scaninfo = item.stat.scan_info

            tmp = item.stat.params
            if tmp:
                self.prepare_resolution_calc(tmp.distl.res.outer, scaninfo)
    # update()

    def load(self, img_file):
        mainframe = self.GetParent().GetParent()

        self.current_img_file = img_file
        self.update()

        if img_file in mainframe.data:
            item = mainframe.data[img_file]
            possible_paths = [os.path.join(os.path.dirname(item.img_file), "_spotfinder",
                                                          os.path.basename(item.img_file)+ext) for ext in (".jpg",".png")]
            tiled_jpg = None
            prefix, num = None, None
            r = re.search("^(.*)_([0-9]+)\.[^0-9]+$", os.path.basename(item.img_file))
            if r:
                prefix, num = r.group(1), int(r.group(2))
                possible_paths.append(os.path.join(os.path.dirname(item.img_file), "_spotfinder",
                                                   "thumb_%s_%.3d" % (prefix, num//1000), os.path.basename(item.img_file)+".jpg"))
                idx = (num-1)//100
                tiled_jpg = os.path.join(os.path.dirname(item.img_file), "_spotfinder",
                                         "thumb_%s" % prefix,
                                         "%s_%.6d-%.6d.jpg" % (prefix, idx*100+1, (idx+1)*100))
                
            img_pics = [f for f in possible_paths if os.path.exists(f)]
            if len(img_pics) > 0:
                self.show_image(img_pics[0], item.stat.thumb_posmag)
            elif os.path.isfile(tiled_jpg):
                thumbw = 600 # MAGIC NUMBER!
                idx2 = (num-1)%100
                x, y = idx2%10, idx2//10
                self.show_image(tiled_jpg, item.stat.thumb_posmag, x*thumbw, y*thumbw, thumbw, thumbw)
            else:
                shikalog.warning("Image for display is not found: %s" % item.img_file)

            scaninfo = item.stat.scan_info
            self.show_scan_info(scaninfo)

            # Decide next and prev buttons available
            data = list(mainframe.data.keys())
            idx = data.index(self.current_img_file)
            is_valid_idx = lambda i: 0<= i < len(data)

            current_fpref = mainframe.plotFrame.decide_fpref(self.current_img_file,
                                                             mainframe.data[self.current_img_file].stat.scan_info)
            if is_valid_idx(idx+1) and self.is_the_same_fpref(data[idx+1], current_fpref):
                self.lbtn.Enable()
            else:
                self.lbtn.Disable()

            if is_valid_idx(idx-1) and self.is_the_same_fpref(data[idx-1], current_fpref):
                self.rbtn.Enable()
            else:
                self.rbtn.Disable()

        else:
            shikalog.error("Not found: " + img_file)

        # Update annotation after showing image.
        try: wx.SafeYield()
        except: pass
        mainframe.plotFrame.annotate(img_file)
        mainframe.plotFrame.peakPanel.select_imgf(img_file)
    # load()

    def onMouseClick(self, ev):
        mainframe = self.GetParent().GetParent()
        if self.current_img_file is not None:
            mainframe.plotFrame.annotate(self.current_img_file)
    # onMouseClick()

    def refresh_image(self):
        self.panel2.Refresh()
    # refresh_image()

    def show_image(self, imgpic, thumb_posmag, sx=0,sy=0,sh=None,sw=None):
        self.panel2.set_image(imgpic, thumb_posmag, sx,sy,sh,sw)
    # show_image()

    def show_scan_info(self, info):
        html = "<b>Scan info</b><br>"  #"<h3>Scan info</h3>"
        html += "<table>\n"
        html += '<tr align="left"><th>scan</th><td>%s</td>\n' % info.filename_template
        html += '<th>date</th><td>%s</td></tr>\n' % (info.date.strftime("%Y/%m/%d %H:%M:%S") if info.date!=0 else "??")

        if info.is_shutterless():
            html += '    <tr align="left"><th>fixed spindle</th><td>%.2f&deg;</td>\n' % info.fixed_spindle
            html += '                     <th>frame rate</th><td>%.2f [Hz]</td></tr>\n' % info.frame_rate
        else:
            html += '    <tr align="left"><th>osc. start</th><td>%.2f&deg;</td>\n' % info.osc_start
            html += '                     <th>osc. step</th><td>%.2f&deg;</td></tr>\n' % info.osc_step
            html += '    <tr align="left"><th>exp. time</th><td>%.2f [sec]</td></tr>\n' % info.exp_time

        html += '    <tr align="left"><th>beam size</th><td>h= %.1f, v= %.1f [um]</td>\n' % (info.beam_hsize, info.beam_vsize)
        html += '                     <th>attenuator</th><td>%s %.1f [um]</td></tr>\n' % info.attenuator
        html += '    <tr align="left"><th>distance</th><td>%.2f [mm]</td>\n' % info.distance
        html += '                     <th>wavelength</th><td>%.4f [A]</td></tr>\n' % info.wavelength
        html += '    <tr align="left"><th>scan points</th><td>v=%d, h=%d</td>\n' % (info.vpoints, info.hpoints)
        html += '                     <th>scan steps</th><td>v=%.2f, h=%.2f [um]</td></tr>\n' % (info.vstep*1000., info.hstep*1000.)
        html += '    <tr align="left"><th>scan direction</th><td>%s</td>\n' % (getattr(info, "scan_direction","None"))
        html += '                     <th>scan path</th><td>%s</td></tr>\n' % (getattr(info, "scan_path", "None"))

        html += '  </table>\n'

        self.view3.SetPage(html)
    # show_scan_info()

    def show_image_info(self, filename, detector, stats, spot_mode):
        n_spots, total_sig, med_sig = stats

        color = "black"

        self.view1.SetPage("""\
<b>Image info</b><br>
<table>
<tr align="left"><th>File name</th><td>%(filename)s</td></tr>
<tr align="left"><th>Detector type</th><td>%(detector)s</td></tr>
<tr align="left"><th>Total integrated signal</th><td><font color="%(col)s">%(total_sig).1f</font></td></tr>
<tr align="left"><th>Median integrated signal</th><td><font color="%(col)s">%(med_sig).1f</font></td></tr>
<tr align="left"><th>N_spots</th><td><font color="%(col)s">%(n_spots)d</font></td></tr>
</table>
""" % dict(filename=filename, detector=" ", total_sig=total_sig, med_sig=med_sig, n_spots=n_spots, col=color)) # XXX detector=", ".join(map(str,detector))
    # show_image_info()

    def prepare_resolution_calc(self, dmin, scaninfo):
        """
        TODO:
        Currently, wavelength and distance are obtained from scaninfo.
        Probably they should be obtained from image header? or should be in params!
        """

        self.r2d = None
        if dmin is None:
            return 

        wavelen = scaninfo.wavelength
        distance = scaninfo.distance

        # conversion factor (pixel -> mm)
        f = 2. * distance / 600. * math.tan(2. * math.asin(wavelen/2./dmin))
        self.r2d = lambda r: wavelen / 2. / math.sin(.5 * math.atan2(f*r,distance))
    # prepare_resolution_calc()

    def onMouseMoveInImage(self, ev):
        width, height = 600, 600
        self.dtxt.SetLabel("")

        if self.r2d is None:
            return

        pt = ev.GetPosition()
        rect = self.panel2.GetClientRect()
        draw_rect = wx.Rect(rect.x, rect.y, width, height)

        if not draw_rect.Contains(pt):
            # Outside
            return

        # beam center
        xc, yc = draw_rect.x + width/2, draw_rect.y + height/2

        # resolution
        r = math.sqrt((pt.x - xc)**2 + (pt.y - yc)**2)
        d = self.r2d(r) if r > 1e-6 else float("inf")
        self.dtxt.SetLabel("d= %6.2f A" % d)
    # onMouseMoveInImage()
# class ImageResultPanel

class MainFrame(wx.Frame):
    class Item(object):
        def __init__(self, img_file):
            self.img_file = img_file # file name (absolute path)
            self.stat = None # distl_stat
    # class Item

    def __init__(self, parent=None, id=wx.ID_ANY, topdir=None, params=None):
        wx.Frame.__init__(self, parent=parent, id=id, title="SHIKA system"+(" (Zoo mode)" if gui_params.mode == "zoo" else ""),
                          size=(1110,950))
        self.adxv_proc = None # subprocess object
        self.adxv_port = 8100 # adxv's default port. overridden later.
        self.adxv_bin = params.adxv
        self.kuma_addr = params.kuma_addr
        self.imgview_host = params.imgview_host
        self.auto_mode = params.auto_mode

        self.Bind(wx.EVT_CLOSE, self.onClose)

        self.splitter = wx.SplitterWindow(self, id=wx.ID_ANY)
        self.ctrlFrame = ControlPanel(self, parent=self.splitter, params=params)

        self.grid = ImageResultPanel(self.splitter)
        self.splitter.SplitVertically(self.ctrlFrame, self.grid)
        self.splitter.SetSashPosition(500)

        self.topdir = topdir

        self.data = collections.OrderedDict() # Data shown in grid

        self.plotFrame = PlotFrame(self)

        self.readonly_mode = params.readonly
        if self.readonly_mode:
            shikalog.info("SHIKA is in read-only mode and will not write any files.")
            self.ctrlFrame.btnUpdate.Disable()

        self.html_maker_thread = ReportHTMLMakerThread(self, dont_work=self.readonly_mode, make_html=params.make_html)
        self.html_maker_thread.start()

        if self.topdir is not None:
            self.ctrlFrame.txtTopDir.SetValue(self.topdir)
            root = self.ctrlFrame.treectrl.AddRoot(self.topdir, image=0)
            self.ctrlFrame.treectrl.SetPyData(root, ".")
            self.ctrlFrame.dic_for_tree[()] = root
            wx.PostEvent(self.ctrlFrame, EventTargetDirChanged(target=self.topdir, fpref=None))

        self.ctrlFrame.cmbTargetDir.SetEditable(False) # This is a current limitation.

        self.grid.panel2.Bind(wx.EVT_LEFT_DCLICK, self.grid_OnDbClick)
        self.grid.panel2.Bind(wx.EVT_RIGHT_DOWN, self.grid_OnRightClick)

        self.Show()

    # __init__()

    def open_img_with_adxv(self, imgfile):
        """
        Start adxv and show image.
        If already started, just update image shown.
        There are maybe two ways to do this.
        1; use -autoload option and update a temporary file. (need two seconds to refresh)
        2; use -socket 8100 option and communicate ('load_image hoge.img').
        Which is better?
        """

        # Hdf5 workaround
        tmp = glob.glob(re.sub("_[0-9]*\.img$", "_master*.h5", imgfile)) # if binned, master_bin*.h5 exists instead.
        print(tmp)
        if tmp and not os.path.isfile(imgfile):
            h5master = tmp[0]
            from yamtbx.dataproc import eiger
            frameno = int(re.search(".*_([0-9]*)\.img$", imgfile).group(1))
            imgfile = os.path.join(tempfile.gettempdir(), "adxvtmp-%s-%s.cbf"%(getpass.getuser(), os.getpid()))
            eiger.extract_to_minicbf(h5master, frameno, imgfile)

        if self.adxv_bin is not None:
            adxv_comm = self.adxv_bin + " -socket %d"
        else:
            adxv_comm = "adxv -socket %d"

        if self.adxv_proc is None or self.adxv_proc.poll() is not None: # None means still running.
            # find available port number
            sock_test = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            sock_test.bind(("localhost", 0))
            self.adxv_port = sock_test.getsockname()[1]
            sock_test.close()
            # start adxv
            self.adxv_proc = subprocess.Popen(adxv_comm%self.adxv_port, shell=True,
                                              cwd=os.path.dirname(imgfile))

        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

        for i in range(4): # try for 2 seconds.
            try:
                sock.connect(("localhost", self.adxv_port))
                break
            except socket.error:
                time.sleep(.5)
                continue

        sent = sock.send("load_image %s\n"%imgfile)
        shikalog.debug("adxv loading %s"%imgfile)
        if sent == 0:
            shikalog.error("adxv load failed! Close adxv and double-click again.")

        sent = sock.send("raise_window Control\n") # raise_window is available from adxv 1.9.9
        sent = sock.send("raise_window Image\n")

        sock.close()
    # open_img_with_adxv

    def open_in_imgview(self, imgfile):
        if self.imgview_host is None:
            shikalog.error("imgview host is not configured!")
            return
        print("Trying",self.imgview_host, 5555)

        import telnetlib
        telnet = telnetlib.Telnet(self.imgview_host, 5555)
        telnet.write("put/video_file/%s\n"%imgfile)
        #print "READ=", telnet.read_all()
        recv = telnet.read_until("/ok", timeout=3)
        if recv == "":
            print("ERROR: imgview not responding!")
            telnet.close()
            return
        #print "READ=", telnet.read_very_eager()
        telnet.write("put/video/disconnect\n")
        print("DONE.")
        telnet.close()
        return

        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.settimeout(3.)
        sock.connect((self.imgview_host, 5555))

        try:
            sock.sendall("put/video_file/%s\n"%imgfile)
        except:
            print("ERROR: imgview load failed!")
            sock.close()
            return

        time.sleep(1)
        recv = sock.recv(4096)
        print("recv=", recv)

        sock.send("put/video/disconnect\n")
        sock.close()

    # open_in_imgview()

    def onClose(self, event):
        self.Destroy()
    # onClose()

    def grid_OnDbClick(self, event):
        """
        Start adxv
        """

        img_file = self.grid.current_img_file
        if img_file is None:
            print("No image")
            return

        self.open_img_with_adxv(img_file)
    # grid_OnDbClick()

    def get_gonio_xyz_phi_from_imgf(self, img_file):
        if img_file not in current_stats:
            shikalog.warning("gonio xyz and phi is unavailable for %s." % img_file)
            return None

        stat = current_stats[img_file]
        gonio_xyz_phi = stat.gonio
        if stat.scan_info.is_shutterless():
            gonio_xyz_phi += (stat.scan_info.fixed_spindle,)
        else:
            gonio_xyz_phi += (stat.scan_info.osc_start,)

        shikalog.info("file, gonio= %s, %s" % (img_file, gonio_xyz_phi))

        if None in gonio_xyz_phi:
            shikalog.warning("None in gonio xyz or phi")

        return gonio_xyz_phi
    # get_gonio_xyz_phi_from_imgf()

    def grid_OnRightClick(self, event):
        img_file = self.grid.current_img_file
        if img_file is None:
            shikalog.error("No image")
            return

        gonio_xyz_phi = self.get_gonio_xyz_phi_from_imgf(img_file)

        if self.kuma_addr is None:
            shikalog.warning("KUMA address (host and port) is not set.")

        menu = wx.Menu()
        menu.Append(0, os.path.basename(img_file))
        menu.Enable(0, False)
        menu.AppendSeparator()
        menu.Append(1, "Let KUMA know")
        menu.Enable(1, None not in (gonio_xyz_phi, self.kuma_addr))
        menu.Append(2, "Let KUMA know (quick)")
        menu.Enable(2, None not in (gonio_xyz_phi, self.kuma_addr))
        menu.Append(3, "Open with adxv")
        menu.Append(4, "Open in imgview")
        menu.Enable(4, self.imgview_host is not None)

        self.Bind(wx.EVT_MENU, lambda e: self.tell_kuma(gonio_xyz_phi, os.path.splitext(os.path.basename(img_file))[0]), id=1)
        self.Bind(wx.EVT_MENU, lambda e: self.tell_kuma(gonio_xyz_phi, os.path.splitext(os.path.basename(img_file))[0], False), id=2)
        self.Bind(wx.EVT_MENU, lambda e: self.open_img_with_adxv(img_file), id=3)
        self.Bind(wx.EVT_MENU, lambda e: self.open_in_imgview(img_file), id=4)
        self.PopupMenu(menu)
        menu.Destroy()

    # grid_OnRightClick()

    def tell_kuma(self, gonio_xyz_phi, comment, with_dialog=True):
        class Dlg(wx.Dialog):
            def __init__(self, parent, gonio_xyz_phi, comment, func):
                wx.Dialog.__init__(self, parent, wx.ID_ANY, "KUMA communicator", size=(250, 100))

                self.gonio_xyz_phi = gonio_xyz_phi
                self.func = func

                vbox = wx.BoxSizer(wx.VERTICAL)
                self.txtComment = wx.TextCtrl(self, wx.ID_ANY, comment, (95, 105))
                hbox = wx.BoxSizer(wx.HORIZONTAL)
                btnOK = wx.Button(self, wx.ID_ANY, 'OK', size=(70, 30))
                btnCancel = wx.Button(self, wx.ID_ANY, 'Cancel', size=(70, 30))
                hbox.Add(btnOK, 1)
                hbox.Add(btnCancel, 1, wx.LEFT, 5)

                vbox.Add(wx.StaticText(self, wx.ID_ANY, "Comment:"))
                vbox.Add(self.txtComment, 1, wx.GROW|wx.LEFT)
                vbox.Add(hbox, 1, wx.ALIGN_CENTER | wx.TOP | wx.BOTTOM, 10)
                self.SetSizer(vbox)

                self.txtComment.SetFocus()

                btnOK.Bind(wx.EVT_BUTTON, self.btnOK_click)
                btnCancel.Bind(wx.EVT_BUTTON, lambda e: self.Destroy())
            # __init__()
            def btnOK_click(self, event):
                try:
                    self.func(self.gonio_xyz_phi, self.txtComment.GetValue(), False)
                finally:
                    self.Destroy()
            # btnOK_click()
        # class Dlg

        if None in gonio_xyz_phi or gonio_xyz_phi is None:
            shikalog.error("Invalid gonio coordinate or phi")
            wx.MessageDialog(None, "Invalid gonio coordinate or phi",
                             "Error", style=wx.OK).ShowModal()
            return


        if with_dialog:
            dlg = Dlg(self, gonio_xyz_phi, comment, self.tell_kuma) # XXX Make sure this function is not called recursively!
            dlg.ShowModal()
            return

        #s = xmlrpclib.ServerProxy('http://192.168.163.2:1920')
        try:
            s = xmlrpc.client.ServerProxy('http://%s'%self.kuma_addr)
            s.append_coords(gonio_xyz_phi, comment)
        except socket.error as e:
            shikalog.error("Cannot communicate with KUMA: %s" % e)
            wx.MessageDialog(None, """\
Cannot communicate with KUMA.

Is KUMA up?
Is this address (%s) correct?
Is network working?

Try restarting KUMA and SHIKA!
"""% self.kuma_addr,
                             "Error", style=wx.OK).ShowModal()

        print(gonio_xyz_phi, comment)
    # tell_kuma()

    def update_result(self, append=False):
        result = current_stats

        for f, stat in list(result.items()):
            if f in self.data:
                self.data[f].stat = stat

        self.plotFrame.set_data(result, append=append)

        # FIXME if track latest, rb_clicked() will be called twice.
        self.ctrlFrame.rb_clicked(None, call_from_runbutton=True, append=append) # this calls plotFrame.rb_clicked()

        if self.ctrlFrame.chkTrackLatest.GetValue():
            wx.CallAfter(self.track_latest_result) # to call it after onTargetDirChanged()

        if gui_params.mode != "zoo":
            with self.html_maker_thread.lock:
                self.html_maker_thread.queue.append((os.path.join(self.ctrlFrame.current_target_dir, "_spotfinder"),
                                                     not append))
            if not self.html_maker_thread.is_running():
                shikalog.debug("html_maker_thread was accidentally stopped. restarting.")
                self.html_maker_thread.start()
    # update_result()

    def load_results(self):
        current_stats.clear()
        

        if self.ctrlFrame.current_target_dir is None:
            return

        dbfile = os.path.join(self.ctrlFrame.current_target_dir, "_spotfinder", "shika.db")
        if not os.path.isfile(dbfile): return

        scanlog = os.path.join(self.ctrlFrame.current_target_dir, "diffscan.log")
        if not os.path.isfile(scanlog):
            shikalog.error("diffscan.log not found in %s" % self.ctrlFrame.current_target_dir)
            return

        slog = bl_logfiles.BssDiffscanLog(scanlog)
        slog.remove_overwritten_scans()
        
        d = wx.lib.agw.pybusyinfo.PyBusyInfo("Loading saved results..", title="Busy SHIKA")

        try: wx.SafeYield()
        except: pass

        try:
            shikalog.info("Loading data: %s" % dbfile)
            startt = time.time()
            result = []
            con = sqlite3.connect(dbfile, timeout=10, isolation_level=None)
            shikalog.debug("Opening db with query_only = ON")
            con.execute('pragma query_only = ON;')
            cur = con.cursor()

            for itrial in range(60):
                try:
                    c = cur.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='status';")
                    if c.fetchone() is None:
                        shikalog.error("No 'status' in %s" % dbfile)
                        return
                    break
                except sqlite3.DatabaseError:
                    shikalog.warning("DB failed. retrying (%d)" % itrial)
                    time.sleep(1)
                    continue

            for itrial in range(60):
                try:
                    c = con.execute("select filename,spots from spots")
                    results = dict([(str(x[0]), pickle.loads(str(x[1]))) for x in c.fetchall()])
                    break
                except sqlite3.DatabaseError:
                    shikalog.warning("DB failed. retrying (%d)" % itrial)
                    time.sleep(1)
                    continue

            exranges = self.ctrlFrame.get_exclude_resolution_ranges()
            if exranges:
                shikalog.info("Applying resolution-range exclusion: %s" % exranges)
                for r in list(results.values()):
                    if not r["spots"]: continue
                    ress = numpy.array([x[3] for x in r["spots"]])
                    test = numpy.zeros(len(r["spots"])).astype(numpy.bool)
                    for rr in exranges: test |= ((min(rr) <= ress) & (ress <= max(rr)))
                    for i in reversed(numpy.where(test)[0]): del r["spots"][i]

            print("DEBUG:: scans=", slog.scans) 
            for scan in slog.scans:
                for imgf, (gonio, gc) in scan.filename_coords:
                    #print imgf, (gonio, gc) 
                    stat = Stat()
                    # extension should be always .img in shika.db if generated from EIGER stream
                    possible_imgfs = (imgf, os.path.splitext(imgf)[0] + ".img",
                                      re.sub("(.*)_0([0-9]{6})\..*$", r"\1_\2.img", imgf), # too dirty fix!! for new bss which writes 7-digits filename..
                                      )
                    imgfs_found = [x for x in possible_imgfs if x in results]
                    if not imgfs_found: continue
                    imgf = imgfs_found[0]
                    snrlist = [x[2] for x in results[imgf]["spots"]]
                    stat.stats = (len(snrlist), sum(snrlist), numpy.median(snrlist) if snrlist else 0)
                    stat.spots = results[imgf]["spots"]
                    stat.gonio = gonio
                    stat.grid_coord = gc
                    stat.scan_info = scan
                    stat.thumb_posmag = results[imgf]["thumb_posmag"]
                    stat.params = results[imgf]["params"]
                    stat.img_file = os.path.join(self.ctrlFrame.current_target_dir, imgf)
                    result.append((stat.img_file, stat))

            delt = time.time() - startt
            shikalog.info("Data loaded: %s (took %f sec)" % (dbfile, delt))

            add_results(result)
        finally:
            d = None

        self.ctrlFrame.onResultsUpdated(EventResultsUpdated(result=result))
    # load_results()

    def track_latest_result(self):
        dic4t = self.ctrlFrame.dic_for_tree
        if self.ctrlFrame.current_target_dir is None:
            return

        key = tuple(os.path.relpath(self.ctrlFrame.current_target_dir, self.topdir).split(os.sep))
        if key in dic4t and self.ctrlFrame.treectrl.GetChildrenCount(dic4t[key]) > 0:
            lastchild = self.ctrlFrame.treectrl.GetLastChild(dic4t[key])
            if not self.ctrlFrame.treectrl.IsSelected(lastchild):
                self.ctrlFrame.treectrl.SelectItem(lastchild)
    # track_latest_result()
# class MainFrame

def run_from_args(argv):
    if "-h" in argv or "--help" in argv:
        print("All parameters:\n")
        iotbx.phil.parse(gui_phil_str).show(prefix="  ", attributes_level=1)
        return

    global gui_params
    cmdline = iotbx.phil.process_command_line(args=argv,
                                              master_string=gui_phil_str)
    gui_params = cmdline.work.extract()
    args = cmdline.remaining_args

    shikalog.config(gui_params.bl)
    shikalog.info("Program started in %s." % os.getcwd())

    topdir = None
    re_addr = re.compile("^[0-9]{,3}\.[0-9]{,3}\.[0-9]{,3}\.[0-9]{,3}:[0-9]+$")
    re_host = re.compile("^[0-9]{,3}\.[0-9]{,3}\.[0-9]{,3}\.[0-9]{,3}$")

    if gui_params.kuma_addr is not None:
        if not re_addr.search(gui_params.kuma_addr):
            shikalog.error("Invalid address definition of KUMA: %s" % gui_params.kuma_addr)
            return
        print("Config: KUMA addr=", gui_params.kuma_addr)


    if gui_params.mode == "zoo":
        gui_params.auto_mode = True

    if gui_params.imgview_host is not None:
        if not re_host.search(gui_params.imgview_host):
            shikalog.error("Invalid host definition of Imgview: %s" % gui_params.imgview_host)
            return

        shikalog.info("Config: imgview host= %s" % gui_params.imgview_host)

    if gui_params.subport:
        shikalog.info("Config: ZMQ SUB port= %s" % gui_params.subport)
        try: control_send.bind("tcp://*:%d"%gui_params.subport)
        except zmq.ZMQError as e:
            shikalog.error("Error in binding SUB port: %s" % e.strerror)
            print("If you don't need to change parameters, try subport=none")
            return
            
    if gui_params.pushport:
        shikalog.info("Config: ZMQ PUSH host= %s" % gui_params.pushport)
        try: ventilator_send.bind("tcp://*:%d"%gui_params.pushport)
        except zmq.ZMQError as e:
            shikalog.error("Error in binding PUSH port: %s" % e.strerror)
            print("If you don't need to recalculate, try pushport=none")
            return

    print("""\

SHIKA (Spot wo Hirotte Ichi wo Kimeru Application) is a spot finder application for diffraction based crystal centering based on Cheetah by Anton Barty et al.
If you found something wrong, please let staff know! We would appreciate your feedback.

""")

    if topdir is None:
        topdir = os.getcwd()

    app = wx.App()

    if gui_params.ask_directory:
        dlg = wx.DirDialog (None, "Choose directory to watch", "",
                            wx.DD_DEFAULT_STYLE | wx.DD_DIR_MUST_EXIST)
        if dlg.ShowModal() == wx.ID_OK:
            topdir = dlg.GetPath()
        else:
            return

    app.TopWindow = MainFrame(parent=None, id=wx.ID_ANY, 
                              topdir=topdir, params=gui_params)
    app.MainLoop()

    shikalog.info("Normal exit.")
# run_from_args()

if __name__ == "__main__":
    run_from_args(sys.argv[1:])
