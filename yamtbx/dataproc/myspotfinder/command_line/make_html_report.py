from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import os
import re
import time
import datetime
import collections
import glob
import sqlite3
import pickle
import numpy
import matplotlib
import matplotlib.figure
matplotlib.use('Agg') # Allow to work without X
from PIL import Image
import iotbx.phil
from yamtbx.util import rotate_file
from yamtbx.dataproc.myspotfinder import shikalog
from yamtbx.dataproc.myspotfinder.command_line.spot_finder_gui import Stat
from yamtbx.dataproc.dataset import re_pref_num_ext
from yamtbx.dataproc import bl_logfiles

master_params_str = """\
target_dir = None
 .type = path
rotate = False
 .type = bool
 .help = backup (rotate) old files
mode = *normal zoo
 .type = choice
plot = *grid circle
 .type = choice
"""

def plot_heatmap(subplot, xs, ys, ds, scaninfo):
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

    p1 = subplot.imshow(zi, origin='lower',
                        extent=[min(xs)-hs/2, max(xs)+hs/2,
                                min(ys)-vs/2, max(ys)+vs/2],
                        interpolation='none', cmap="YlOrRd")#PuRd

    return p1
# plot_heatmap()

def plot_circles(subplot, xs, ys, ds, zero_xs, zero_ys):
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

    p1 = subplot.scatter(xs, ys, s=normalize_max(ds), c=ds, alpha=0.5)
    p2 = subplot.scatter(zero_xs, zero_ys, s=50, marker="x", c=[0]*len(zero_xs), alpha=0.5)
    return p1, p2
# plot_circles()

def prepare_plot(plot_data, f, kind, wdir, rotate=False, plot_grid=True):
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
    for imgf, stat in plot_data[f]:
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
    #p = ax.scatter(xs, ys, s=normalize_max(ds), c=ds, alpha=0.5) # s in points^2

    scaninfo = plot_data[f][0][1].scan_info

    if plot_grid:
        p = plot_heatmap(ax, xs, ys, ds, scaninfo)
    else:
        p, _ = plot_circles(ax, xs, ys, ds, zero_xs, zero_ys)


    if max(ds) - min(ds) > 1e-5:
        fig.colorbar(p)
    ax.scatter(zero_xs, zero_ys, s=50, marker="x", c=[0]*len(zero_xs), alpha=0.5)
    ax.set_xlabel("horizontal [um]")
    ax.set_ylabel("vertical [um]")

    if scaninfo is not None:
        vp, hp = scaninfo.vpoints, scaninfo.hpoints
        vs, hs = scaninfo.vstep*1000., scaninfo.hstep*1000.

        if 1 in (vp, hp) or len(plot_data[f]) <= hp:
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
        vs, hs = 5, 5

    canvas = matplotlib.backends.backend_agg.FigureCanvasAgg(fig)
    canvas.print_figure(pngout+".tmp", dpi=80, format="png")
    img_width = fig.get_figwidth() * 80
    img_height = fig.get_figheight() * 80

    map_str = '<map name="%smap">\n' % scan_prefix
    for x, y, imgf in zip(xs, ys, imgfs):
        if plot_grid:
            tx1, ty1 = ax.transData.transform((x-hs/2.,y-vs/2.))
            tx2, ty2 = ax.transData.transform((x+hs/2.,y+vs/2.))
            map_str += '  <area shape="rect" coords="%.2f,%.2f,%.2f,%.2f" title="%s" onClick=\'plotClick("%s", "%s")\'>\n' % (tx1, img_height-ty1, tx2, img_height-ty2, os.path.basename(imgf), scan_prefix, os.path.basename(imgf))
        else:
            tx, ty = ax.transData.transform((x,y))
            map_str += '  <area shape="circle" coords="%.2f,%.2f,10" title="%s" onClick=\'plotClick("%s", "%s")\'>\n' % (tx, img_height-ty, os.path.basename(imgf), scan_prefix, os.path.basename(imgf))

    map_str += "</map>"
    return pngout, map_str
# prepare_plot()

def make_html_report(current_stats, wdir, htmlout, zoo_mode, rotate=False, plot_grid=True):
    #plot_data = self.plotFrame.data
    shikalog.info("Making HTML report for %s"%wdir)
    startt = time.time()

    plot_data = collections.OrderedDict()
    for f, stat in list(current_stats.items()):
        if stat is None: continue
        fpref = decide_fpref(f, stat.scan_info)
        plot_data.setdefault(fpref, []).append((f, stat))

    #if gui_params.mode == "zoo": htmlout = os.path.join(wdir, "report_zoo.html")
    #else: htmlout = os.path.join(wdir, "report.html")

    if rotate: rotate_file(htmlout)

    if zoo_mode: assert len(plot_data) <= 1

    kinds = ("total_integrated_signal", "median_integrated_signal", "n_spots")
    plots=""
    pngs = []
    for f in plot_data:
        scan_prefix = f[:f.index(" ")] if " (phi=" in f else f
        info = plot_data[f][0][1].scan_info

        if info is None: info = bl_logfiles.ScanInfo() # Empty info
        plots += '<table border=0 style="margin-bottom:0px">\n  <tr><td>\n'

        if zoo_mode:
            try:
                im = Image.open(os.path.join(wdir, "../../../before.ppm"))
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

        for i, kind in enumerate(kinds):
            pngout, mapstr = prepare_plot(plot_data, f, kind, wdir, rotate, plot_grid)
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
    dbspots = dict([(str(x[0]), pickle.loads(x[1])) for x in c.fetchall()])
    spot_data = "var spot_data = {"
    for i, (f, stat) in enumerate(result):
        if stat is None: continue
        bf = os.path.basename(f)
        spots = dbspots[bf]["spots"]
        thumb_posmag = dbspots[bf]["thumb_posmag"]
        r = re.search(r"^(.*)_([0-9]+)\.[^0-9]+$", bf)
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
            r = re.search(r"^(.*)_([0-9]+)\.[^0-9]+$", os.path.basename(res[0]))
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
# make_html_report()

def load_results(target_dir):
    current_stats = collections.OrderedDict()

    dbfile = os.path.join(target_dir, "_spotfinder", "shika.db")
    if not os.path.isfile(dbfile):
        shikalog.error("%s not found." % dbfile)
        return

    scanlog = os.path.join(target_dir, "diffscan.log")
    if not os.path.isfile(scanlog):
        shikalog.error("diffscan.log not found in %s" % target_dir)
        return

    slog = bl_logfiles.BssDiffscanLog(scanlog)
    slog.remove_overwritten_scans()

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
            results = dict([(str(x[0]), pickle.loads(x[1])) for x in c.fetchall()])
            break
        except sqlite3.DatabaseError:
            shikalog.warning("DB failed. retrying (%d)" % itrial)
            time.sleep(1)
            continue

    print("DEBUG:: scans=", slog.scans) 
    for scan in slog.scans:
        for imgf, (gonio, gc) in scan.filename_coords:
            #print imgf, (gonio, gc) 
            stat = Stat()
            # extension should be always .img in shika.db if generated from EIGER stream
            possible_imgfs = (imgf, os.path.splitext(imgf)[0] + ".img",
                              re.sub(r"(.*)_0([0-9]{6})\..*$", r"\1_\2.img", imgf), # too dirty fix!! for new bss which writes 7-digits filename..
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
            stat.img_file = os.path.join(target_dir, imgf)
            result.append((stat.img_file, stat))

    delt = time.time() - startt
    shikalog.info("Data loaded: %s (took %f sec)" % (dbfile, delt))

    for f, stat in result: current_stats[f] = stat
    return current_stats
# load_results()

def decide_fpref(f, scaninfo):
    fpref = re_pref_num_ext.search(os.path.basename(f)).group(1)
    if scaninfo is not None:
        if scaninfo.is_shutterless():
            fpref += " (phi=%.2f)" % (scaninfo.fixed_spindle)
        else:
            fpref += " (phi=%.2f)" % (scaninfo.osc_start)
    return fpref
# decide_fpref()

def run(params):
    wdir = os.path.abspath(params.target_dir)
    target_dir = os.path.normpath(os.path.join(wdir, ".."))
    current_stats = load_results(target_dir)

    zoo_mode = params.mode == "zoo"
    htmlout = os.path.join(wdir, "report_zoo.html" if zoo_mode else "report.html")
    make_html_report(current_stats, wdir, htmlout, zoo_mode, params.rotate, params.plot=="grid")
# run()

if __name__ == "__main__":
    shikalog.config(None)

    import sys
    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args

    if not params.target_dir and len(args) >= 1:
        params.target_dir = args[0]

    run(params)
