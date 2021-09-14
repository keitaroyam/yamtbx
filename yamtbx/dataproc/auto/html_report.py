"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import division
from __future__ import unicode_literals
import sys
import os
import time
import pickle
import iotbx.phil
import io

from yamtbx.dataproc.xds import correctlp
from yamtbx.dataproc.xds import idxreflp
from yamtbx.dataproc.xds import integratelp
from yamtbx.dataproc.xds import xdsstat
from yamtbx.dataproc.xds import get_xdsinp_keyword
from yamtbx import util

master_params_str = """
topdir = None
 .type = path
htmlout = None
 .type = path
"""

#amcharts_root = "/oys/xtal/yamtbx/yamtbx/dataproc/auto/js/amcharts_3.12.0/amcharts"
amcharts_root = "http://www.amcharts.com/lib/3"

def make_individual_report(xds_wd, html_wd):
    problems = []
    xds_inp = os.path.join(xds_wd, "XDS.INP")
    idxref_lp = os.path.join(xds_wd, "IDXREF.LP")
    xparm_xds = os.path.join(xds_wd, "XPARM.XDS")
    spot_xds = os.path.join(xds_wd, "SPOT.XDS")
    integrate_lp = os.path.join(xds_wd, "INTEGRATE.LP")
    correct_lp = os.path.join(xds_wd, "CORRECT.LP")
    xdsstat_lp = os.path.join(xds_wd, "XDSSTAT.LP")
    stats_pkl = os.path.join(xds_wd, "merging_stats.pkl")
    shika_log = os.path.join(xds_wd, "shika.log")

    if not os.path.exists(html_wd):
        os.makedirs(html_wd)

    html = """\
<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<script src="%(amcharts_root)s/amcharts.js" charset="utf-8"></script>
<script src="%(amcharts_root)s/serial.js" charset="utf-8"></script>
<script src="%(amcharts_root)s/xy.js" charset="utf-8"></script>
</head>
<body>
<h1>%(title)s</h1>
created on %(cdate)s
""" % dict(title=os.path.abspath(xds_wd), cdate=time.strftime("%Y-%m-%d %H:%M:%S"), amcharts_root=amcharts_root)

    if os.path.isfile(shika_log):
        lines = open(shika_log).readlines()
        if len(lines) > 1:
            data = [(x.split()[0], x.split()[2]) for x in lines[1:]]
            html += "<h1>Low resolution spots</h1>\n"
            html += """\
<div id="chartdiv_shika" style="width: 640px; height: 400px;"></div>
<script>
var chart_idxref_ibf = AmCharts.makeChart("chartdiv_shika", {
    "type": "serial",
    "theme": "none",
    "legend": {
        "useGraphSettings": true
    },
    "dataProvider": [%s],
    "valueAxes": [{
        "stackType": "regular",
        "gridAlpha": 0.07,
        "position": "left",
        "title": "Number of located spots by SHIKA"
    }],
    "graphs": [{
        "balloonText": "n_spots: [[value]]",
        "fillAlphas": 0.6,
        "lineAlpha": 0.4,
        "title": "number of spots",
        "valueField": "nspt"
    }],
    "chartCursor": {
        "cursorPosition": "mouse"
    },
    "categoryField": "frame",
    "categoryAxis": {
        "minorGridEnabled": true,
        "title": "frame"
    }
});
</script>
""" % (",".join(['{"frame": %s, "nspt": %s}'%x for x in data]),
       )

    html += "<h1>IDXREF (Indexing)</h1>\n"
    if os.path.isfile(idxref_lp):
        SX = idxreflp.SpotXds(spot_xds)
        if os.path.isfile(xparm_xds): SX.set_xparm(xparm_xds)
        else: SX.set_xdsinp(xds_inp)

        i_ui_on_det = SX.indexed_and_unindexed_on_detector(with_resolution=True)
        
        html += """\
<h3>Indexed/Unindexed spots</h3>
<h4>by frame</h4>
<div id="chartdiv_idxref_indexed_by_frame" style="width: 640px; height: 400px;"></div>
<script>
var chart_idxref_ibf = AmCharts.makeChart("chartdiv_idxref_indexed_by_frame", {
    "type": "serial",
    "theme": "none",
    "legend": {
        "useGraphSettings": true
    },
    "dataProvider": [%s],
    "valueAxes": [{
        "stackType": "regular",
        "gridAlpha": 0.07,
        "position": "left",
        "title": "Number of located spots"
    }],
    "graphs": [{
        "balloonText": "indexed: [[value]]",
        "fillAlphas": 0.6,
        "lineAlpha": 0.4,
        "title": "Indexed spots",
        "valueField": "ind"
    }, {
        "balloonText": "unindexed: [[value]]",
        "fillAlphas": 0.6,
        "lineAlpha": 0.4,
        "title": "Unindexed spots",
        "valueField": "und"
    }],
    "chartCursor": {
        "cursorPosition": "mouse"
    },
    "categoryField": "frame",
    "categoryAxis": {
        "minorGridEnabled": true,
        "title": "frame"
    }
});
</script>
<h4>on detector surface</h4>
<div id="chartdiv_idxref_indexed_on_det" style="width: 640px; height: 640px;"></div>
<script>
var chart_idxref_iod = AmCharts.makeChart("chartdiv_idxref_indexed_on_det", {
    "type": "xy",
    "theme": "none",
    "legend": {
        "useGraphSettings": true
    },
    "dataProvider": [%s],
    "valueAxes": [{
        "position":"bottom",
        "axisAlpha": 0,
        "dashLength": 1,
        "title": "Detector X"
    }, {
        "axisAlpha": 0,
        "dashLength": 1,
        "position": "left",
        "title": "Detector Y"
    }],
    "graphs": [{
        "balloonText": "x:[[x]] y:[[y]]; d=[[description]]",
        "bullet": "round",
        "bulletSize": 4,
        "lineAlpha": 0,
        "xField": "ix",
        "yField": "iy",
        "descriptionField": "d",
        "title": "Indexed spots",
        "hidden": true,
        "lineColor": "#FF6600",
		"fillAlphas": 0
    }, {
        "balloonText": "x:[[x]] y:[[y]]; d=[[description]]",
        "bullet": "round",
        "bulletSize": 4,
        "lineAlpha": 0,
        "xField": "ux",
        "yField": "uy",
        "descriptionField": "d",
        "title": "Unindexed spots",
        "hidden": true,
        "lineColor": "#FCD202",
		"fillAlphas": 0
    }],
    "chartCursor": {
        "cursorPosition": "mouse"
    },
});
</script>
""" %  (",".join(['{"frame": %d, "ind": %s, "und": %s}'%(x[0],x[1][0],x[1][1]) for x in SX.indexed_and_unindexed_by_frame()]),
        ",".join(['{"ix":%.1f,"iy":%.1f,"d":%.1f}'%x for x in i_ui_on_det["indexed"]] + ['{"ux":%.1f,"uy":%.1f,"d":%.1f}'%x for x in i_ui_on_det["unindexed"]]),
        )

        idexref_has_problem = False
        lp = idxreflp.IdxrefLp(idxref_lp)
        intness = lp.cluster_integerness()

        for h in intness:
            if sum(h) > 0 and h[0]/sum(h) < 0.9:
                idexref_has_problem = True
        
        html += "<h3>Clusters and subtrees</h3>"
        html += """
<h4>Cluster integerness</h4>
<div id="chartdiv_idxref_cluster_int" style="width: 640px; height: 400px;"></div>
<script>
var chart_idxref_cluster_int = AmCharts.makeChart("chartdiv_idxref_cluster_int", {
	"type": "serial",
	"categoryField": "f",
	"categoryAxis": {
		"gridPosition": "start",
		"position": "left",
		"title": "fractional part of indices",
	},
    "legend": {
        "useGraphSettings": true
    },
	"dataProvider": [%s],
	"graphs": [{
       		"balloonText": "h: [[value]]",
		"fillAlphas": 0.8,
		"lineAlpha": 0.2,
		"title": "h",
		"type": "column",
		"valueField": "h"
		}, {
       		"balloonText": "k: [[value]]",
		"fillAlphas": 0.8,
		"lineAlpha": 0.2,
		"title": "k",
		"type": "column",
		"valueField": "k"
		}, {
       		"balloonText": "l: [[value]]",
		"fillAlphas": 0.8,
		"lineAlpha": 0.2,
		"title": "l",
		"type": "column",
		"valueField": "l"
		}
	],
	"guides": [],
	"valueAxes": [
		{
			"position": "top",
			"axisAlpha": 0,
			"title": "frequency"
		}
	],
});
</script>
""" % (",".join(['{"f":.%d,"h":%d,"k":%d,"l":%d}'%x for x in [tuple([x,]+[intness[y][x] for y in range(3)]) for x in range(6)]])
       )

        html += "<h4>Subtree populations</h4>\n<pre>"
        accum_pop = 0
        for i, pp in enumerate([float(x)/sum(lp.subtree_population) for x in lp.subtree_population]):
            html += "%s subtree: %.2f%%\n" % (util.num_th_str(i+1), pp*100.)
            if i==0 and pp < 0.9: idexref_has_problem = True

            accum_pop += pp
            if accum_pop > 0.9:
                break
        html += "</pre>"

    if not os.path.isfile(xparm_xds): idexref_has_problem = True

    if idexref_has_problem:
        problems.append("IDXREF")

    html += "<h1>INTEGRATE (Integration)</h1>\n"
    if os.path.isfile(integrate_lp):
        lp = integratelp.IntegrateLp(integrate_lp)

        if any([float(x) > 1.5 for x in lp.sigmars]):
            problems.append("INTEGRATE")
        
        html += """\
<h3>Mosaicity and beam divergence</h3>
<div id="chartdiv_integrate_sigmas" style="width: 640px; height: 400px;"></div>

<h3>rotations off from initial orientation</h3>
<div id="chartdiv_integrate_rotxyz" style="width: 640px; height: 400px;"></div>
<script>
var chart1 = AmCharts.makeChart("chartdiv_integrate_sigmas", {
    "type": "serial",
    "theme": "none",
    "legend": {
        "useGraphSettings": true
    },
    "dataProvider": [%s],
    "valueAxes": [{
        "id":"v1",
        "axisColor": "#FF6600",
        "axisThickness": 2,
        "gridAlpha": 0,
        "axisAlpha": 1,
        "position": "left"
    }, {
        "id":"v2",
        "axisColor": "#FCD202",
        "axisThickness": 2,
        "gridAlpha": 0,
        "axisAlpha": 1,
        "position": "right"
    }],
    "graphs": [{
        "balloonText": "&sigma;<sub>R</sub>: [[value]]",
        "valueAxis": "v1",
        "lineColor": "#FF6600",
        "bullet": "round",
        "bulletBorderThickness": 1,
        "hideBulletsCount": 30,
        "title": "mosaicity",
        "valueField": "sigmar",
		"fillAlphas": 0
    }, {
        "balloonText": "&sigma;<sub>B</sub>: [[value]]",
        "valueAxis": "v2",
        "lineColor": "#FCD202",
        "bullet": "square",
        "bulletBorderThickness": 1,
        "hideBulletsCount": 30,
        "title": "beam divergence",
        "valueField": "sigmad",
		"fillAlphas": 0
    }],
    "chartCursor": {
        "cursorPosition": "mouse"
    },
    "categoryField": "frame",
    "categoryAxis": {
        "minorGridEnabled": true,
        "title": "frame"
    }
});

var chart2 = AmCharts.makeChart("chartdiv_integrate_rotxyz", {
    "type": "serial",
    "theme": "none",
    "legend": {
        "useGraphSettings": true
    },
    "dataProvider": [%s],
    "valueAxes": [{
        "id":"v1",
        "axisColor": "#FF6600",
        "axisThickness": 2,
        "gridAlpha": 0,
        "axisAlpha": 1,
        "position": "left"
    }],
    "graphs": [{
        "valueAxis": "v1",
        "lineColor": "#FF6600",
        "bullet": "round",
        "bulletBorderThickness": 1,
        "hideBulletsCount": 30,
        "title": "rot x",
        "valueField": "rotx",
		"fillAlphas": 0
    }, {
        "valueAxis": "v2",
        "lineColor": "#FCD202",
        "bullet": "square",
        "bulletBorderThickness": 1,
        "hideBulletsCount": 30,
        "title": "rot y",
        "valueField": "roty",
		"fillAlphas": 0
    }, {
        "valueAxis": "v3",
        "lineColor": "#B0DE09",
        "bullet": "triangleUp",
        "bulletBorderThickness": 1,
        "hideBulletsCount": 30,
        "title": "rot z",
        "valueField": "rotz",
		"fillAlphas": 0
    }],
    "chartCursor": {
        "cursorPosition": "mouse"
    },
    "categoryField": "frame",
    "categoryAxis": {
        "minorGridEnabled": true,
        "title": "frame"
    }
});

</script>
""" % (",".join(['{"frame": %d, "sigmar": %s, "sigmad": %s}'%x for x in zip(lp.frames,lp.sigmars,lp.sigmads)]),
       ",".join([",".join(['{"frame": %d, "rotx": %s, "roty": %s, "rotz": %s}'%((y,)+tuple(x[1].get("rotation", ["NaN"]*3))) for y in x[0]]) for x in list(lp.blockparams.items())]),
       )

    html += "<h1>CORRECT (Scaling)</h1>\n"
    if os.path.isfile(correct_lp):
        lp = correctlp.CorrectLp(correct_lp)

        if not lp.is_ISa_valid() or lp.get_ISa() < 10:
            problems.append("CORRECT")

        html += """\
<table border=0 cellpadding=5 style="text-align:right;">
<tr><th>Space group</th><td colspan=6 style="text-align:left;">%(sg)s</td></tr>
""" % dict(#sgstr=lp.space_group.info(), sgnum=lp.space_group.type().number(),
            sg=lp.space_group.info().symbol_and_number() if lp.space_group is not None else "n/a")

        if lp.unit_cell is not None:
            html += """\
<tr><th>Unit cell</th><td>%(a).3f</td><td>%(b).3f</td><td>%(c).3f</td><td>%(alpha).3f</td><td>%(beta).3f</td><td>%(gamma).3f</td></tr>
<tr><th>e. s. d.</th><td>%(a_e).3f</td><td>%(b_e).3f</td><td>%(c_e).3f</td><td>%(alpha_e).3f</td><td>%(beta_e).3f</td><td>%(gamma_e).3f</td></tr>
</table>
""" % dict(a=lp.unit_cell[0], b=lp.unit_cell[1], c=lp.unit_cell[2], 
           alpha=lp.unit_cell[3], beta=lp.unit_cell[4], gamma=lp.unit_cell[5],
           a_e=lp.unit_cell_esd[0], b_e=lp.unit_cell_esd[1], c_e=lp.unit_cell_esd[2],
           alpha_e=lp.unit_cell_esd[3], beta_e=lp.unit_cell_esd[4], gamma_e=lp.unit_cell_esd[5],
           )
        else: # in case parameter refinement failed
            html += """\
<tr><th>Unit cell</th><td colspan=6>n/a</td></tr>
<tr><th>e. s. d.</th><td colspan=6>n/a</td></tr>
</table>
"""

        html += "<h3>ISa</h3>\n<pre>%s</pre>\n" % lp.snippets.get("ISa","")
        html += "<h3>Statistics of all data</h3>\n"
        if os.path.isfile(stats_pkl):
            tmp = pickle.load(open(stats_pkl, "rb"))
            if "stats" in tmp:
                sio = io.StringIO()
                tmp["stats"].show(out=sio, header=False)
                stats_str = sio.getvalue()
            else:
                stats_str = tmp["stats_str"]
            html += "<pre>%s</pre>\n" % stats_str.replace("<","&lt;").replace(">","&gt;")
        else:
            html += "<pre>%s</pre>\n" % lp.snippets.get("table1","")

        if lp.error_table != {}:
            tmp = ['{"d":%.3f,"ios":%.2f,"nrej":%d}'%x for x in zip(lp.error_table["dmin"], 
                                                                         lp.error_table["ios"], 
                                                                         lp.error_table["nrej"])]
            tmp = tmp[:-1] # last is about all data
            html += """\
<h3>By resolution</h3>
<div id="chartdiv_correct_error" style="width: 640px; height: 400px;"></div>
<script>
var chart_correct_er = AmCharts.makeChart("chartdiv_correct_error", {
    "type": "serial",
    "theme": "none",
    "legend": {
        "useGraphSettings": true
    },
    "dataProvider": [%s],
    "valueAxes": [{
        "id":"v1",
        "axisColor": "#FF6600",
        "axisThickness": 2,
        "gridAlpha": 0,
        "axisAlpha": 1,
        "position": "left"
    }, {
        "id":"v2",
        "axisColor": "#FCD202",
        "axisThickness": 2,
        "gridAlpha": 0,
        "axisAlpha": 1,
        "position": "right"
    }],
    "graphs": [{
        "balloonText": "Mn(I/sd): [[value]]",
        "valueAxis": "v1",
        "lineColor": "#FF6600",
        "bullet": "round",
        "bulletBorderThickness": 1,
        "hideBulletsCount": 30,
        "title": "Mean(I/Sigma)",
        "valueField": "ios",
		"fillAlphas": 0
    }, {
        "balloonText": "rejected: [[value]]",
        "valueAxis": "v2",
        "lineColor": "#FCD202",
        "bullet": "square",
        "bulletBorderThickness": 1,
        "hideBulletsCount": 30,
        "title": "Rejected reflections",
        "valueField": "nrej",
		"fillAlphas": 0
    }],
    "chartCursor": {
        "cursorPosition": "mouse"
    },
    "categoryField": "d",
    "categoryAxis": {
        "minorGridEnabled": true,
        "title": "resolution"
    }
});
</script>
""" % (",".join(tmp))

    if os.path.isfile(xdsstat_lp):
        lp = xdsstat.XdsstatLp(xdsstat_lp)
        if set(("frame","iobs","sigma","ios","rmeas","nmisfits")).issubset(list(lp.by_frame.keys())):
            html += """\
<h3>By frame</h3>
<div id="chartdiv_xdsstat_ios" style="width: 640px; height: 400px;"></div>
<script>
var chart_xdsstat_ios = AmCharts.makeChart("chartdiv_xdsstat_ios", {
    "type": "serial",
    "theme": "none",
    "legend": {
        "useGraphSettings": true
    },
    "dataProvider": [%s],
    "valueAxes": [{
        "id":"v1",
        "axisThickness": 2,
        "gridAlpha": 0,
        "axisAlpha": 1,
        "position": "left"
    }, {
        "id":"v2",
        "axisThickness": 2,
        "gridAlpha": 0,
        "axisAlpha": 1,
        "position": "right"
    }],
    "graphs": [{
        "balloonText": "Mn(I): [[value]]",
        "valueAxis": "v1",
        "bullet": "round",
        "bulletBorderThickness": 1,
        "hideBulletsCount": 30,
        "title": "Mean(I)",
        "valueField": "i",
		"fillAlphas": 0
    }, {
        "balloonText": "Mn(sd): [[value]]",
        "valueAxis": "v1",
        "bullet": "round",
        "bulletBorderThickness": 1,
        "hideBulletsCount": 30,
        "title": "Mean(Sigma)",
        "valueField": "sig",
		"fillAlphas": 0
    }, {
        "balloonText": "Mn(I/sd): [[value]]",
        "valueAxis": "v2",
        "bullet": "square",
        "bulletBorderThickness": 1,
        "hideBulletsCount": 30,
        "title": "Mean(I/Sigma)",
        "hidden": true,
        "valueField": "ios",
		"fillAlphas": 0
    }],
    "chartCursor": {
        "cursorPosition": "mouse"
    },
    "categoryField": "frame",
    "categoryAxis": {
        "minorGridEnabled": true,
        "title": "frame"
    }
});
</script>

<div id="chartdiv_xdsstat_r" style="width: 640px; height: 400px;"></div>
<script>
var chart_xdsstat_r = AmCharts.makeChart("chartdiv_xdsstat_r", {
    "type": "serial",
    "theme": "none",
    "legend": {
        "useGraphSettings": true
    },
    "dataProvider": [%s],
    "valueAxes": [{
        "id":"v1",
        "axisColor": "#FF6600",
        "axisThickness": 2,
        "gridAlpha": 0,
        "axisAlpha": 1,
        "position": "left"
    }, {
        "id":"v2",
        "axisThickness": 2,
        "axisColor": "#FCD202",
        "gridAlpha": 0,
        "axisAlpha": 1,
        "position": "right"
    }],
    "graphs": [{
        "balloonText": "R-meas: [[value]]",
        "valueAxis": "v1",
        "lineColor": "#FF6600",
        "bullet": "round",
        "bulletBorderThickness": 1,
        "hideBulletsCount": 30,
        "title": "R-meas",
        "valueField": "r",
		"fillAlphas": 0
    }, {
        "balloonText": "Misfits: [[value]]",
        "valueAxis": "v2",
        "lineColor": "#FCD202",
        "bullet": "square",
        "bulletBorderThickness": 1,
        "hideBulletsCount": 30,
        "title": "#Misfits",
        "valueField": "mis",
		"fillAlphas": 0
    }],
    "chartCursor": {
        "cursorPosition": "mouse"
    },
    "categoryField": "frame",
    "categoryAxis": {
        "minorGridEnabled": true,
        "title": "frame"
    }
});
</script>
""" % (",".join(['{"frame":%d,"i":%.2f,"sig":%.2f,"ios":%.2f}'%x for x in zip(lp.by_frame["frame"], lp.by_frame["iobs"], lp.by_frame["sigma"], lp.by_frame["ios"])]),
       ",".join(['{"frame":%d,"r":%.4f,"mis":%d}'%x for x in zip(lp.by_frame["frame"], lp.by_frame["rmeas"], lp.by_frame["nmisfits"])])
       )

    html += "</body></html>"
    htmlout = os.path.join(html_wd, "report.html")
    open(htmlout, "w").write(html)
    return htmlout, problems
# make_individual_report()

def find_problems(xds_wd):
    problems = []
    idxref_lp = os.path.join(xds_wd, "IDXREF.LP")
    xparm_xds = os.path.join(xds_wd, "XPARM.XDS")
    integrate_lp = os.path.join(xds_wd, "INTEGRATE.LP")
    correct_lp = os.path.join(xds_wd, "CORRECT.LP")
    xdsstat_lp = os.path.join(xds_wd, "XDSSTAT.LP")

    if os.path.isfile(idxref_lp):
        idexref_has_problem = False
        lp = idxreflp.IdxrefLp(idxref_lp)
        intness = lp.cluster_integerness()

        for h in intness:
            if sum(h) > 0 and h[0]/sum(h) < 0.9:
                idexref_has_problem = True
        
        if len(lp.subtree_population) > 0:
            pp0 = float(lp.subtree_population[0])/sum(lp.subtree_population)
            if pp0 < 0.9: idexref_has_problem = True
    else:
        idexref_has_problem = True

    if not os.path.isfile(xparm_xds): idexref_has_problem = True

    if idexref_has_problem:
        problems.append("IDXREF")

    if os.path.isfile(integrate_lp):
        lp = integratelp.IntegrateLp(integrate_lp)

        if any([float(x) > 1.5 for x in lp.sigmars]):
            problems.append("INTEGRATE")
        
    if os.path.isfile(correct_lp):
        lp = correctlp.CorrectLp(correct_lp)

        if not lp.is_ISa_valid() or lp.get_ISa() < 10:
            problems.append("CORRECT")

    if os.path.isfile(xdsstat_lp):
        lp = xdsstat.XdsstatLp(xdsstat_lp)

    return problems
# find_problems()

def make_kamo_report(bssjobs, topdir, htmlout):
    report_html = """\
<html>
<head>
<style>
.dataset_table {
    font-family: "Trebuchet MS", Arial, Helvetica, sans-serif;
    width: 100%%;
    border-collapse: collapse;
}

.dataset_table td, .dataset_table th {
    font-size: 1em;
    border: 1px solid #98bf21;
    padding: 3px 7px 2px 7px;
}

.dataset_table th {
    font-size: 1.1em;
    text-align: left;
    padding-top: 5px;
    padding-bottom: 4px;
    background-color: #A7C942;
    color: #ffffff;
}

.dataset_table tr.alt td {
    color: #000000;
    background-color: #EAF2D3;
}

h1, h2 { text-align: center; text-indent: 0px; font-weight: bold; hyphenate: none;  
      page-break-before: always; page-break-inside: avoid; page-break-after: avoid; }
h3, h4, h5, h6 { text-indent: 0px; font-weight: bold; 
      hyphenate:  none; page-break-inside: avoid; page-break-after: avoid; }

</style>

</head>
<body>
<h1>KAMO (auto data processing) report</h1>
<div align="right">
top dir: %s<br>
created on %s
</div>
""" % (os.path.abspath(topdir), time.strftime("%Y-%m-%d %H:%M:%S"))

    reports = {}

    for jobkey in sorted(bssjobs.jobs):
        wd = bssjobs.get_xds_workdir(jobkey)
        relwd = os.path.relpath(wd, os.path.dirname(htmlout))
        job = bssjobs.jobs[jobkey]

        dsname = os.path.basename(wd).replace("xds_", "")
        correct_lp = os.path.join(wd, "CORRECT.LP")
        spot_xds = os.path.join(wd, "SPOT.XDS")
        idxref_lp = os.path.join(wd, "IDXREF.LP")
        xds_inp = os.path.join(wd, "XDS.INP")
        stats_pkl = os.path.join(wd, "merging_stats.pkl")

        sampleid = "%s(%.2d)" % job.sample if job.sample is not None else "?"
        wavelen = job.wavelength
        deltaphi = job.osc_step

        lattp, ISa, resn, cmpl, sg = float("nan"), float("nan"), float("nan"), float("nan"), "?"

        problems = find_problems(wd)
        indiv_html = os.path.join(relwd, "report.html")
        problems_str = "OK"
        if len(problems) > 0:
            problems_str = ",".join(problems)

        totalphi = job.osc_end - job.osc_start

        if os.path.isfile(idxref_lp):
            lpobj = idxreflp.IdxrefLp(idxref_lp)
            lattp = lpobj.first_subtree_fraction()*100.

        lp = bssjobs._load_if_chached("correctlp", correct_lp)
        if lp is not None:
            ISa = lp.get_ISa() if lp.is_ISa_valid() else float("nan")
            sg = lp.space_group_str()
            cmpl = float(lp.table["all"]["cmpl"][-1]) if "all" in lp.table else float("nan")
        
        resn = bssjobs._load_if_chached("resn", correct_lp)
        if resn is None: resn = bssjobs._load_if_chached("resn", spot_xds)
        if resn is None: resn = float("nan")

        tmp = '<tr>\n <td><a href="%s">%s</a></td> <td>%s</td> <td>%.4f</td> <td>%.1f</td> <td>%.3f</td> <td>%.1f</td> <td>%s</td> <td>%.1f</td> <td>%.0f</td> <td>%.2f</td>  <td>%s</td>\n</tr>\n' % (indiv_html, dsname, sampleid, wavelen, totalphi, deltaphi, lattp, sg,
                                                                                                                                                               resn, cmpl, ISa, problems_str)
        reports.setdefault(os.path.dirname(relwd), []).append(tmp)

    for wd in sorted(reports):
        report_html += "<h3>%s</h3>\n" % wd
        report_html += """\
<table class="dataset_table">
<tr>
 <th>Dataset</th> <th>Sample ID</th> <th>&lambda; (&Aring;)</th> <th>Total &phi; (&deg;)</th> <th>&Delta;&phi; (&deg;)</th> <th>IdxLatt (%)</th> <th>SG</th> <th>Resn (&Aring;)</th> <th>Cmpl (%)</th> <th>ISa</th> <th>Problems</th>
</tr>
"""
        report_html += "".join(reports[wd])
        report_html += "</table>\n"

    report_html += """
</body>
</html>
"""

    open(htmlout, "w").write(report_html)

# run()
