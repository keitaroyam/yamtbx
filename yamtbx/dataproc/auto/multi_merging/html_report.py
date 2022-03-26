"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import unicode_literals
from yamtbx.dataproc.auto.html_report import amcharts_root
from yamtbx import util
from yamtbx.util.xtal import CellConstraints
from yamtbx.dataproc.xds import xscalelp

import libtbx.phil
import libtbx.load_env
from cctbx.array_family import flex
from cctbx import sgtbx

import time
import os
import shutil

class HtmlReportMulti(object):
    def __init__(self, root):
        self.root = root

        jsdir = os.path.join(self.root, "js")
        d3path = os.path.join(util.yamtbx_module_root(), "dataproc/auto/js/d3-3.5.10")
        if not d3path:
            raise Exception("Can't find d3-3.5.10 location. Check installation of yamtbx.")

        shutil.copytree(d3path, jsdir)

        self.html_head = """\
<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<script src="%(amcharts_root)s/amcharts.js" charset="utf-8"></script>
<script src="%(amcharts_root)s/serial.js" charset="utf-8"></script>
<script src="%(amcharts_root)s/xy.js" charset="utf-8"></script>
<script src="js/d3.min.js"></script>
<script src="js/index.js"></script>

<script>
  var toggle_show = function(obj_id) {
    var trg = document.getElementById(obj_id);
    if (trg.style.display === 'block' || trg.style.display === '')
      trg.style.display = 'none';
    else
      trg.style.display = 'block';
  }

  var toggle_show2 = function(caller, obj_id) {
    var trg = document.getElementById(obj_id);
    if (trg.style.display === 'block' || trg.style.display === '') {
      trg.style.display = 'none';
      trg.style.padding = '0px';
      caller.innerHTML= '&#x25bc;';
    } else {
      trg.style.display = '';
      trg.style.padding = '7px';
      caller.innerHTML= '&#x25b2;';
    }
  }

  var show_or_hide = function(caller_id, obj_id, flag) { 
    var trg = document.getElementById(obj_id);
    if (flag) { 
      trg.style.display = 'none';  
      trg.style.padding = '0px';   
      document.getElementById(caller_id).innerHTML= '&#x25bc;';   
    } else {    
      trg.style.display = '';      
      trg.style.padding = '7px';   
      document.getElementById(caller_id).innerHTML= '&#x25b2;';   
    }
  }
</script>

<style>
pre {
    font-family: Consolas, 'Courier New', Courier, Monaco, monospace;
    font-size: 1.1em;
}

.cells_table {
    font-family: "Trebuchet MS", Arial, Helvetica, sans-serif;
    width: 100%%;
    border-collapse: collapse;
}

.cells td, .dataset_table th,
.merge td {
    font-size: 1em;
    #border: 1px solid #98bf21;
    padding: 4px 7px 4px 7px;
}

.cells th,
.merge th {
    font-size: 1.1em;
    text-align: center;
    padding: 4px;
    background-color: #A7C942;
    color: #ffffff;
}

/*
.cells tr.alt td {
    color: #000000;
    background-color: #EAF2D3;
}
*/

.merge tr:nth-child(4n+3),
.merge tr:nth-child(4n),
.cells tr:nth-child(odd) {
    color: #000000;
    background-color: #EAF2D3;
}
.merge tr:nth-child(4n+1),
.merge tr:nth-child(4n+2),
.cells tr:nth-child(even) {
    #color: #f8fbf1;
    background-color: #f8fbf1;
}

.node circle {
  fill: #fff;
  stroke: steelblue;
  stroke-width: 1.5px;
}

.node {
 font: 15px sans-serif;
}

.link {
  fill: none;
  stroke: #ccc;
  stroke-width: 1.5px;
}

.d3-tip {
  line-height: 1;
  font-weight: bold;
  padding: 12px;
  background: rgba(0, 0, 0, 0.8);
  color: #fff;
  border-radius: 2px;
}

.node--source {
  fill: #2ca02c;
}

.node--target {
  fill: #d62728;
}

.link--source,
.link--target {
  stroke-opacity: 1;
  stroke-width: 2px;
}

.link--source {
  stroke: #d62728;
}

.link--target {
  stroke: #2ca02c;
}

/* Histogram */
.bar rect {
  fill: steelblue; 
  shape-rendering: crispEdges;
}

.bar text {
  font: 10px sans-serif;
  fill: #fff;
}

.axis path, .axis line {
  fill: none;
  stroke: #000;
  shape-rendering: crispEdges;
}

</style>
</head>
<body>
<h1>KAMO.MULTI_MERGE report</h1>
<div align="right">
workdir: %(wd)s<br />
created on %(cdate)s
</div>
<hr>
""" % dict(wd=self.root, cdate=time.strftime("%Y-%m-%d %H:%M:%S"), amcharts_root=amcharts_root)

        self.html_params = ""
        self.html_inputfiles = ""
        self.html_clustering = ""
        self._html_clustering_details = ""
        self.html_merge_results = []
        self.html_merge_plot_data = []
        self.params = None
        self.cells = None
    # __init__()

    def add_params(self, params, master_params_str):
        self.params = params
        self.html_params = """
<h2>Parameters</h2>
<a href="#" onClick="toggle_show('div-params'); return false;">Show/Hide</a>
<div id="div-params" style="display:none;">
<pre>
%s
</pre>
</div>
""" % libtbx.phil.parse(master_params_str).format(params).as_str(prefix=" ")
        self.write_html()
    # add_params()

    def add_cells_and_files(self, cells, symm_str):
        self.cells = cells
        # Table
        table_str = ""
        for idx, xac in enumerate(cells):
            cell = cells[xac]
            table_str += "<tr>\n"
            table_str += " <td>%.4d</td><td>%s</td>" % (idx+1, xac) # idx, file
            table_str += "".join(["<td>%.2f</td>"%x for x in cell])
            table_str += "\n</tr>\n"

        # Hist
        cellconstr = CellConstraints(sgtbx.space_group_info(symm_str).group())
        show_flags = (True, not cellconstr.is_b_equal_a(), not cellconstr.is_c_equal_a_b(),
                      not cellconstr.is_angle_constrained("alpha"),
                      not cellconstr.is_angle_constrained("beta"),
                      not cellconstr.is_angle_constrained("gamma"))
        names = ("a", "b", "c", "&alpha;", "&beta;", "&gamma;")

        hist_str = ""
        label1 = ""
        for i, (name, show) in enumerate(zip(names, show_flags)):
            tmp = ""
            if i in (0,3): tmp += "<tr>"
            if show: tmp += "<th>%s</th>" % name
            if i in (2,5): tmp += "</tr>"

            if i < 3: hist_str += tmp
            else: label1 += tmp

        hist_str += "\n<tr>\n"

        for idx, (name, show) in enumerate(zip(names, show_flags)):
            if idx==3: hist_str += "</tr>" + label1 + "<tr>"
            if not show: continue
            vals = flex.double([x[idx] for x in list(cells.values())])
            if len(vals) == 0: continue
            nslots = max(30, int((max(vals) - min(vals)) / 0.5))
            hist = flex.histogram(vals, n_slots=nslots)
            x_vals = [hist.data_min() + hist.slot_width() * (i+.5) for i in range(len(hist.slots()))]
            y_vals = hist.slots()
            hist_str += """
<td>
<div id="chartdiv_cell%(idx)d" style="width: 500px; height: 400px;"></div>
<script>
 var chart_cell%(idx)d = AmCharts.makeChart("chartdiv_cell%(idx)d", {
    "type": "serial",
    "theme": "none",  
    "legend": {
        "useGraphSettings": true,
        "markerSize":12,
        "valueWidth":0,
        "verticalGap":0
    },
    "dataProvider": [%(data)s],
    "valueAxes": [{
        "minorGridAlpha": 0.08,
        "minorGridEnabled": true,
        "position": "top",
        "axisAlpha":0
    }],
    "graphs": [{
        "balloonText": "[[category]]: [[value]]",
        "title": "%(name)s",
        "type": "column",
        "fillAlphas": 0.8,
        "valueField": "yval"
    }],
    "rotate": false,
    "categoryField": "xval",
    "categoryAxis": {
        "gridPosition": "start",
        "title": ""
    }
});
</script>
</td>
""" % dict(idx=idx, name=name,
           data=",".join(['{"xval":%.2f,"yval":%d}'%x for x in zip(x_vals,y_vals)])
           )

        hist_str += "</tr>"

        self.html_inputfiles = """
<h2>Input files</h2>
%d files for merging in %s symmetry

<h3>Unit cell histogram</h3>
<table>
%s
</table>

<h3>Files</h3>
<a href="#" onClick="toggle_show('div-input-files'); return false;">Show/Hide</a>
<div id="div-input-files" style="display:none;">
<table class="cells">
<tr>
 <th>idx</th> <th>file</th> <th>a</th> <th>b</th> <th>c</th> <th>&alpha;</th> <th>&beta;</th> <th>&gamma;</th>
</tr>
%s
</table>
</div>
""" % (len(cells), symm_str, hist_str, table_str)
        self.write_html()
    # add_files()

    def add_clutering_result(self, clusters, method):
        assert method in ("cc_clustering", "blend", "cumulative")

        if method == "cumulative":
            return
        elif method == "cc_clustering":
            clsdat = os.path.join(self.root, "cc_clustering", "cc_cluster_summary.dat")
        else:
            clsdat = os.path.join(self.root, "blend", "blend_cluster_summary.dat")
            
        header = None
        data = ""
        for l in open(clsdat):
            if header is None and not l.startswith("#"):
                header = "".join(["<th>%s</th>" % x for x in l.split()])
            elif header is not None:
                data += "<tr>%s</tr>" % "".join(["<td>%s</td>" % x for x in l.split()])

        treejson = os.path.join(self.root, method, "dendro.json")
        treedata = open(treejson).read() if os.path.isfile(treejson) else ""

        cluster_descr, file_descr = [], []
        IDs_set = set()
        for vv in clusters:
            if method == "cc_clustering":
                clno, IDs, clh, cmpl, redun, acmpl, aredun, ccmean, ccmin = vv
                cluster_descr.append('"%s":"Cluster_%.4d (%d files)<br />ClH: %5.2f<br />Cmpl= %5.1f%%, Redun=%5.1f<br />ACmpl= %5.1f%%, ARedun=%5.1f<br />CCmean=%.4f, CCmin=%.4f"' % (clno, clno, len(IDs), clh, cmpl, redun, acmpl, aredun, ccmean, ccmin))
            else:
                clno, IDs, clh, cmpl, redun, acmpl, aredun, LCV, aLCV = vv
                cluster_descr.append('"%s":"Cluster_%.4d (%d files)<br />ClH: %5.2f, LCV: %5.2f%%, aLCV: %5.2f &Aring;<br />Cmpl= %5.1f%%, Redun=%5.1f<br />ACmpl= %5.1f%%, ARedun=%5.1f"' % (clno, clno, len(IDs), clh, LCV, aLCV, cmpl, redun, acmpl, aredun))

            IDs_set.update(IDs)

        xac_files = list(self.cells.keys())
        for idx in IDs_set:
            file_descr.append('%s:"%s"' % (idx, xac_files[idx-1]))

        if method == "cc_clustering":
            self.html_clustering = "<h2>CC-based clustering</h2>"
        else:
            self.html_clustering = "<h2>Cell-based clustering by BLEND</h2>"

        self.html_clustering += self._html_clustering_details

        self.html_clustering += """
<a href="%(method)s/tree.png">See original cluster dendrogram</a>
<div id="tree-svg-div">
<script type="application/json" id="treedata">%(treedata)s</script>
<script>
  
  var data = document.getElementById('treedata').innerHTML;
  root = JSON.parse(data);
  
  var merged_clusters = [%(merged_clusters)s];
  var cluster_descr = {%(cluster_descr)s};
  var file_descr = {%(file_descr)s};

  var width = 1500, height = 600;
  var cluster = d3.layout.cluster()
      .size([width-20, height-40]);

  var diagonal = d3.svg.diagonal()
      .projection(function(d) { return [d.x, d.y]; });

  var svg = d3.select("#tree-svg-div").append("svg")
    .attr("width", width)
    .attr("height", height)
    .append("g")
    .attr("transform", "translate(0,20)");
  
  var nodes = cluster.nodes(root),
  links = cluster.links(nodes);
  var link = svg.selectAll(".link")
  .data(links)
  .enter().append("path")
  .attr("class", "link")
  .attr("d", diagonal);
  
  var node = svg.selectAll(".node")
  .data(nodes)
  .enter().append("g")
  .attr("class", "node")
  .on("mouseover", mouseovered)
  .on("mouseout", mouseouted)
  .attr("transform", function(d) { return "translate(" + d.x + "," + d.y + ")"; });

  node.append("circle")
  .attr("r", function(d) { return d.children ? 4.5 : 3; });

  var tip = d3.tip()
   .attr('class', 'd3-tip')
   .offset([-10, 0]);
  svg.call(tip);

  node.append("text")
   .attr("dx", function(d) { return d.children ? 12 : 0; })
   .attr("dy", function(d) { return d.children ? 3 : 12; })
   .attr("text-anchor", "middle")
   .style("font-size", function(d) { return d.children ? "14px" : "10px"; })
   .style("font-weight", function(d) { return d.children ? "bold" : ""; })
   .style("fill", function(d) { return d.children && merged_clusters.indexOf(d.name)!=-1 ? "red" : ""; })
   .text(function(d) { return d.name; });
  
// http://bl.ocks.org/mbostock/7607999
function mouseovered(d) {
  node
      .each(function(n) { n.target = n.source = false; });
  link
     .classed("link--source", function(l) { if (l.source === d) return l.target.target = true; });
  if (d.children)
     d.children.forEach(function(n) {
        link.classed("link--source", function(n) { if (n.source.target || n.target.target) return n.target.target = true; }); node.classed("node--source", function(n) { if (!n.children && n.target) return true; }) 
     });
  if (d.children) tip.html(cluster_descr[d.name]);
  else tip.html("file number " + d.name + "<br />" +file_descr[d.name]);
  tip.show();
}

 // hint? http://mbostock.github.io/d3/talk/20111018/tree.html

function mouseouted(d) {
  link
      .classed("link--target", false)
      .classed("link--source", false);

  node
      .classed("node--target", false)
      .classed("node--source", false);
 tip.hide();
}

  d3.select(self.frameElement).style("height", height + "px");
</script></div>
<a href="#" onClick="toggle_show('div-cc-clusters'); return false;">Show/Hide cluster list</a>
<div id="div-cc-clusters" style="display:none;">
<table class="cells">
<tr>%(header)s</tr>
%(data)s
</table>
</div>
""" % dict(treedata=treedata, merged_clusters=",".join(['"%s"'%x[0] for x in clusters]),
           file_descr=",".join(file_descr), cluster_descr=",".join(cluster_descr),
           header=header, data=data, method=method)
        self.write_html()
    # add_clutering_result()

    def add_cc_clustering_details(self, cc_data):
        """
        cc_data = [(i,j,cc,nref),...]
        """

        ccvalues = ",".join(["%.4f"%x[2] for x in [y for y in cc_data if y[2]==y[2]]])
        ccdata = ",".join(["{i:%d,j:%d,cc:%.4f,n:%d}"%(x[0]+1,x[1]+1,x[2],x[3]) for x in cc_data])
        ccdata = ccdata.replace("nan", "NaN")
        ncols = max([max(x[0],x[1]) for x in cc_data]) + 1

        # Histogram
        self._html_clustering_details = """
<div id="cchist"></div>
<script>
// Reference: http://bl.ocks.org/mbostock/3048450
var ccvalues = [%(ccvalues)s];
var margin = {top: 10, right: 30, bottom: 30, left: 30},
    width = 800 - margin.left - margin.right,
    height = 400 - margin.top - margin.bottom;

var x = d3.scale.linear()
    .domain([-1., 1.])
    .range([0, width]);

// Generate a histogram using twenty uniformly-spaced bins.
var data = d3.layout.histogram()
    .bins(x.ticks(30)).range([-1.,1.])
    (ccvalues);

var y = d3.scale.linear()
    .domain([0, d3.max(data, function(d) { return d.y; })])
    .range([height, 0]);
var xAxis = d3.svg.axis()
    .scale(x)
    .orient("bottom");

var svg = d3.select("#cchist").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
  .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

var bar = svg.selectAll(".bar")
    .data(data)
  .enter().append("g")
    .attr("class", "bar")
    .attr("transform", function(d) { return "translate(" + x(d.x) + "," + y(d.y) + ")"; });

bar.append("rect")
    .attr("x", 1)
    .attr("width", x(-1+data[0].dx) - 1)
    .attr("height", function(d) { return height - y(d.y); });

bar.append("text")
    .attr("dy", ".75em")
    .attr("y", 6)
    .attr("x", x(-1+data[0].dx) / 2)
    .attr("text-anchor", "middle")
    .text(function(d) { return d3.format(",.0f")(d.y); });
svg.append("g")
    .attr("class", "x axis")
    .attr("transform", "translate(0," + height + ")")
    .call(xAxis)
    .append('text')
     .text("CC")
     .attr('transform','translate('+(width+10)+',0)');
</script>
""" % dict(ccvalues=ccvalues)

    # Matrix heatmap
        self._html_clustering_details += """
<div id="ccheatmap"></div>
<script>
  // Reference: Day / Hour Heatmap http://bl.ocks.org/tjdecke/5558084
  //            Days-Hours Heatmap http://bl.ocks.org/oyyd/859fafc8122977a3afd6
  //            Method: Hierarchical clustering with SciPy and visualization in D3.js http://blog.nextgenetics.net/?e=44
  var ccdata = [%(ccdata)s];
  var ncols = %(ncols)d;
  var margin = { top: 50, right: 0, bottom: 100, left: 30 },
      width = 800 - margin.left - margin.right,
      height = 800 - margin.top - margin.bottom,
      gridSize = Math.floor(width / ncols);

  var svg = d3.select("#ccheatmap").append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  //define a color scale using the min and max expression values
  var colorScale = d3.scale.linear()
        .domain([-1, 0, 1])
        .range(["blue", "white", "red"]); // ['#f6faaa','#FEE08B','#FDAE61','#F46D43','#D53E4F','#9E0142']

  var cards = svg.selectAll(".i")
        .data(ccdata);

  // http://stackoverflow.com/questions/17776641/fill-rect-with-pattern
  svg
    .append('defs')
    .append('pattern')
      .attr('id', 'diagonalHatch')
      .attr('patternUnits', 'userSpaceOnUse')
      //.attr('width', gridSize)
      .attr('width', 4)
      //.attr('height', gridSize)
      .attr('height', 4)
    .append('path')
      //.attr('d', 'M0,'+gridSize+' l'+gridSize+',-'+gridSize)
      .attr('d', 'M0,4 l4,-4')
      .attr('stroke', '#000000')
      .attr('stroke-width', 1);

  cards.enter().append("rect")
      .attr('width',gridSize)
      .attr('height',gridSize)
      .attr('x', function(d) {
            return (d.i * gridSize);// + 25;
           })
      .attr('y', function(d) {
            return (d.j * gridSize);// + 50;
           })
      .style('fill',function(d) {
            return d.cc==d.cc ? colorScale(d.cc) : 'url(#diagonalHatch)';
           });

  var xAxis = d3.svg.axis()
        .orient('top')
        .ticks(10),
      yAxis = d3.svg.axis()
        .orient('left')
        .ticks(10);

  //render axises
  xAxis.scale(d3.scale.linear().range([0,width]).domain([0, ncols-1]));
  yAxis.scale(d3.scale.linear().range([0,width]).domain([0, ncols-1]));
  svg.append('g')
    .attr('class','x axis')
    .call(xAxis);

  svg.append('g')
    .attr('class','y axis')
    .call(yAxis);

  var tipmat = d3.tip()
       .attr('class', 'd3-tip')
       .offset([-10, 100]);
  svg.call(tipmat);

   cards
       .on('mouseover', function(d) {
            d3.select(this)
               .attr('stroke-width',1)
               .attr('stroke','black')

            tipmat.html(d.i + ' vs ' + d.j +' (nref: ' + d.n + ')<br />' + 'CC= ' + d.cc + '<br />');
            tipmat.show();
      })
      .on('mouseout', function(d,i,j) {
         d3.select(this)
            .attr('stroke-width',0)
            .attr('stroke','none')
         tipmat.hide();
      });
</script>

""" % dict(ccdata=ccdata, ncols=ncols)
    # add_cc_clustering_details()

    def _make_merge_table_framework(self):
        return """
<h2>Merging summary</h2>

<table class="merge">
 <tr>
  <th rowspan="2" colspan="2">cluster</th><th rowspan="2" title="Cluster height">ClH</th><th rowspan="2" title="Linear Cell Variation defined in BLEND">LCV</th><th rowspan="2" title="absolute Linear Cell Variation defined in BLEND">aLCV</th><th rowspan="2" title="Number of all datasets in the cluster">#DS<br />all</th><th rowspan="2" title="Number of actually merged datasets in the cluster">#DS<br />used</th>
  <th colspan="5">Overall</th>
  <th colspan="5">Outer shell</th>
  <th colspan="7">Inner shell</th>
  <th rowspan="2" title="ML estimate of isotropic Wilson B-factor by phenix.xtriage"><i>B</i><sub>Wilson</sub></th>
  <th colspan="2" title="Anisotropic resolution cutoffs (based on CC1/2=0.5) for best and worst directions">Aniso resol</th>
  <th rowspan="2" title="Resolution estimate based on CC1/2"><i>d</i><sub>min</sub><br>est.</th>
 </tr>
 <tr>
  <th>Cmpl</th><th>Redun</th><th>I/&sigma;(I)</th><th><i>R</i><sub>meas</sub></th><th>CC<sub>1/2</sub></th>
  <th>Cmpl</th><th>Redun</th><th>I/&sigma;(I)</th><th><i>R</i><sub>meas</sub></th><th>CC<sub>1/2</sub></th>
  <th>Cmpl</th><th>Redun</th><th>I/&sigma;(I)</th><th><i>R</i><sub>meas</sub></th><th>CC<sub>1/2</sub></th><th>SigAno</th><th>CC<sub>ano</sub></th>
  <th>best</th><th>wrst</th>
 </tr>
 %s
</table>
"""
    # _make_merge_table_framework()

    def _make_merge_plot_framework(self):
        axis_opts = "ClH   LCV aLCV ds.all ds.used  Cmpl Redun I/sigI Rmeas CC1/2 Cmpl.ou Red.ou I/sig.ou Rmeas.ou CC1/2.ou Cmpl.in Red.in I/sig.in Rmeas.in CC1/2.in SigAno.in CCano.in WilsonB aniso.best aniso.worst dmin.est".split()

        axis_opts_x, axis_opts_y = [], []
        for a in axis_opts:
            es = ' selected="selected"' if a=="Redun" else ""
            axis_opts_x.append('<option value="%s"%s>%s</option>'%(a,es,a))
            es = ' selected="selected"' if a=="CC1/2" else ""
            axis_opts_y.append('<option value="%s"%s>%s</option>'%(a,es,a))
            
        return """
<form name="merge_plot_selector">
<table>
<tr>
<td>
<div id="chartdiv_merge" style="width: 800px; height: 600px;"></div>
<script>
 var chart_merge = AmCharts.makeChart("chartdiv_merge", {
    "type": "xy",
    "theme": "none",  
    "legend": {
        "useGraphSettings": true,
    },
    "dataProvider": [%(data)s],
    "valueAxes": [{
        "position":"bottom",
        "axisAlpha": 0,
        "dashLength": 1,
        "title": "Redun"
    }, { 
        "axisAlpha": 0,
        "dashLength": 1,
        "position": "left",
        "title": "CC1/2"
    }],
    "graphs": [{
        "balloonText": "[[description]]: x=[[x]], y=[[y]]",
        "bullet": "round",
        "bulletSize": 8,
        "lineAlpha": 0,
        "xField": "Redun",
        "yField": "CC1/2",
        "descriptionField": "cls",
        "title": "Clusters",
        "hidden": false,
        "lineColor": "#FF6600",
        "fillAlphas": 0
    }],
    "chartScrollbar": {},
});
</script>
</td>
<td>X-axis:<br \>
<select name="xaxis" style="width: 150px" size="24" onChange="chart_merge.graphs[0].xField=this.value; chart_merge.valueAxes[0].title=this.value; chart_merge.validateNow(); chart_merge.validateData();">
%(axis_opts_x)s
</select>
</td>
<td>Y-axis:<br \>
<select name="xaxis" style="width: 150px" size="24" onChange="chart_merge.graphs[0].yField=this.value; chart_merge.valueAxes[1].title=this.value; chart_merge.validateNow(); chart_merge.validateData();">
%(axis_opts_y)s
</select>
</td>
</tr>
</table>
</form>
""" % dict(axis_opts_x="\n".join(axis_opts_x), axis_opts_y="\n".join(axis_opts_y), data="%s")
    # _make_merge_plot_framework()

    def add_merge_result(self, workdir, clh, LCV, aLCV, xds_files, num_files, stats):
        axis_opts = "cls ClH   LCV aLCV ds.all ds.used  Cmpl Redun I/sigI Rmeas CC1/2 Cmpl.ou Red.ou I/sig.ou Rmeas.ou CC1/2.ou Cmpl.in Red.in I/sig.in Rmeas.in CC1/2.in SigAno.in CCano.in WilsonB aniso.best aniso.worst dmin.est".split()

        cls = os.path.relpath(workdir, self.params.workdir)
        tmps = "%12s %5.2f %4.1f %4.1f %6d %7d %5.1f %5.1f %6.2f %5.1f %5.1f %7.1f %6.1f % 8.2f % 8.1f %8.1f %7.1f %6.1f % 8.2f % 8.1f %8.1f %9.1f %8.1f %7.2f %.2f %.2f %.2f"
        tmps = tmps % (cls, clh, LCV, aLCV,
                       len(xds_files), num_files,
                       stats["cmpl"][0],
                       stats["redundancy"][0],
                       stats["i_over_sigma"][0],
                       stats["r_meas"][0],
                       stats["cc_half"][0],
                       stats["cmpl"][2],
                       stats["redundancy"][2],
                       stats["i_over_sigma"][2],
                       stats["r_meas"][2],
                       stats["cc_half"][2],
                       stats["cmpl"][1],
                       stats["redundancy"][1],
                       stats["i_over_sigma"][1],
                       stats["r_meas"][1],
                       stats["cc_half"][1],
                       stats["sig_ano"][1],
                       stats["cc_ano"][1],
                       stats["xtriage_log"].wilson_b,
                       #stats["xtriage_log"].anisotropy,
                       stats["aniso"]["d_min_best"],
                       stats["aniso"]["d_min_worst"],
                       stats["dmin_est"],
                       )

        tmptmp = tmps.replace("nan",'"nan"').split()
        tmptmp[0] = '"%s"' % tmptmp[0]
        self.html_merge_plot_data.append("{%s}"%",".join(['"%s":%s'%tuple(x) for x in zip(axis_opts, tmptmp)]))

        tmps = "".join(["<td>%s</td>"%x for x in tmps.split()])
        idno = len(self.html_merge_results)
        if self.params.program == "xscale":
            table_snip = xscalelp.snip_symm_and_cell(stats["lp"]) + "\n"
            table_snip += xscalelp.snip_stats_table(stats["lp"])
            if stats["aniso"]:
                table_snip += "\nAnisotropy:\n"
                if stats["aniso"]["has_anisotropy"]:
                    if stats["aniso"]["aniso_cutoffs"]:
                        lab_maxlen = max(len("direction"), max([len(x[1]) for x in stats["aniso"]["aniso_cutoffs"]]))
                        table_snip += ("%"+str(lab_maxlen)+"s B_eigen Resol(CC1/2=0.5)\n")%"direction" # XXX if not 0.5?
                        for _, lab, reso, eig in stats["aniso"]["aniso_cutoffs"]:
                            table_snip += ("%"+str(lab_maxlen)+"s %7.2f %.2f\n") % (lab, eig, reso)
                    else:
                        table_snip += " Anisotropy analysis failed. Check the logfile.\n"
                else:
                    table_snip += " No anisotropy in this symmetry.\n"
                    
        else:
            table_snip = ""
        tmps2 = """ <tr><td onClick="toggle_show2(this, 'merge-td-%d');" id="merge-td-mark-%d"">&#x25bc;</td>%s</tr>\n""" %(idno,idno,tmps)
        tmps2 += """ <tr><td style="padding: 0px;"><td colspan="27" style="display:none;padding:0px;" id="merge-td-%d"><pre>%s</pre></td></tr>""" % (idno, table_snip)

        self.html_merge_results.append(tmps2)
        
    # add_merge_result()

    def write_html(self):
        ofs = open(os.path.join(self.root, "report.html"), "w")
        ofs.write(self.html_head)
        ofs.write(self.html_params)
        ofs.write(self.html_inputfiles)
        ofs.write(self.html_clustering)

        # merging table
        ofs.write(self._make_merge_table_framework()%"".join(self.html_merge_results))
        tmps = ";".join(["show_or_hide('merge-td-mark-%d', 'merge-td-%d', 0)"%(x,x) for x in range(len(self.html_merge_results))])
        ofs.write("""<a href="#" onClick="%s;return false;">Expand all</a> / """ % tmps)
        tmps = ";".join(["show_or_hide('merge-td-mark-%d', 'merge-td-%d', 1)"%(x,x) for x in range(len(self.html_merge_results))])
        ofs.write("""<a href="#" onClick="%s;return false;">Collapse all</a>\n""" % tmps)

        # merging plot
        ofs.write(self._make_merge_plot_framework()%",".join(self.html_merge_plot_data))

        ofs.write("\n</body></html>")
        ofs.close()
    # write_html()

# class HtmlReportMulti
