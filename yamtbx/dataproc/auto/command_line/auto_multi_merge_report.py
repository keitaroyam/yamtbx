"""
(c) RIKEN 2017. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals
import os
import glob
import time
import re
from yamtbx.dataproc.xds import get_xdsinp_keyword
from yamtbx.dataproc.auto.command_line import auto_multi_merge
from yamtbx.dataproc.xds import xscalelp
from yamtbx.util import replace_forbidden_chars
from yamtbx.dataproc.command_line import beam_direction_plot
from libtbx.utils import null_out

def total_deg_from_xds_inp(xdsinp):
    kwds = dict(get_xdsinp_keyword(xdsinp))
    dr = list(map(int, kwds["DATA_RANGE"].split()))
    osc = float(kwds["OSCILLATION_RANGE"])
    return (dr[1]-dr[0]+1)*osc
# total_deg_from_xds_inp()

def run(csvin, prefix, rootdir, datadir=None):
    html_head = """\
<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<script>
  var toggle_show = function(caller, obj_id) {
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
</script>

<style>
pre {
    font-family: Consolas, 'Courier New', Courier, Monaco, monospace;
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

</style>

</head>

<body>
<h1>KAMO.AUTO_MULTI_MERGE report</h1>
<div align="right">
workdir: %(wd)s<br />
created on %(cdate)s
</div>
<hr>

<table class="merge">
<tr>
 <th colspan="2">Sample</th>
 <th>#Collected</th>
 <th>#Processed</th>
 <th>#Mergeable</th>
 <th>Symmetry</th>
 <th>Unit cell</th>
 <th>Cmpl(all)</th>
 <th>Mult(all)</th>
 <th><i>d</i><sub>min</sub></th>
</tr>
""" % dict(wd=rootdir, cdate=time.strftime("%Y-%m-%d %H:%M:%S"))

    ofs = open(os.path.join(rootdir, "report.html"), "w")
    ofs.write(html_head)

    samples = auto_multi_merge.read_sample_info(csvin, datadir)
    
    for name in samples:
        workdir_rel = "%s%s" % (prefix,
                                replace_forbidden_chars(name).replace(" ","_"))
        workdir = os.path.join(rootdir, workdir_rel)
        print(workdir)
        assert os.path.isdir(workdir)
        topdirs = samples[name][0]
        ncol,dcol,npro,dpro,nmrg,dmrg = 0,0,0, 0,0,0
        symm, cell = "?", "?"
        rephtml = "#"
        cmpl, redun, dmin = 0, 0, float("nan")
        symm_group_info = "n/a" # now useless
        deg_dict = {}
        for topdir in topdirs:
            for root, dirnames, filenames in os.walk(topdir):
                if "XDS.INP" not in filenames: continue
                print("Checking", name, root)
                deg = total_deg_from_xds_inp(os.path.join(root, "XDS.INP"))
                deg_dict[root] = deg
                ncol += 1
                dcol += deg
                if "XDS_ASCII.HKL_noscale" in filenames:
                    npro += 1
                    dpro += deg

        mrg_lst = os.path.join(workdir, "formerge.lst")
        if os.path.isfile(mrg_lst):
            mrg_dirs = [os.path.dirname(x.strip()) for x in open(mrg_lst).readlines()]
            nmrg = len(mrg_dirs)
            dmrg = sum([deg_dict.get(x, 0) for x in mrg_dirs])

            beam_plot_png = os.path.join(workdir, "beam_plot.png")
            #if not os.path.isfile(beam_plot_png):
            #    beam_direction_plot.run_from_args([mrg_lst, 'plot_out="%s"'%beam_plot_png])

        mrg_log = os.path.join(workdir, "multi_merge.log")
        if os.path.isfile(mrg_log):
            flag_read = False
            flag_first = True
            symm_group_info = ""
            for l in open(mrg_log):
                if flag_first and "members:" in l:
                    flag_read = True
                    flag_first = False
                elif flag_read and ("members:" in l or l.strip()==""):
                    flag_read = False
                elif flag_read and "Members=" not in l:
                    symm_group_info += l
                    
                if "group_choice=" in l:
                    symm, cell = re.search("symmetry= ([^\(]+) \((.+)\)", l).groups()

        mrg_dirs = glob.glob(os.path.join(workdir, "*final/"))
        best_result_loc, best_table_snip = "N/A", ""
        if mrg_dirs:
            mrg_dir = mrg_dirs[0]
            dmin = float(re.search("_([0-9\.]+)A_final/", mrg_dir).group(1))
            cls_dat = glob.glob(os.path.join(mrg_dir, "*", "*_cluster_summary.dat"))[0]
            tmp = open(cls_dat).readlines()[3].split()
            cmpl, redun = float(tmp[3]), float(tmp[4])
            rephtml = os.path.relpath(os.path.join(mrg_dir, "report.html"), rootdir)
            summary_dat = os.path.join(mrg_dir, "cluster_summary.dat")
            best_result = auto_multi_merge.choose_best_result(summary_dat, null_out())
            if best_result:
                best_result_loc = os.path.dirname(best_result)
                lp = os.path.join(best_result_loc, "XSCALE.LP")
                best_table_snip = xscalelp.snip_symm_and_cell(lp) + "\n" + xscalelp.snip_stats_table(lp)
                best_result_loc = os.path.relpath(best_result_loc, rootdir) # to show in html

            
        html_tr = """\
<tr>
 <td onClick="toggle_show(this, 'sample-td-%(name)s');" id="sample-td-mark-%(name)s"">&#x25bc;</td>
 <td><a href="%(rephtml)s">%(name)s</a></td>
 <td>%(ncol)d (%(dcol).0f&deg;)</td>
 <td>%(npro)d (%(dpro).0f&deg;)</td>
 <td>%(nmrg)d (%(dmrg).0f&deg;)</td>
 <td>%(symm)s</td>
 <td>%(cell)s</td>
 <td>%(cmpl).1f</td>
 <td>%(redun).1f</td>
 <td>%(dmin).1f</td>
</tr>
<tr>
 <td style="padding: 0px;"></td>
  <td colspan="9" style="display:none;padding:0px;" id="sample-td-%(name)s">
  BEST RESULT: <a href="%(best_result_loc)s">%(best_result_loc)s</a> <br><br>
<pre>
 %(best_table_snip)s
</pre>
  
 </td>
</tr>

""" % locals()
        ofs.write(html_tr)
        ofs.flush()
        #break

    ofs.write("\n</table>\n")
    ofs.write("\n</body></html>\n")
    ofs.close()

    print("Done!")
    print("firefox", os.path.join(rootdir, "report.html"))
# run()

def run_from_args(argv):
    if len(argv) == 2:
        csvin, prefix = argv
        datadir = None
    elif len(argv) == 3:
        csvin, prefix,datadir = argv

    run(csvin, prefix, os.getcwd(), datadir)

# run_from_args()

if __name__ == "__main__":
    import sys
    run_from_args(sys.argv[1:])
