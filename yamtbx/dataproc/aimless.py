"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals
from yamtbx import util
from yamtbx.util import batchjob
import os
import re
import shutil
from libtbx import easy_mp
import pickle as pickle

aimless_comm = "aimless"
def run_aimless(mtzin, wdir, anomalous=False, d_min=None, prefix=None, add_stdin=None):
    if prefix is None:
        prefix = os.path.splitext(os.path.basename(mtzin))[0] + "_aimless"

    args = """\
HKLIN %(mtzin)s \
HKLOUT %(prefix)s.mtz \
SCALES %(prefix)s.scala \
ROGUES %(prefix)s_rogues.log \
NORMPLOT %(prefix)s_normplot.xmgr \
ANOMPLOT %(prefix)s_anomplot.xmgr \
PLOT %(prefix)s_surface_plot.plt \
CORRELPLOT %(prefix)s_correlplot.xmgr \
ROGUEPLOT %(prefix)s_rogueplot.xmgr \
""" % dict(mtzin=mtzin, prefix=prefix)

    stdin = ""
    if d_min is not None:
        stdin += "resolution high %.3f\n" % d_min
    if anomalous:
        stdin += "anomalous on\n"
    if add_stdin is not None:
        stdin += add_stdin

    exitcode, output, err = util.call(aimless_comm, args,
                                      stdin=stdin, wdir=wdir)

    open(os.path.join(wdir, prefix+".log"), "w").write(output)

# run_aimless

def snip_summary(login):
    ret = ""
    read_flag = False
    for l in open(login):
        if "Summary data for" in l:
            read_flag = True
        elif read_flag:
            if l.startswith("$$ <!--SUMMARY_END-->"): break
            ret += l
    return ret
# snip_summary()

def read_summary(login):
    # Total, InnerShell, OuterShell
    read_3 = lambda l: (l[38:38+10], l[38+10:38+20], l[38+20:38+30])
    re_read3_1 = re.compile("([-0-9]+\.[0-9]) *([-0-9]+\.[0-9]) *([-0-9]+\.[0-9])")
    re_read3_2 = re.compile("([-0-9]+\.[0-9]{2}) *([-0-9]+\.[0-9]{2}) *([-0-9]+\.[0-9]{2})")
    re_read3_3 = re.compile("([-0-9]+\.[0-9]{3}) *([-0-9]+\.[0-9]{3}) *([-0-9]+\.[0-9]{3})")
    read_3_1 = lambda l: re_read3_1.search(l[38:]).groups()
    read_3_2 = lambda l: re_read3_2.search(l[38:]).groups()
    read_3_3 = lambda l: re_read3_3.search(l[38:]).groups()

    ret = {}
    for l in snip_summary(login).splitlines():
        if l.startswith("Low resolution limit"):
            ret["lowres"] = list(map(float, read_3_2(l)))
        elif l.startswith("High resolution limit"):
            ret["highres"] = list(map(float, read_3_2(l)))
        elif l.startswith("Rmerge  (all I+ and I-)"):
            ret["r_merge"] = list(map(float, read_3_3(l)))
        elif l.startswith("Rmeas (all I+ & I-)"):
            ret["r_meas"] = list(map(float, read_3_3(l)))
        elif l.startswith("Rpim (all I+ & I-)"):
            ret["r_pim"] = list(map(float, read_3_3(l)))
        elif l.startswith("Total number unique"):
            ret["nuniq"] = list(map(float, read_3(l)))
        elif l.startswith("Mean((I)/sd(I))"):
            ret["i_over_sigma"] = list(map(float, read_3_1(l)))
        elif l.startswith("Mn(I) half-set correlation CC(1/2)"):
            ret["cc_half"] = list(map(float, read_3_3(l)))
        elif l.startswith("Completeness"):
            ret["cmpl"] = list(map(float, read_3_1(l)))
        elif l.startswith("Multiplicity"):
            ret["redundancy"] = list(map(float, read_3_1(l)))
        elif l.startswith("DelAnom correlation between half-sets"):
            ret["cc_ano"] = list(map(float, read_3_3(l)))

    return ret
# read_summary()



def _calc_cchalf_by_removing_worker_1(tmpdir, mtzin, batch_info, inpfiles, anomalous_flag, d_min, iex, nproc=None):
    print("Doing", iex)
    if not os.path.exists(tmpdir): os.mkdir(tmpdir)

    inp_str = ""
    files = inpfiles[:iex] + inpfiles[iex+1:]
    for i, f in enumerate(files):
        brange = batch_info[f]
        inp_str += "RUN %3d BATCH %4d to %4d\n" % (i+1, brange[0], brange[1])

    if nproc is not None:
        inp_str += "PARALLEL %d\n" % nproc
    
    run_aimless(mtzin=os.path.relpath(mtzin, tmpdir),
                wdir=tmpdir,
                anomalous=anomalous_flag, d_min=d_min, prefix="aimless",
                add_stdin=inp_str)
# _calc_cchalf_by_removing_worker_1()

def _calc_cchalf_by_removing_worker_2(wdir, tmpdir, stat_bin, iex):
    assert stat_bin in ("total", "outer")
    print("Doing", iex)

    aimless_log = os.path.join(tmpdir, "aimless.log")

    table = read_summary(aimless_log)
    i_stat = 0 if stat_bin == "total" else 2
    cchalf_exi = table["cc_half"][i_stat]
    nuniq = table["nuniq"][i_stat]

    # backup .INP and .LP, and then remove directory.
    os.rename(aimless_log, os.path.join(wdir, "aimless.log.ex%.3d"%iex))
    shutil.rmtree(tmpdir)

    return iex, cchalf_exi, nuniq
# _calc_cchalf_by_removing_worker_2()


def calc_cchalf_by_removing(wdir, mtzin, batch_info, inpfiles, anomalous_flag, d_min, 
                            with_sigma=False, stat_bin="total",
                            nproc=1, nproc_each=None, batchjobs=None):
    assert not with_sigma # Not supported now
    assert stat_bin in ("total", "outer")

    if not os.path.exists(wdir): os.makedirs(wdir)

    datout = open(os.path.join(wdir, "cchalf.dat"), "w")
    datout.write("idx exfile cc1/2(%s) Nuniq\n" % stat_bin)

    cchalf_list = [] # (i_ex, CC1/2, Nuniq)

    tmpdirs = [os.path.join(wdir, "work.%.3d"%iex) for iex in range(len(inpfiles))]

    # Run Aimless 
    if batchjobs is not None:
        pkltmp = os.path.join(wdir, "tmp.pkl")
        pickle.dump(dict(mtzin=mtzin, batch_info=batch_info, inpfiles=inpfiles, 
                         anomalous_flag=anomalous_flag, d_min=d_min, nproc=nproc_each),
                    open(pkltmp, "wb"), 2)
        jobs = []
        for i, tmpdir in enumerate(tmpdirs):
            if not os.path.exists(tmpdir): os.mkdir(tmpdir)

            job = batchjob.Job(tmpdir, "aimless.sh", nproc=nproc_each)
            jobstr = """\
cd %s
yamtbx.dev.python - << +
import pickle
from yamtbx.dataproc.aimless import _calc_cchalf_by_removing_worker_1
params = pickle.load(open("%s/tmp.pkl", "rb"))
params["tmpdir"] = "%s"
params["iex"] = %d
_calc_cchalf_by_removing_worker_1(**params)
+
""" % (os.getcwd(), wdir, tmpdir, i)
            job.write_script(jobstr)
            batchjobs.submit(job)
            jobs.append(job)

        batchjobs.wait_all(jobs)
        os.remove(pkltmp)
    else:
        easy_mp.pool_map(fixed_func=lambda x: _calc_cchalf_by_removing_worker_1(x[1], mtzin, batch_info, inpfiles, anomalous_flag, d_min, x[0], nproc_each),
                         args=enumerate(tmpdirs),
                         processes=nproc)

    # Finish runs
    cchalf_list = [_calc_cchalf_by_removing_worker_2(wdir, x[1], stat_bin, x[0]) for x in enumerate(tmpdirs)]

    for iex, cchalf_exi, nuniq in cchalf_list:
        datout.write("%3d %s %.4f %d\n" % (iex, inpfiles[iex], cchalf_exi, nuniq))


    cchalf_list.sort(key=lambda x: -x[1])
    print()
    print("# Sorted table")
    for idx, cch, nuniq in cchalf_list:
        print("%3d %-.4f %4d %s" % (idx, cch, nuniq, inpfiles[idx]))

    return cchalf_list
# calc_delta_cchalf()
