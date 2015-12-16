#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""


import os
import sys
from yamtbx.dataproc.bl_logfiles import BssDiffscanLog
import time
import math

def remove_overwritten(scans):
    del_idx = []
    for i, scan in enumerate(scans):
        for scan2 in scans[i+1:]:
            if scan.filename_template == scan2.filename_template:
                if i not in del_idx:
                    del_idx.append(i)
                    #print "Delete it!", scan.filename_template, i

    for i in reversed(del_idx):
        del scans[i]
    return scans

def get_mean_time(files):
    if len(files) < 2:
        return float("nan"), float("nan"), float("nan")
    img_date = []

    for f in files:
        if os.path.isfile(f):
            img_date.append((f, os.path.getmtime(f)))
        elif os.path.isfile(f+".bz2"):
            img_date.append((f, os.path.getmtime(f+".bz2")))
    img_date.sort(key=lambda x:x[1])

    if len(img_date) < 2:
        return float("nan"), float("nan"), float("nan")
    length = img_date[-1][1] - img_date[0][1]
    mean = length / (len(img_date)-1)
    if len(img_date)-2 > 0:
        sd = math.sqrt(sum([(img_date[i][1] - img_date[i-1][1] - mean)**2 for i in xrange(1,len(img_date))])/(len(img_date)-2))
    else:
        sd = 0
    return length, mean, sd

if __name__ == "__main__":
    if len(sys.argv) > 1:
        parentdir = sys.argv[1]
    else:
        parentdir = os.getcwd()

    print "date file shutter exp framerate sec.total sec.mean sec.sd vpoints hpoints lines num"
    for root, dirnames, filenames in os.walk(parentdir):
        if "diffscan.log" in filenames:
            scanlog = BssDiffscanLog(os.path.join(root, "diffscan.log"))
            print "#", os.path.relpath(root, parentdir)
            for scan in remove_overwritten(scanlog.scans):
                files = [os.path.join(root,f) for f, c in scan.filename_coords]
                #etime = time.mktime(scan.date.timetuple())
                print '"%s"' % scan.date.strftime("%Y/%m/%d %H:%M:%S"),
                print scan.filename_template,
                print "shutterless" if scan.is_shutterless() else "shutter",
                print scan.exp_time, scan.frame_rate,
                print "%4.2f %4.2f %.3f" % get_mean_time(files),
                print "%4d %4d %4d %4d" %(scan.vpoints, scan.hpoints, min(scan.vpoints, scan.hpoints), len(files))
            print "#"
