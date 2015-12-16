"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from yamtbx.dataproc import XIO
from yamtbx.dataproc.dataset import template_to_filenames, find_data_sets

def timestamps(img_files):
    img_date = []

    for f in img_files:
        im = XIO.Image(f)
        img_date.append((f, im.header["DateSeconds"]))

    for i in xrange(len(img_files)):
        print img_files[i],
        if i == 0:
            print 0
        else:
            print img_date[i][1] - img_date[i-1][1]


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        wdir = sys.argv[1]
    else:
        wdir = os.getcwd()

    for img_template, min_frame, max_frame in find_data_sets(wdir):
        img_files = template_to_filenames(img_template, min_frame, max_frame)
        timestamps(img_files)
