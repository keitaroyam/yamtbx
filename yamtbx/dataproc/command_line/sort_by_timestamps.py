#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

"""
Usage:
  sort_by_timestamps.py scan1_*img

And you will get:

DateStr                    TimeDiff Filename
Mon Jan 21 04:47:24 2013          0 scan1_000001.img
Mon Jan 21 04:47:27 2013          3 scan1_000002.img
Mon Jan 21 04:47:30 2013          6 scan1_000003.img
....
"""

from yamtbx.dataproc import XIO
import time
import os

if __name__ == "__main__":
    import sys
    import optparse

    parser = optparse.OptionParser(usage="usage: %prog [options] img_files")
    parser.add_option("--file-timestamp", "-f", action="store_true", dest="file_timestamp", 
                      help="Use file system time stamp instead of detector header")

    opts, args = parser.parse_args(sys.argv)

    img_date = []

    for f in args[1:]:
        if opts.file_timestamp:
            img_date.append((f, os.path.getmtime(f)))
        else:
            im = XIO.Image(f)
            img_date.append((f, im.header["DateSeconds"]))
        
    img_date.sort(key=lambda x:x[1])

    print "DateStr                    TimeDiff Filename"

    for i in xrange(len(img_date)):
        print time.ctime(img_date[i][1]),
        print "%10d %s" % (img_date[i][1] - img_date[0][1], 
                           img_date[i][0])

            
            
