"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals
from yamtbx.dataproc.pointless import parse_pointless_output_for_runs
from yamtbx.dataproc.pointless import parse_pointless_output_for_input_files
import collections

def run(log):
    logstr = open(log).read()
    runinfo = parse_pointless_output_for_runs(logstr)
    xds_files = parse_pointless_output_for_input_files(logstr)

    assert len(runinfo) == len(xds_files)

    batch_info = collections.OrderedDict([(x[0], (x[1][1:3])) for x in zip(xds_files, runinfo)])

    print("file run batch.s batch.e")
    maxlen = max([len(x) for x in batch_info])
    for i, (xdsfile, (sb, eb)) in enumerate(batch_info.items()):
        print(("%-"+str(maxlen)+"s %4d %8d %8d") % (xdsfile, i+1, sb, eb))
        
    return batch_info
# run()

if __name__ == "__main__":
    import sys
    run(sys.argv[1]) # give pointless log
