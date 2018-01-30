# LIBTBX_SET_DISPATCHER_NAME yamtbx.easy_anode
"""
(c) RIKEN 2018. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

from yamtbx.dataproc.command_line import easy_anode

if __name__ == "__main__":
    import sys

    if len(sys.argv) != 3:
        print "Usage: yamtbx.easy_anode hkl_in pdb_in"
        quit()

    hklin = sys.argv[1]
    pdbin = sys.argv[2]
    easy_anode.run(hklin, pdbin)
