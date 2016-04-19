# LIBTBX_SET_DISPATCHER_NAME yamtbx.eiger_splitdata_after
"""
(c) RIKEN 2016. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

from yamtbx.dataproc.command_line import eiger_splitdata_after

if __name__ == "__main__":
    import sys
    eiger_splitdata_after.run(sys.argv[1], int(sys.argv[2]))
