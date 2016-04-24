# LIBTBX_SET_DISPATCHER_NAME yamtbx.eiger_binning
"""
(c) RIKEN 2016. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

from yamtbx.dataproc.command_line import eiger_binning

if __name__ == "__main__":
    import sys
    eiger_binning.run_from_args(sys.argv[1:])
