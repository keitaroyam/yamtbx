# LIBTBX_SET_DISPATCHER_NAME yamtbx.eiger_bslz4_after
"""
(c) RIKEN 2016. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

import dxtbx.format # to set HDF5_PLUGIN_PATH in phenix environment
from yamtbx.dataproc.command_line import eiger_bslz4_after

if __name__ == "__main__":
    import sys
    eiger_bslz4_after.run_from_args(sys.argv[1:])
