# LIBTBX_SET_DISPATCHER_NAME yamtbx.run_xscale
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

from yamtbx.dataproc.xds.command_line import run_xscale

if __name__ == "__main__":
    import sys
    run_xscale.run_from_args(sys.argv[1:])

