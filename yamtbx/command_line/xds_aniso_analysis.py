# LIBTBX_SET_DISPATCHER_NAME yamtbx.xds_aniso_analysis
"""
(c) RIKEN 2017. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

from yamtbx.dataproc.xds.command_line import xds_aniso_analysis

if __name__ == "__main__":
    import sys
    xds_aniso_analysis.run_from_args(sys.argv[1:])

