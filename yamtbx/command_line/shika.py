# LIBTBX_SET_DISPATCHER_NAME yamtbx.shika
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

from yamtbx.dataproc.command_line import spot_finder_gui

if __name__ == "__main__":
    import sys
    spot_finder_gui.run_from_args(sys.argv[1:])
