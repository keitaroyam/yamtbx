# LIBTBX_SET_DISPATCHER_NAME yamtbx.beam_direction_plot
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

from yamtbx.dataproc.command_line import beam_direction_plot

if __name__ == "__main__":
    import sys
    beam_direction_plot.run_from_args(sys.argv[1:])
