# LIBTBX_SET_DISPATCHER_NAME kamo
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

from yamtbx.dataproc.auto.command_line import auto_data_proc_gui

if __name__ == "__main__":
    import sys
    auto_data_proc_gui.run_from_args(sys.argv[1:])
