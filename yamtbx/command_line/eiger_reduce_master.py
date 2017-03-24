# LIBTBX_SET_DISPATCHER_NAME yamtbx.eiger_reduce_master
"""
(c) RIKEN 2017. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

from yamtbx.dataproc.command_line import eiger_reduce_master

if __name__ == "__main__":
    import sys
    eiger_reduce_master.run_from_args(sys.argv[1:], command_name="yamtbx.eiger_reduce_master")
