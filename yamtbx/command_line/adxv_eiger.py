# LIBTBX_SET_DISPATCHER_NAME yamtbx.adxv_eiger
"""
(c) RIKEN 2016. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

from yamtbx.dataproc.command_line import adxv_eiger

if __name__ == "__main__":
    import sys
    adxv_eiger.run(sys.argv[1:])
