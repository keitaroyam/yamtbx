# LIBTBX_SET_DISPATCHER_NAME kamo.decide_resolution_cutoff
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

from yamtbx.dataproc.auto.command_line import decide_resolution_cutoff

if __name__ == "__main__":
    import sys
    decide_resolution_cutoff.run_from_args(sys.argv[1:])

