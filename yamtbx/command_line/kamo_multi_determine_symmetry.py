# LIBTBX_SET_DISPATCHER_NAME kamo.multi_determine_symmetry
"""
(c) RIKEN 2017. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

from yamtbx.dataproc.auto.command_line import multi_determine_symmetry

if __name__ == "__main__":
    import sys
    multi_determine_symmetry.run_from_args(sys.argv[1:])
