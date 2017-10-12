# LIBTBX_SET_DISPATCHER_NAME kamo.auto_multi_merge
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

from yamtbx.dataproc.auto.command_line import auto_multi_merge

if __name__ == "__main__":
    import sys
    auto_multi_merge.run_from_args(sys.argv[1:])
