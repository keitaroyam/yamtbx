# LIBTBX_SET_DISPATCHER_NAME kamo.multi_merge
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

from yamtbx.dataproc.auto.command_line import multi_merge

if __name__ == "__main__":
    import sys
    multi_merge.run_from_args(sys.argv[1:])
