"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from yamtbx.dataproc.auto.multi_merging.resolve_reindex import ReferenceBased

if __name__ == "__main__":
    import sys
    lst,ref = sys.argv[1:3]
    xac_files = map(lambda x:x.strip(), open(lst))

    rb = ReferenceBased(xac_files, ref, log_out=sys.stdout)
    rb.assign_operators()
    rb.debug_write_mtz()
    rb.modify_xds_ascii_files()

        
