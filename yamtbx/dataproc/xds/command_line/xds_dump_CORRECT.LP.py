#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals

from yamtbx.dataproc.xds import correctlp

def run(f):
    lp = correctlp.CorrectLp(f)
    print(lp.snippets["ISa"])
    print(lp.snippets["table1"]) 
    return lp
# run()

if __name__ =="__main__":
    import sys
    for f in sys.argv[1:]:
        print(f)
        print("====================")
        run(f)
        print()
