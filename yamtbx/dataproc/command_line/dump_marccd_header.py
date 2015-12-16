#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

import sys
from yamtbx.dataproc import marccd

if __name__ == "__main__":
    imgin = sys.argv[1]
    mccd = marccd.MarCCD(imgin)
    mccd.read_header()
    
