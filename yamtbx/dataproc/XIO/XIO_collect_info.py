#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __init__.py
# Maintained by P.Legrand
# 14th April 2005
# 

from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals
import os, sys
from . import  __init__

def _print_autopar(parDict):
    K = list(parDict.keys())
    K.sort()
    for k in K:
        print("%s = %s" % (k, parDict[k]))

def _test1(filename):
    line_fmt = "%(Wavelength)8.4f\t%(Distance).2f\t%(Nimages)4d\t%(PhiWidth).2f\t%(ImageNameTemplate)s"
    try:
        
        collect = __init__.Collect(filename)
        collect.interpretImage()
        CollectData = collect.image.info()
        CollectData["ImageRange"] = collect.lookup_imageRanges(forceCheck=False)
        CollectData["Nimages"] = CollectData["ImageRange"][-1][-1] - \
                                 CollectData["ImageRange"][0][0] + 1
        CollectData["ImageNameTemplate"] = os.path.split(collect.xdsTemplate)[-1]
        print(line_fmt % CollectData)
    except __init__.XIOError as _mess:
        print(_mess)
        print("\nError: Can't access to file(s) %s.\nStop." % inputf)
        sys.exit(2)
        #dc = __init__.Collect(filename)
    #im.info()
    #im.info2()
    #im.data_info()

def _export(filename, format):
    try:
        datacoll = __init__.Collect(filename)
        datacoll.interpretImage()
        datacoll.lookup_imageRanges()
    except __init__.XIOError:
        print("\nError while trying to acceess %s.\nSorry." % filename)
        sys.exit()
    newPar = datacoll.export(format)
    if format == "adp":
        newPar["anomalous"] = "off"
        newPar["spg"] = 0
    print(datacoll.export_template(format))

if __name__ == "__main__":

    if len(sys.argv) > 1:
        if "-xds" in sys.argv:
            sys.argv.remove("-xds")
            filename = sys.argv[1]
            _export(filename, "xds")
        elif "-mos" in sys.argv:
            sys.argv.remove("-mos")
            filename = sys.argv[1]
            _export(filename, "mosflm")
        elif "-adp" in sys.argv:
            sys.argv.remove("-adp")
            filename = sys.argv[1]
            _export(filename, "adp")
        else:
            filename = sys.argv[1]
            _test1(filename)
