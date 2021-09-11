#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals

from yamtbx.dataproc.xds import idxreflp
from yamtbx.dataproc.xds import get_xdsinp_keyword
from yamtbx.dataproc import dataset
from yamtbx.dataproc.adxv import Adxv

import iotbx.phil

import time
import os
import socket
import subprocess

master_params_str = """
adxv_bin = adxv
 .type = str
image = 1
 .type = int
xds_dir = "."
 .type = path
interactive = True
 .type = bool
"""

def run(params):
    xds_inp = os.path.join(params.xds_dir, "XDS.INP")
    spot_xds = os.path.join(params.xds_dir, "SPOT.XDS")

    spots = idxreflp.SpotXds(spot_xds).indexed_and_unindexed_by_frame_on_detector()
    inputs = get_xdsinp_keyword(xds_inp)

    filename_template = dict(inputs).get("NAME_TEMPLATE_OF_DATA_FRAMES", "")
    if filename_template == "":
        print("Error! Can't find filename from XDS.INP")
        return

    # Start adxv
    adxv = Adxv(params.adxv_bin)
    adxv.start(params.xds_dir)
    type_indexed = adxv.define_spot("blue")
    type_unindexed = adxv.define_spot("red")
    
    num = params.image

    while True:
        print("Showing image %d" % num)

        img_file = dataset.template_to_filenames(filename_template, num, num)[0]
        adxv.open_image(img_file)

        uninds = [[x[0], x[1], type_unindexed] for x in spots["unindexed"].get(num, [])]
        inds = [[x[0], x[1], type_indexed] for x in spots["indexed"].get(num, [])]

        print("Showing %d Indexed spots (blue)" % len(inds))
        print("Showing %d Unindexed spots (red)" % len(uninds))

        adxv.load_spots(inds+uninds)

        time.sleep(1) # wait until adxv finishes process..
        num = int(input("Next image number?: "))

# run()

if __name__ == "__main__":
    import sys

    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args
    
    run(params)
