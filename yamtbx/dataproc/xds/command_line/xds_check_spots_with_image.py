#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

from yamtbx.dataproc.xds import idxreflp
from yamtbx.dataproc.xds import get_xdsinp_keyword
from yamtbx.dataproc import dataset

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

class Adxv:
    def __init__(self, adxv_bin):
        self.adxv_comm = adxv_bin + " -socket %d"
        self.sock = None
        self.adxv_proc = None
        self.spot_type_counter = -1
    # __init__()

    def start(self, wdir):
        # find available port number
        sock_test = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock_test.bind(("localhost", 0))
        adxv_port = sock_test.getsockname()[1]
        sock_test.close()

        # start adxv
        self.adxv_proc = subprocess.Popen(self.adxv_comm % adxv_port, shell=True,
                                          cwd=wdir)

        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

        for i in xrange(4): # try for 2 seconds.
            try:
                self.sock.connect(("localhost", adxv_port))
                break
            except socket.error:
                time.sleep(.5)
                continue
    # start()

    def define_spot(self, color, radius=0):
        self.spot_type_counter += 1
        sent = self.sock.send("define_type %d color %s radius %d\n"%(self.spot_type_counter, color, radius))
        print sent
        if sent == 0:
            print "define_spot failed!"

        sent = self.sock.send("box 20 20\n")
        return self.spot_type_counter
    # define_spot()

    def open_file(self, imgfile):
        sent = self.sock.send("load_image %s\n"%imgfile)
        if sent == 0:
            print "adxv load failed!"

        #sock.close()
    # open_file()

    def load_spots(self, spots):
        if len(spots) == 0:
            return

        sent = self.sock.send("load_spots %d\n" % len(spots))

        for x, y, t in spots:
            sent = self.sock.send("%.2f %.2f %d\n" % (x, y, t))
        
        sent = self.sock.send("end_of_pack\n")
    # load_spots()

# class Adxv

def run(params):
    xds_inp = os.path.join(params.xds_dir, "XDS.INP")
    spot_xds = os.path.join(params.xds_dir, "SPOT.XDS")

    spots = idxreflp.SpotXds(spot_xds).indexed_and_unindexed_by_frame_on_detector()
    inputs = get_xdsinp_keyword(xds_inp)

    filename_template = dict(inputs).get("NAME_TEMPLATE_OF_DATA_FRAMES", "")
    if filename_template == "":
        print "Error! Can't find filename from XDS.INP"
        return

    # Start adxv
    adxv = Adxv(params.adxv_bin)
    adxv.start(params.xds_dir)
    type_indexed = adxv.define_spot("blue")
    type_unindexed = adxv.define_spot("red")
    
    num = params.image

    while True:
        print "Showing image %d" % num

        img_file = dataset.template_to_filenames(filename_template, num, num)[0]
        adxv.open_file(img_file)

        uninds = map(lambda x: [x[0], x[1], type_unindexed], spots["unindexed"].get(num, []))
        inds = map(lambda x: [x[0], x[1], type_indexed], spots["indexed"].get(num, []))

        print "Showing %d Indexed spots (blue)" % len(inds)
        print "Showing %d Unindexed spots (red)" % len(uninds)

        adxv.load_spots(inds+uninds)

        time.sleep(1) # wait until adxv finishes process..
        num = int(raw_input("Next image number?: "))

# run()

if __name__ == "__main__":
    import sys

    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args
    
    run(params)
