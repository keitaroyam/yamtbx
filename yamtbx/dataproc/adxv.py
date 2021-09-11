"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals
from builtins import range
from builtins import object
import socket
import subprocess
import time
import os
import getpass
import tempfile

class Adxv(object):
    def __init__(self, adxv_bin=None, no_adxv_beam_center=True):
        self.adxv_bin = adxv_bin
        self.no_adxv_beam_center = no_adxv_beam_center
        if self.adxv_bin is None: self.adxv_bin = "adxv"

        self.adxv_proc = None # subprocess object
        self.adxv_port = 8100 # adxv's default port. overridden later.
        self.sock = None

        self.spot_type_counter = -1
    # __init__()

    def start(self, cwd=None):
        adxv_comm = self.adxv_bin + " -socket %d"
        if self.no_adxv_beam_center: adxv_comm += " -no_adxv_beam_center"

        if not self.is_alive():
            # find available port number
            sock_test = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            sock_test.bind(("localhost", 0))
            self.adxv_port = sock_test.getsockname()[1]
            sock_test.close()
            # start adxv
            self.adxv_proc = subprocess.Popen(adxv_comm%self.adxv_port, shell=True, cwd=cwd)

            for i in range(10): # try for 5 seconds.
                try: 
                    self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM) # On OSX(?), need to re-create object when failed
                    self.sock.connect(("localhost", self.adxv_port))
                    break
                except socket.error:
                    time.sleep(.5)
                    continue

    # start()

    def is_alive(self):
        return self.adxv_proc is not None and self.adxv_proc.poll() is None  # None means still running.
    # is_alive()

    def open_image(self, imgfile, raise_window=True):
        self.start(cwd=os.path.dirname(imgfile))


        sent = self.sock.send(b"load_image %s\n"%imgfile.encode("utf-8"))
        if sent == 0:
            raise RuntimeError("adxv load failed! Close adxv and double-click again.")

        if raise_window:
            sent = self.sock.send("raise_window Control\n") # raise_window is available from adxv 1.9.9
            sent = self.sock.send("raise_window Image\n")

        #sock.close()
    # open_image()

    def open_hdf5(self, h5file, frameno_or_path, tmpdir=None, raise_window=True, binning=1):
        from yamtbx.dataproc import eiger

        if tmpdir is None:
            tmpdir = "/dev/shm" if os.path.isdir("/dev/shm") else tempfile.gettempdir()

        imgfile = os.path.join(tmpdir, "adxvtmp-%s-%s.cbf"%(getpass.getuser(), os.getpid()))
        eiger.extract_to_minicbf(h5file, frameno_or_path, imgfile, binning=binning)
        self.open_image(imgfile, raise_window=raise_window)
    # open_hdf5()

    def define_spot(self, color, radius=0, box=0):
        self.spot_type_counter += 1
        sent = self.sock.send("box %d %d\n" % (box,box)) # seems ignored?
        sent = self.sock.send("define_type %d color %s radius %d\n"%(self.spot_type_counter, color, radius))
        print(sent)
        if sent == 0:
            print("define_spot failed!")

        sent = self.sock.send("box 20 20\n")
        return self.spot_type_counter
    # define_spot()

    def load_spots(self, spots):
        if len(spots) == 0:
            return

        sent = self.sock.send("load_spots %d\n" % len(spots))

        for x, y, t in spots:
            sent = self.sock.send("%.2f %.2f %d\n" % (x, y, t))
        
        sent = self.sock.send("end_of_pack\n")
    # load_spots()

# class Adxv
