"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
import socket
import subprocess
import time
import os

class Adxv:
    def __init__(self, adxv_bin=None, no_adxv_beam_center=True):
        self.adxv_bin = adxv_bin
        self.no_adxv_beam_center = no_adxv_beam_center
        if self.adxv_bin is None: self.adxv_bin = "adxv"

        self.adxv_proc = None # subprocess object
        self.adxv_port = 8100 # adxv's default port. overridden later.

    def start(self, cwd=None):
        adxv_comm = self.adxv_bin + " -socket %d"
        if self.no_adxv_beam_center: adxv_comm += " -no_adxv_beam_center"

        if self.adxv_proc is None or self.adxv_proc.poll() is not None: # None means still running.
            # find available port number
            sock_test = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            sock_test.bind(("localhost", 0))
            self.adxv_port = sock_test.getsockname()[1]
            sock_test.close()
            # start adxv
            self.adxv_proc = subprocess.Popen(adxv_comm%self.adxv_port, shell=True, cwd=cwd)
    # start()

    def open_image(self, imgfile):
        self.start(cwd=os.path.dirname(imgfile))

        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

        for i in xrange(10): # try for 5 seconds.
            try:
                sock.connect(("localhost", self.adxv_port))
                break
            except socket.error:
                time.sleep(.5)
                continue

        sent = sock.send("load_image %s\n"%imgfile)
        if sent == 0:
            raise RuntimeError("adxv load failed! Close adxv and double-click again.")

        sent = sock.send("raise_window Control\n") # raise_window is available from adxv 1.9.9
        sent = sock.send("raise_window Image\n")

        sock.close()
    # open_image()
# class Adxv
