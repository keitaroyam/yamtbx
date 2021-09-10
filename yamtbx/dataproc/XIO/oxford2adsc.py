#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals
import sys
import struct

#sys.path.append("/Users/plegrand/work/xdsme/XIO")
from . import  XIO

header_fmt = """{
HEADER_BYTES= 1024;
BEAM_CENTER_X=%(beamx).3f;
BEAM_CENTER_Y=%(beamy).3f;
BIN=2x2;
BYTE_ORDER=little_endian;
DATE=%(date)s;
DETECTOR_SN=D2;
DIM=2;
DISTANCE=%(distance).3f;
Data_type=unsigned short int;
KAPPA=0;
OMEGA=0.0000;
OSC_AXIS=PHI;
OSC_RANGE=%(osc_range).4f;
OSC_START=%(osc_start).4f;
PHI=0.0000;
PIXEL_SIZE=%(pixelx).6f;
TIME=0.5;
TWOTHETA=0;
TYPE=unsigned_short;
WAVELENGTH=0.732;
CREV=1;
CCD=TH7899;
ADC=slow;
SIZE1=%(height)d;
SIZE2=%(width)d;
CCD_IMAGE_SATURATION=65535;
}

"""

def _export(filename, format):
    try:
        datacoll = XIO.Collect(filename)
        datacoll.interpretImage()
        #datacoll.lookup_imageRanges()
    except XIO.XIOError:
        print("\nError while trying to acceess %s.\nSorry." % filename)
        sys.exit()
    return datacoll.export(format)
    #print datacoll.export_template(format)


for oxf in sys.argv[1:]:
    datacoll = XIO.Collect(oxf)
    datacoll.interpretImage()
    header = header_fmt % datacoll.export("diffdump")
    
    datacoll.nDigits = 4
    datacoll.setTemplates()
    num = datacoll.getNumber(oxf)
    new_name = datacoll.pythonTemplate % num

    print("%s --> %s" % (oxf, new_name))    
    new = open(new_name,"w")
    data = datacoll.image.getData(clipping=True)
    #new.write("%-1024s" % header)
    new.write(struct.pack("<"+len(data)*"H", *data))
    new.close()
