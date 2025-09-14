# -*- coding: utf-8 -*-

""" XIO plugin for the MarCCD image format.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

__version__ = "0.4.3"
__author__ = "Pierre Legrand (pierre.legrand \at synchrotron-soleil.fr)"
__date__ = "23-09-2012"
__copyright__ = "Copyright (c) 2005-2012 Pierre Legrand"
__license__ = "New BSD, http://www.opensource.org/licenses/bsd-license.php"

import time
import struct
import sys

def extr_time(time_str):
    "from str return tupple"
    try:
        return time.strptime(time_str[12:27].decode(), "%m%d%H%M%Y.%S")
    except ValueError as err:
        print("Warning:", err)
        print("... Using time.localtime() instead.")
        return time.localtime()

def date_seconds(time_str):
    "from tupple return seconds"
    try:
        return time.mktime(extr_time(time_str))
    except ValueError as err:
        print("Warning:", err)
        print("... Using time.time() instead.")
        return time.time()

def date_time(time_str):
    "from str return standard str: 'Wed Oct 28 16:42:12 2009'"
    return time.ctime(date_seconds(time_str))

GETDIST = lambda x, y: (x or y)/1e3
DIVE3 = lambda x: float(x)/1e3
DIVE5 = lambda x: float(x)/1e5
DIVE6 = lambda x: float(x)/1e6
AXIS_CODE = {0: "twotheta", 1:"omega", 2:"chi", 3:"kappa", 4: "phi"}

def get_serial(comment):
    "Try to find the serial string in comments."
    for line in comment.splitlines():
        if line.lower().count(b"serial"):
            try:
                return line.split()[-1].decode()
            except:
                return 'unknown'

def get_edge_resolution_marccd(pixel_x, width, distance,
                                    start_distance, wavelength):
    "Calculate EdgeResolution"
    from math import sin, atan
    distance = distance or start_distance
    if abs(float(distance)) > 0.0:
        rad = 0.5 * float(pixel_x)/1e6 * int(width)
        return float(wavelength/1e5)/(2*sin(0.5*atan(rad/float(distance)*1e3)))
    else:
        return 0.

def get_phi_end(phi_end, phi_start, phi_width):
    "Calculate phiEnd form phiStart and phiWidth"
    if not int(phi_end):
        return (float(phi_start) + float(phi_width))/1000.
    else:
        return float(phi_end)/1000.

HEADER_STRUCTURE = [
('tiff_stuff','1024s'),         # 
('header_type','I'),       # flag for header type (can be used as magic number)
('header_name','16s'),     # header name (MMX)
('header_major_version','I'),   # header_major_version (n.) 
('header_minor_version','I'),   # header_minor_version (.n)
('header_byte_order','I'), # BIG_ENDIAN (Motorola,MIPS); LITTLE_ENDIAN (Intel)
('data_byte_order','I'),   # BIG_ENDIAN (Motorola,MIPS); LITTLE_ENDIAN (Intel)
('header_size','I'),       # in bytes
('frame_type','I'),        # flag for frame type
('magic_number','I'),      # to be used as a flag - usually to indicate new file
('compression_type','I'),  # type of image compression
('compression1','I'),      # compression parameter 1
('compression2','I'),
('compression3','I'),
('compression4','I'),
('compression5','I'),
('compression6','I'),
('nheaders','I'),         # total number of headers
('pixel_count_x','I'),    # number of pixels in one line
('pixel_count_y','I'),    # number of lines in image
('bytes_per_pixel','I'),  # number of bytes per pixel
('bytes_per_line','I'),   # number of pixels between succesive rows
('bits_per_pixel','I'),   # true depth of data, in bits
('data_type','I'),        # (signed,unsigned,float...)
('saturated_pixel','I'),  # value marks pixel as saturated
('sequence','I'),         # TRUE or FAL
('nimages','I'),  # total num. of images - size of each is nfast*(nslow/nimages)
('origin_location','I'),
('detector_orientation','I'),
('view_direction','I'),
('overflow_location','I'),
('overflow_8_bits','I'),
('overflow_16_bits','I'),
('junk3','500s'),
('distance','I'),
('beam_x','I'),
('beam_y','I'),
('integration_time','I'),
('exposure_time','I'),
('readout_time','I'),
('nreads','I'),
('two_theta','i'),
('start_omega','i'),   # 1000*omega angle
('start_chi','i'),     # 1000*chi angle
('start_kappa','i'),   # 1000*kappa angle
('phi_start','i'),
('start_delta','i'),   # 1000*delta angle
('start_gamma','i'),   # 1000*gamma angle
('start_xtal_to_detector','i'), # 1000*distance in mm (dist in um)
('end_twotheta','i'),  # 1000*two_theta angle
('end_omega','i'),     # 1000*omega angle
('end_chi','i'),       # 1000*chi angle
('end_kappa','i'),     # 1000*kappa angle
('phi_end','i'),
('junk7','12s'),
('axis_code','i'),
('phi_width','i'),
('detector_rotx','i'),
('detector_roty','i'),
('detector_rotz','i'),
('total_dose','i'),
('junk8','12s'),
('detector_type','i'),
('pixel_size_x','i'),
('pixel_size_y','i'),
('mean_bias','i'),
('junk9','124s'),
('lambda','I'),
('junk10','100s'),
('filetitle','128s'),
('filepath','128s'),
('filename','64s'),
('acquire_timestamp','32s'),
('header_timestamp','32s'),
('save_timestamp','32s'),
('file_comment','512s'),
('junk11','1132s')]

class Interpreter(object):
    "Dummy class, container for standard Dict and Function."

    HTD = {
    # The marccd Header Translator Dictionary.
    # To add:
    # + SerialNumber (for the rules) or other unique identifier

    'ExposureTime':(['exposure_time'],DIVE3),
    'BeamX':(['beam_x'], DIVE3),
    'BeamY':(['beam_y'], DIVE3),
    'Distance':(['distance','start_xtal_to_detector'], GETDIST),
    'Wavelength':(['lambda'], DIVE5),
    'PixelX':(['pixel_size_x'], DIVE6),
    'PixelY':(['pixel_size_y'], DIVE6),
    'Width':(['pixel_count_x'], int),
    'Height':(['pixel_count_y'], int),
    'Message':(['MESSAGE'], lambda x: x.split(';')),
    'PhiStart':(['phi_start'], DIVE3),
    'PhiEnd':(['phi_end','phi_start','phi_width'], get_phi_end),
    'PhiWidth':(['phi_width'], DIVE3),
    'EdgeResolution':(['pixel_size_x','pixel_count_x','distance',
                       'start_xtal_to_detector','lambda'],
                        get_edge_resolution_marccd),

    # Added keys from Graeme's convention.
    'TwoTheta':(['two_theta'], DIVE3),
    'SerialNumber':(['file_comment'], get_serial),
    'EndianType':(['EndianType'], str),
    'HeaderSize':(['HEADER_BYTES'], int),
    'OscAxis':(['axis_code'], lambda x: AXIS_CODE[x]),
    # Date and time
    'DateStr':(['acquire_timestamp'], date_time),
    'DateSeconds':(['acquire_timestamp'], date_seconds),
    }

    Identifiers = {
    # ESRF info found at
    # http://www.esrf.eu/UsersAndScience/Experiments/MX/Software/PXSOFT/Denzo
    # Based on Serial Number. Contains (Synchrotron,BLname,DetectorType)
    '38':('SLS','X06SA','MarCCD 165'),
    '4':('ESRF','BM14','MarCCD 225'),
    '5':('ESRF','ID23-2','MarCCD 225'),
    '10':('EMBL_HAMBURG','???','MarCCD 225'),
    '12':('SLS','X06SA','MarCCD 225'),
    '21':('EMBL_HAMBURG','???','MarCCD 225'),
    '28':('BESSY','BL14.1','MarCCD 225'),
    '20':('SSRL4','BL11-1','MarCCD 325'),
    }

    SpecialRules = {
    # No special rules for now
    }

    def __init__(self):
        self.raw_head_dict = None

    def getRawHeadDict(self, raw_head, verbose=False):
        "Intepret the binary structure of the marccd header."
        # Get the header endian type
        if struct.unpack('<I', raw_head[1052:1056])[0] == 1234:
            endian_type = "<"
        else:
            endian_type = ">"

        header_structure_fmt = endian_type + \
                            "".join([a[1] for a in HEADER_STRUCTURE])
        header_structure_keys = [a[0] for a in HEADER_STRUCTURE]

        # uncoding from header_structure_fmt
        read_size = struct.calcsize(header_structure_fmt)
        read_unp = list(struct.unpack(header_structure_fmt,
                                         raw_head[:read_size]))

        self.raw_head_dict = {}
        for line in range(len(header_structure_keys)):
            _key = header_structure_keys[line]
            self.raw_head_dict[_key] = read_unp[line]
            if verbose and not _key.count("junk"):
                print("%s ->%s<-" % (_key, read_unp[line]))
        self.raw_head_dict.update({'MESSAGE': '', 'HEADER_BYTES': 4096,
                                           'EndianType': endian_type})
        self.raw_head_dict = dict(((k.encode("utf-8"), self.raw_head_dict[k]) for k in self.raw_head_dict))
        return self.raw_head_dict

def test1():
    "Get the raw header keys from one image"
    image = open(sys.argv[1], "rb")
    raw_head = image.read(9000)
    header = Interpreter().getRawHeadDict(raw_head)

    header_keys =  [a[0] for a in HEADER_STRUCTURE]
    for k in header_keys:
        print("%s:\t%s" % (k, header[k]))

if __name__ == "__main__":
    test1()
