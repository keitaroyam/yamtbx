from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

# FIXME: serial number...
# reading the data: unpack the data using the ratio.

__version__ = "0.2.0"
__author__ = "Pierre Legrand (pierre.legrand@synchrotron-soleil.fr)"
__date__ = "1-12-2009"
__copyright__ = "Copyright (c) 2009 Pierre Legrand"
__license__ = "New BSD, http://www.opensource.org/licenses/bsd-license.php"

import sys
import struct
import time
import re

PI = 3.1415926535897931
R2D = 180./PI

def date_seconds(time_str):
    "from tupple return seconds"
    try:
        return time.mktime(time.strptime(time_str))
    except ValueError as err:
        print("Warning:", err)
        print("... Using time.time() instead.")
        return time.time()

def pixel_size(detector_type, binning):
    "Determine pixel size from detector_type str."
    det, ver = detector_type.split()
    if det == "SAPPHIRE" and ver[0] == '3':
        return 0.03*binning[0]
    elif det == "RUBY":
        return 0.048*binning[0]
    else:
        return 0.06*binning[0]

# UINT32 unsigned int   = I          = 4 bytes
# INT32 int             = i          = 4 bytes 
# UINT16 unsigned short = H          = 2 bytes
# INT16 short           = h          = 2 bytes
# char                  = c          = 1 byte
# string                = s or p     = 1 byte
# float                 = f          = 4 bytes
# double                = d          = 8 bytes

class Interpreter(object):

    HTD = {
    # The marccd Header Translator Dictionary.
    # To add:
    # + SerialNumber (for the rules) or other unique identifier
           
    'HeaderSize':(['NHEADER'], int),
    'Width':(['NX'], int),
    'Height':(['NY'], int),
    'DateStr':(['TIME'], str),
    'DateSeconds':(['TIME'], date_seconds),
    'Manufacturer':(['MANUFACTURER'], str),
    'DetectorType':(['DETECTOR_TYPE'], str),
    'ImageFormat':(['IMAGE_FORMAT'], str),
    'PixelX':(['DETECTOR_TYPE', 'BINNING'], pixel_size),
    'PixelY':(['DETECTOR_TYPE', 'BINNING'], pixel_size),

    'ExposureTime':(['EXPOSURE_TIME'],lambda x: 60.*float(x)),
    'BeamX':(['BEAM_X'], float),
    'BeamY':(['BEAM_Y'], float),
    'Distance':(['DISTANCE'], float),
    'Wavelength':(['WAVELENGTH'], float),
#    'Message':(['MESSAGE'], lambda x: x.split(';')),
    'PhiStart':(['OSC_START'], float),
    'PhiEnd':(['OSC_END'], float),
    'PhiWidth':(['OSC_WIDTH'], float),
    'OscAxis':(['OSC_AXIS'], lambda x: x.lower()),
    #'EdgeResolution':(['Pixel size X (mm)','Number of pixel (X)',
    #                   'Camera length (mm)','Wavelength'],
    #                    getEdgeResolutionRAXIS),

    # Added keys from Graeme's convention.
    'TwoTheta':(['THETA'], float),
    'SerialNumber':(['SERIAL'], str),  # _VERIFY_
    'EndianType':(['ENDIANTYPE'], str),
    }

    SpecialRules = {
    # No special rules for now
    }

    def getRawHeadDict(self, rawHead):
        EndianType = ">"
        detector_type = rawHead[3:16]
        ascii_head = rawHead[18:205].replace("\r\n"," ")
        re_xdsPar = r"([^ ]+[=])"
        rec_xdsPar = re.compile(re_xdsPar)
        
        #print ascii_head
        lpar = []
        l_s = rec_xdsPar.split(ascii_head)
        len_s = len(l_s)
        if len_s > 1 and len_s % 2:
            for i in range(1,len_s,2):
                lpar.append((l_s[i][:-1],l_s[i+1].strip()))
        RawHeadDict = {}
        for key, val in lpar:
            if key not in RawHeadDict:
                RawHeadDict[key] = val
        RawHeadDict.update({'MESSAGE':'',
                            'TWO_THETA':'0',
                            'MANUFACTURER': 'OXFORD-DIFFRACTION',
                            'DETECTOR_TYPE': detector_type,
                            'IMAGE_FORMAT': 'CRYSALIS'})

        sec1 = int(RawHeadDict['NG'])
        sec2 = int(RawHeadDict['NS'])
        sec3 = int(RawHeadDict['NK'])
        total_r = sec1 + sec2 + sec3
        sec12_d8, sec12_d4, sec1_d8 = (sec1+sec2)/8, (sec1+sec2)/4, sec1/8
        
        tmpstr = rawHead[256:256+total_r]
        tmpint = struct.unpack(total_r/4*'i',tmpstr)
        tmpshort = struct.unpack(total_r/2*'h',tmpstr)
        tmpdble = struct.unpack(total_r/8*'d',tmpstr)

        RawHeadDict['BEAM_X'] = tmpdble[sec12_d8+83]
        RawHeadDict['BEAM_Y'] = tmpdble[sec12_d8+84]
        RawHeadDict['DISTANCE'] = tmpdble[sec12_d8+89]
        RawHeadDict['EXPOSURE_TIME'] = tmpdble[sec1_d8+60]
        RawHeadDict['WAVELENGTH'] = tmpdble[sec1_d8+73]
        # Binning factors (1,1), (2,2) or (4,4)
        RawHeadDict['BINNING'] = tmpshort[:2]
        
        fact1 = tmpdble[sec12_d8+46]*R2D
        angles_order = "omega", "theta", "kappa", "phi", "omega'", "theta'"
        start_angles = tmpint[sec12_d4+71:sec12_d4+77]
        end_angles = tmpint[sec12_d4+81:sec12_d4+87]
        for iii in range(6):
            angle_name = angles_order[iii].upper()
            RawHeadDict[angle_name] = start_angles[iii]*fact1
            if start_angles[iii] != end_angles[iii]:
                RawHeadDict['OSC_AXIS'] = angle_name
                RawHeadDict['OSC_START'] = start_angles[iii]*fact1
                RawHeadDict['OSC_END'] = end_angles[iii]*fact1
                RawHeadDict['OSC_WIDTH'] = RawHeadDict['OSC_END'] - \
                                           RawHeadDict['OSC_START']
        
        RawHeadDict['ENDIANTYPE'] = '>'
        RawHeadDict['SERIAL'] = 'N/A'
        import pprint
        #pprint.pprint(RawHeadDict)
        return  RawHeadDict            

    def getData(self):
        """Read the image bytes. For now only support the 16bits unsigned, and
        uncompressed internaly. Can read compressed file directly (.gz or .Z)."""

        if not self.intCompression:
            _dataSize = self.header['Width']*self.header['Width']
                        
            _image = self.open()
            # Jump over the header
            _image.read(self.header['HeaderSize'])
            # Read the remaining bytes
            _data = _image.read()
            # unpack
            _fmt = self.header['EndianType'] + "H" * _dataSize
            _data = struct.unpack(_fmt, _data)
            
            _image.close()
            assert _dataSize == len(_data)
            
            return _data
        else:
            raise XIOError("Sorry, this image is internaly compressed.")

if __name__ == "__main__":
    from pprint import pprint
    _image = open(sys.argv[1])
    rawHead = _image.read(9000)
    h = getRawHeadDict(rawHead)
   
    headKeys =  [a[0] for a in headerStructure]

    l = ["beam_x","beam_y","two_theta","lambda","distance","phi_start",
            "phi_width","phi_end","pixel_size_x","pixel_size_y","pixel_count_x","pixel_count_y"]
    for k in headKeys: #l:
        print("%s:\t%s" % (k,h[k]))
