
__version__ = "0.3.4"
__author__ = "Pierre Legrand (pierre.legrand@synchrotron-soleil.fr)"
__date__ = "23-09-2012"
__copyright__ = "Copyright (c) 2007-2012 Pierre Legrand"
__license__ = "New BSD License http://www.opensource.org/licenses/bsd-license.php"

import time
from math import sin, atan

def _datetime(timestr):
    return time.strptime(timestr)

def _dateseconds(timestr):
    return time.mktime(time.strptime(timestr, "%d-%b-%Y %H:%M:%S"))

def getEdgeResolution(DetectorSize, Distance, Wavelength):
       """Calculate EdgeResolution. Rq: doesn't take
       into account the position of direct-beam and the two_theta
       angle of the detector."""
       Distance = float(Distance.split()[-1])
       radius = float(DetectorSize.split()[-1])/2.
       Wavelength = float(Wavelength)
       if abs(Distance) > 0.0:
          return Wavelength/(2*sin(0.5*atan(radius/Distance)))
       else:
          return 0.

def endian(code):
    if code == 'big_endian': return '>'
    else: return '<'

class Interpreter:

    HTD = {
    # The adsc Header Translator Dictionary.
    # Potential problems:
    # - There are multiple SIZE1, SIZE2 instances.
    # = The orientation of SIZE1 and SIZE2 is unknown
    #     Not a problem as long as SIZE1 = SIZE2..

    'ExposureTime':(['SCAN_DET_RELZERO'], lambda x: float(x.split()[1])),
    'BeamX':(['CCD_SPATIAL_BEAM_POSITION'], lambda x: float(x.split()[0])),
    'BeamY':(['CCD_SPATIAL_BEAM_POSITION'], lambda x: float(x.split()[1])),
    'Distance':(['CCD_GONIO_VALUES'], lambda x: float(x.split()[-1])),
    'Wavelength':(['SCAN_WAVELENGTH'], float),
    'PixelX':(['CCD_SPATIAL_DISTORTION_INFO'], lambda x: float(x.split()[-2])),
    'PixelY':(['CCD_SPATIAL_DISTORTION_INFO'], lambda x: float(x.split()[-1])),
    'Width':(['SIZE1'], int),
    'Height':(['SIZE2'], int),
    'Message':(['COMMENT'], str),
    'PhiStart':(['ROTATION'], lambda x: float(x.split()[0])),
    'PhiEnd':(['ROTATION'], lambda x: float(x.split()[1])),
    'PhiWidth':(['ROTATION'], lambda x: float(x.split()[2])),
    'EdgeResolution':(['CCD_DETECTOR_SIZE','CCD_GONIO_VALUES','SCAN_WAVELENGTH'],
                                                            getEdgeResolution),
    'HeaderSize':(['HEADER_BYTES'], int),
    'EndianType':(['BYTE_ORDER'], endian),
    # Added keys from Graeme's convention.
    'TwoTheta':(['CCD_GONIO_VALUES'], lambda x: float(x.split()[1])),   # _FIXME_ Not really here now...
    'SerialNumber':(['CCD_SERIAL_NUMBER'], str),
    'EndianType':(['BYTE_ORDER'], endian),
    'OscAxis':(['SCAN_ROTATION_AXIS_NAME'], lambda x: x.lower()),
    # Date and time
    #'DateTime':(['DATE'], _datetime),
    'DateStr':(['DTREK_DATE_TIME'], str),
    'DateSeconds':(['DTREK_DATE_TIME'], _dateseconds),
    }
    #'EdgeResolution':(['PIXEL_SIZE','SIZE1','DISTANCE','WAVELENGTH'],
    #    getEdgeResolution),

    # Added keys from Graeme's convention.
    #'TwoTheta':(['TWO_THETA'], float),   # _FIXME_ Not really here now...
    #'SerialNumber':(['DETECTOR_SN'], str),

    SpecialRules = {
    # No special rules for now
    }
    
    Identifiers = {
    # Based on Serial Number. Contains (Synchrotron,BLname,DetectorType)
    }
    
    def getRawHeadDict(self, rawHead):
        _lis = rawHead[2:].split("}")[0].split(";\n")[:-1]
        RawHeadDict = dict([par.split("=") for par in _lis])
        RawHeadDict.update({'MESSAGE':'','TWO_THETA':'0'}) # _FIXME_ Not really here now...
        return RawHeadDict
