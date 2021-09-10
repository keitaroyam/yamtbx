from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

__version__ = "0.3.3"
__author__ = "Pierre Legrand (pierre.legrand@synchrotron-soleil.fr)"
__date__ = "23-09-2012"
__copyright__ = "Copyright (c) 2005-2012 Pierre Legrand"
__license__ = "New BSD, http://www.opensource.org/licenses/bsd-license.php"

import sys
import struct
import time

VERBOSE = False
divE3 = lambda x: float(x)/1e3
    
def getEdgeResolutionMAR345(PixelX, Width, Distance, Wavelength):
            "Calculate EdgeResolution:"
            from math import sin, atan
            if abs(Distance) > 0.0:
                r = 0.5 * PixelX/1.e3 * Width
                return Wavelength/(2*sin(0.5*atan(r/Distance)))
            else:
                return 0.
def get_osc_axis(omega_osc, phi_osc):
    if omega_osc != 0.:
        return "omega"
    elif phi_osc != 0.:
        return "phi"
    else:
        return "undefined"

def date_seconds(time_str):
    "from tupple return seconds"
    try:
        return time.mktime(time.strptime(time_str))
    except ValueError as err:
        print("Warning:", err)
        print("... Using time.time() instead.")
        return time.time()

# Contains (parameterName,binaryEncodedType [,defaultsValue])
# Default value is optional

binHeaderStructure = [
# Initial binaray coded structure
('byteorder','I',1234),         # Always = 1234
('size','i',0),                 # No. of pixels in 1 dimension
('high','I',0),                 # No. high intensity pixels
('format','1s',1),              # Image format
('mode','1s',1),                # Exposure mode 
('pixels','I',0),               # No. of pixels in image
('pixel_length','f',150.),      # Length of 1 pixel
('pixel_height','f',150.),      # Height of 1 pixel
('wave','f',1.541789),          # Wavelength [Ang.]
('dist','f',70.0),              # Distance [mm]
('phibeg','f',0.),              # Starting PHI
('phiend','f',0.),              # Ending   PHI
('omebeg','f',0.),              # Starting Omega
('omeend','f',0.),              # Ending   Omega
('chi','f',0.),                 # Chi
('theta','f',0.),               # Two theta

# Scanner specific things  
('version','8s'),               # Program version
('program','16s'),              # Program name
('scanner','i',1),              # Scanner serial no
('adc_A','I',-1),               # Offset from channel A of ADC
('adc_B','I',-1),               # Offset from channel B of ADC
('add_A','I',0),                # ADD to channel A of ADC
('add_B','I',0),                # ADD to channel B of ADC
('gap','8s'),                   # GAP1+2 position seen by controller
('multiplier','f',1.0),         # Multiplication factor
('xcen','f',600.),              # Center x of transf. image
('ycen','f',600.),              # Center y of transf. image
('roff','f',0.0),               # Radial offset
('toff','f',0.0),               # Tangential offset
('gain','f',1.0),               # Gain of detector

# Experimental conditions for this image
('time','f',0.0),               # Exposure time in secs
('dosebeg','f',0.0),            # Dose at start of expose
('doseend','f',0.0),            # Dose at end   of expose
('dosemin','f',0.0),            # Min. dose during expose
('dosemax','f',0.0),            # Max. dose during expose
('doseavg','f',0.0),            # Avg. dose during expose
('dosesig','f',0.0),            # Sig. dose during expose
('resol','f'),                  # Max. resolution
('phiosc','I',0),               # Phi oscillations
('omeosc','I',0),               # Omega oscillations
('dosen','I',0),                # No. of X-ray readings  
# Generator settings
#('',''),
]

asciiHeaderStructure = {
'program':('PROGRAM', lambda x: x[0]),
'version':('PROGRAM', lambda x: x[2]),
'roff':('OFFSET', lambda x: float(x[x.index('ROFF')+1])),
'toff':('OFFSET', lambda x: float(x[x.index('TOFF')+1])),
'gap':('GAP', lambda x: int(x[0])),
'adc_A':('ADC', lambda x: int(x[x.index('A')+1])),
'adc_B':('ADC', lambda x: int(x[x.index('B')+1])),
'add_A':('ADC', lambda x: int(x[x.index('ADD_A')+1])),
'add_B':('ADC', lambda x: int(x[x.index('ADD_B')+1])),
'multiplier':('MULTIPLIER', lambda x: float(x[0])),
'gain':('GAIN', lambda x: float(x[0])),
'doseavg':('COUNTS', lambda x: float(x[x.index('AVE')+1])),
'dosesig':('COUNTS', lambda x: float(x[x.index('SIG')+1])),
'dosebeg':('COUNTS', lambda x: float(x[x.index('START')+1])),
'doseend':('COUNTS', lambda x: float(x[x.index('END')+1])),
'dosen':('COUNTS', lambda x: float(x[x.index('NMEAS')+1])),
'dosemin':('COUNTS', lambda x: float(x[x.index('MIN')+1])),
'dosemax':('COUNTS', lambda x: float(x[x.index('MAX')+1])),
'mode':('MODE', lambda x: int(x[0]=="TIME")),
'dist':('DISTANCE', lambda x: float(x[0])),
'pixel_length':('PIXEL', lambda x: float(x[x.index('LENGTH')+1])),
'pixel_height':('PIXEL', lambda x: float(x[x.index('HEIGHT')+1])),
'scanner':('SCANNER', lambda x: int(x[0])),
'high':('HIGH', lambda x: int(x[0])),
'date':('DATE', lambda x: " ".join(x)),
'remark':('REMARK', lambda x: " ".join(x)),
'format':('FORMAT', lambda x: int(x[1]=="PCK")+2*int(x[1]=="SPIRAL")),
'size':('FORMAT', lambda x: int(x[0])),
'pixels':('FORMAT', lambda x: int(x[2])),
'wave':('WAVELENGTH', lambda x: float(x[0])),
'polar':('MONOCHROMATOR', lambda x: float(x[x.index('POLAR')+1])),
'filter':('MONOCHROMATOR', lambda x: x[0]),
'phibeg':('PHI', lambda x: float(x[x.index('START')+1])),
'phiend':('PHI', lambda x: float(x[x.index('END')+1])),
'phiosc':('PHI', lambda x: int(x[x.index('OSC')+1])),
'omebeg':('OMEGA', lambda x: float(x[x.index('START')+1])),
'omeend':('OMEGA', lambda x: float(x[x.index('END')+1])),
'omeosc':('OMEGA', lambda x: int(x[x.index('OSC')+1])),
'theta':('TWOTHETA', lambda x: float(x[0])),
'chi':('CHI', lambda x: float(x[0])),
'resol':('RESOLUTION', lambda x: float(x[0])),
'time':('TIME', lambda x: float(x[0])),
'xcen':('CENTER', lambda x: float(x[x.index('X')+1])),
'ycen':('CENTER', lambda x: float(x[x.index('Y')+1])),
'slitx':('COLLIMATOR', lambda x: float(x[x.index('WIDTH')+1])),
'slity':('COLLIMATOR', lambda x: float(x[x.index('HEIGHT')+1])),
'mA':('GENERATOR', lambda x: float(x[x.index('mA')+1])),
'kV':('GENERATOR', lambda x: float(x[x.index('kV')+1])),
'source':('GENERATOR', lambda x: x[0]),
'valmin':('INTENSITY', lambda x: float(x[x.index('MIN')+1])),
'valmax':('INTENSITY', lambda x: float(x[x.index('MAX')+1])),
'valavg':('INTENSITY', lambda x: float(x[x.index('AVE')+1])),
'valsig':('INTENSITY', lambda x: float(x[x.index('SIG')+1])),
'histmax':('HISTOGRAM', lambda x: int(x[x.index('MAX')+1])),
'histbeg':('HISTOGRAM', lambda x: int(x[x.index('START')+1])),
'histend':('HISTOGRAM', lambda x: int(x[x.index('END')+1])),
}

class Interpreter(object):

    HTD = {
    # The mar345 Header Translator Dictionary.

    'ExposureTime':(['time'],float),
    'BeamX':(['xcen','pixel_length'], lambda x, y: x*y/1.e3),
    'BeamY':(['ycen','pixel_height'], lambda x, y: x*y/1.e3),
    'Distance':(['dist'], float),
    'Wavelength':(['wave'], float),
    'PixelX':(['pixel_length'], divE3),
    'PixelY':(['pixel_height'], divE3),
    'Width':(['size'], int),
    'Height':(['size'], int),
    'Message':(['MESSAGE'], lambda x: x.split(';')),
    'PhiStart':(['phibeg'], float),
    'PhiEnd':(['phiend'], float),
    'PhiWidth':(['phiend','phibeg'], lambda x, y: x-y),
    'EdgeResolution':(['pixel_length','size','dist','wave'],
        getEdgeResolutionMAR345),

    # Added keys from Graeme's convention.
    'TwoTheta':(['theta'], float),
    'SerialNumber':(['scanner'], str),
    'OscAxis':(['omeosc','phiosc'], get_osc_axis),
    'HeaderSize':(['HEADER_BYTES'], int),
    'EndianType':(['EndianType'], str),
    #'DateTime':(['DATE'], date_time),
    'DateStr':(['date'], str),
    'DateSeconds':(['date'], date_seconds),
    }

    Identifiers = {
    # ESRF info found at
    # http://www.esrf.eu/UsersAndScience/Experiments/MX/Software/PXSOFT/Denzo
    # Based on Serial Number. Contains (Synchrotron,BLname,DetectorType)
    '50':('EMBL_HAMBURG','BW7A','Mar345'),
    '92':('MARSEILLE','AFMB','Mar345'),
    }

    SpecialRules = {
    # No special rules for now
    }

    def getRawHeadDict(self, rawHead):

        # Initialise using optional default values
        initvalues = [(a[0],a[2]) for a in binHeaderStructure if len(a) == 3]
        RawHeadDict = dict(initvalues)
        
        # Get the header endian type
        if struct.unpack('<I',rawHead[:4])[0] == 1234:
            EndianType = "<"
        else:
            EndianType = ">"
        
        # Try to interpret the first 128 bits containing only 
        # 16 * 4bits encoded values. 
        fours = [rawHead[4*a:4*(a+1)] for a in  range(32)]
        for n in range(16):
            fmt = EndianType + binHeaderStructure[n][1]
            packed = fours[n][:struct.calcsize(fmt)]
            RawHeadDict[binHeaderStructure[n][0]] = struct.unpack(fmt,packed)[0]
        
        _dic = {}
        _lis = rawHead[192:].split("END OF HEADER")[0].split("\n")[:-1]
        for a in [a for a in [x.split() for x in _lis] if a]:
            if a[0] not in _dic: _dic.update({a[0]:a[1:]})
            else: _dic[a[0]].extend(a[1:])
        
        for k in asciiHeaderStructure:
            keyw, func = asciiHeaderStructure[k]
            if keyw in _dic:
                try:
                    RawHeadDict[k] = func(_dic[keyw])
                except ValueError:
                    if VERBOSE:
                        print("WARNING: Can't interpret header key: %s" % k)
        
        # Add other needed information to RawHeadDict
        RawHeadDict.update({'MESSAGE':'','HEADER_BYTES':4096,
                                         'EndianType':EndianType})
        
        return RawHeadDict
