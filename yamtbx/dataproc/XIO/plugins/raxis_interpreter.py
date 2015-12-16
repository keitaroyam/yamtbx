
# FIXME: serial number...
# reading the data: unpack the data using the ratio.

__version__ = "0.1.1"
__author__ = "Pierre Legrand (pierre.legrand@synchrotron-soleil.fr)"
__date__ = "23-09-2012"
__copyright__ = "Copyright (c) 2005-2012 Pierre Legrand"
__license__ = "New BSD, http://www.opensource.org/licenses/bsd-license.php"

import sys
import struct


def getEdgeResolutionRAXIS(PixelX, Width, Distance, Wavelength):
            "Calculate EdgeResolution"
            from math import sin, atan
            if abs(float(Distance)) > 0.0:
                r = 0.5 * float(PixelX) * int(Width)
                return float(Wavelength)/(2*sin(0.5*atan(r/float(Distance))))
            else:
                return 0.

def getPhiEnd(phiEnd,phiStart,phiWidth):
    if not int(phiEnd):
        return (float(phiStart)+float(phiWidth))/1000.
    else:
        return float(phiEnd)/1000.


headerStructure = [
('Device name', '10s'),         # Type of instrument
('Versions', '10s'),
('Crystal name', '20s'),
('Crystal system', '12s'),
('Cell Parameter a', 'f'),
('Cell Parameter b', 'f'),
('Cell Parameter c', 'f'),
('Cell Parameter alpha', 'f'),
('Cell Parameter beta', 'f'),
('Cell Parameter gamma', 'f'),
('Space group', '12s'),         # Space group symbol
('Mosaic spread', 'f'),         # Mosaic spread 
('Memo', '80s'),                # Memo, comments
('res1', '84s'),                # Reserved space for future use

('Date', '12s'),                # Date of measurement
('Operator name', '20s'),       # Username of account collecting image
('X-Ray target', '4s'),         # Type of X-ray target (Cu, Mo, etc.) 
('Wavelength', 'f'),            # X-ray wavelength
('Monochromator', '20s'),       # Type of monochromator
('Monochromator 2theta', 'f'),  # Monochromator 2theta (deg)
('Colliomator', '20s'),         # Collimator size and type
('Filter', '4s'),               # Filter type (Ni, etc.)
('Camera length (mm)', 'f'),    # Crystal-to-detector distance
('Tube Voltage  (KV)', 'f'),    # Generator voltage (kV)
('Tube Current  (mA)', 'f'),    # Generator current (mA)
('X-Ray Focus', '12s'),         # Focus info 
('X-Ray Optics', '80s'),
('Detector Shape Code', 'l'),   # IP shape, 0=flat,1=cylinder
('Weissenberg Coupling', 'f'),  # Weissenberg oscillation 1
('res2', '56s'),                # Reserved space for future use */

('Mounted axis', '4s'),         # Crystal mount axis closest to spindle axis   
('X-ray beam axis', '4s'),      # Crystal mount axis closest to beam axis
('Phi 0', 'f'),                 # Datum phi angle (deg)
('Phi start', 'f'),             # Phi oscillation start angle (deg)
('Phi end', 'f'),               # Phi oscillation end angle (deg)
('Oscillation number ', 'l'),   # Frame number
('Exposure time (min)', 'f'),   # Exposure time (min)
('Direct beam X (pixel)', 'f'), # Direct beam X position
('Direct beam Z (pixel)', 'f'), # Direct beam Z position
('Omega', 'f'),                 # Goniostat angle omega
('Chi', 'f'),                   # Goniostat angle chi
('Theta', 'f'),                 # Goniostat angle 2theta
('Mu', 'f'),                    # Spindle inclination angle
('res3', '204s'),               # Reserved space for future use
                                # This space is now used for storing the scan
                                # templates information - tlh, 01 Feb 1999

('Number of pixel (X)', 'l'),   # number of pixels in X direction 
('Number of pixel (Z)', 'l'),   # number of pixels in Z direction
('Pixel size X (mm)', 'f'),     # size of pixel in X direction (mm) 
('Pixel size Z (mm)', 'f'),     # size of pixel in Z direction (mm)
('Record length (byte)', 'l'),  # Header record length (bytes)
('Number of records', 'l'),     # number of records (lines) in inage
('Start line (Z)', 'l'),        # starting line number
('IP number', 'l'),             # IP number 
('Ratio', 'f'),                 # photomultiplier output hi/lo ratio
('Fading function No.1', 'f'),  # fading time, end of exposure to start of read
('Fading function No.2', 'f'),  # fading time, end of exposure to end of read
('Host Type', '10s'),           # type of computer (IRIS, VAX) => endian
('IP Type', '10s'),             # type of IP 
('Image Direction X', 'l'),     # horizontal scanning code: 0=left->right, 1=>right->left
('Image Direction Z', 'l'),     # vertical scanning code: 0=down->up, 1=up->down
('Image Direction XZ', 'l'),    # front/back scanning code: 0=front, 1=back

('shft', 'f'),                  # Pixel shift, R-AXIS V 
('ineo', 'f'),                  # Intensity ratio E/O R-AXIS V 
('majc', 'l'),                  # Magic number to indicate next values are legit
('naxs', 'l'),                  # Number of goniometer axes
('gvec', '15f'),                # Goniometer axis vectors
('gst', '5f'),                  # Start angles for each of 5 axes
('gend', '5f'),                 # End angles for each of 5 axes 
('goff', '5f'),                 # Offset values for each of 5 axes
('saxs', 'l'),                  # Which axis is the scan axis?
('gnom', '40s'),                # Names of the axes (space or comma separated?) 

# Most of below is be program dependent.  Different programs use
# this part of the header for different things.  So it is essentially 
# a big "common block" area for dumping transient information.

('file', '16s'),                # 
('cmnt', '20s'),                # 
('smpl', '20s'),                # 
('iext', 'l'),                  # 
('reso', 'l'),                  # 
('save', 'l'),                  # 
('dint', 'l'),                  # 
('byte', 'l'),                  # 
('init', 'l'),                  # 
('ipus', 'l'),                  # 
('dexp', 'l'),                  # 
('expn', 'l'),                  # 
('posx', '20l'),                # 
('posy', '20l'),                # 
('xray', 'i'),                  # 
('res5', '768s'),               # Reserved space for future use
]

class Interpreter:

    HTD = {
    # The marccd Header Translator Dictionary.
    # To add:
    # + SerialNumber (for the rules) or other unique identifier
                
    'ExposureTime':(['Exposure time (min)'],lambda x: 60.*float(x)),
    'BeamX':(['Direct beam X (pixel)', 'Pixel size X (mm)', 'Number of pixel (X)'], lambda x, y, z: (z-x)*y),
    'BeamY':(['Direct beam Z (pixel)', 'Pixel size Z (mm)', 'Number of pixel (Z)'], lambda x, y, z: (z-x)*y),
    'Distance':(['Camera length (mm)'], float),
    'Wavelength':(['Wavelength'], float),
    'PixelX':(['Pixel size X (mm)'], float),
    'PixelY':(['Pixel size Z (mm)'], float),
    'Width':(['Number of pixel (X)'], int),
    'Height':(['Number of pixel (X)'], int),
    'Message':(['MESSAGE'], lambda x: x.split(';')),
    'PhiStart':(['Phi start'], float),
    'PhiEnd':(['Phi end'], float),
    'PhiWidth':(['Phi end','Phi start'], lambda x, y: x-y),
    'EdgeResolution':(['Pixel size X (mm)','Number of pixel (X)',
                       'Camera length (mm)','Wavelength'],
                        getEdgeResolutionRAXIS),

    'TwoTheta':(['Theta'], float),
    'SerialNumber':(['IP number'], str),  # _VERIFY_
    'EndianType':(['EndianType'], str),
    'HeaderSize':(['Record length (byte)'], int),
    'OscAxis':(['OSC_AXIS'], lambda x: x.lower()),
    }

    SpecialRules = {
    # No special rules for now
    }

    def getRawHeadDict(self, rawHead):
        
        # It seems that
        EndianType = ">"
        headerStructureFmt = EndianType + \
                             "".join([fmt for key,fmt in headerStructure])
        headerStructureKeys = [key for key,fmt in headerStructure]
        
        # uncoding from headerStructureFmt
        read_size = struct.calcsize(headerStructureFmt)
        read_unp = struct.unpack(headerStructureFmt, rawHead[:read_size])

        RawHeadDict = dict(zip(headerStructureKeys, read_unp))
        import pprint
        #pprint.pprint(RawHeadDict)
        RawHeadDict.update({'MESSAGE':'',
	                    'EndianType':EndianType,
			    'OSC_AXIS': "PHI"})
        return RawHeadDict

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
            raise XIOError, "Sorry, this image is internaly compressed."

if __name__ == "__main__":
    from pprint import pprint
    _image = open(sys.argv[1])
    rawHead = _image.read(9000)
    h = getRawHeadDict(rawHead)
   
    headKeys =  [a[0] for a in headerStructure]

    l = ["beam_x","beam_y","two_theta","lambda","distance","phi_start",
            "phi_width","phi_end","pixel_size_x","pixel_size_y","pixel_count_x","pixel_count_y"]
    for k in headKeys: #l:
        print "%s:\t%s" % (k,h[k])
