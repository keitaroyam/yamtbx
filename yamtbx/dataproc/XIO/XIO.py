#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Maintained by P.Legrand
 28th March 2005

 TODO:
 - if given, get the goniometer angles (list of 4 angles omega, kappa/chi, phi, theta).
 - What needs to be unicode compatible (type(val) == str or unicode)???

New BSD License http://www.opensource.org/licenses/bsd-license.php
"""

__version__ = "0.4.5"
__author__ = "Pierre Legrand (pierre.legrand \at synchrotron-soleil.fr)"
__date__ = "7-12-2009"
__copyright__ = "Copyright (c) 2005-2009 Pierre Legrand"
__license__ = "New BSD License www.opensource.org/licenses/bsd-license.php"

#
# Standard modules
#

import os
import sys
import struct
import re
import time
import datetime

PLUGIN_DIR_NAME = "plugins"
VERBOSE = 0

# Defines General regexp for the complete image names:
# dirPath + imageName + externalCompression

RE_DIRPATH =       r"(.*/+)?"               # opt. directory path
RE_IMAGE_NAME1 =    r"(.*)_(\d+)\.(\w+)"    # image name with prefix_num.ext
RE_IMAGE_NAME2 =    r"(.*)\.(\d{3,5})"      # image name with prefixnum.ext
RE_IMAGE_NAME3 =    r"(.*)(\d{3,5})\.(\w+)" # image name with prefixnum.ext
RE_EXT_COMPRESSION = r"(\.gz|\.Z|\.bz2)?\Z"   # opt. ending compression ext
RE_EXTCOMPRESSED =   RE_EXT_COMPRESSION.replace('?','+')  # compression
RE_FULLIMAGENAME1 =  RE_DIRPATH + RE_IMAGE_NAME1 + RE_EXT_COMPRESSION
RE_FULLIMAGENAME2 =  RE_DIRPATH + RE_IMAGE_NAME2 + RE_EXT_COMPRESSION
RE_FULLIMAGENAME3 =  RE_DIRPATH + RE_IMAGE_NAME3 + RE_EXT_COMPRESSION

REC_DIRPATH =         re.compile(RE_DIRPATH)
REC_EXTCOMPRESSION =  re.compile(RE_EXT_COMPRESSION)
REC_EXTCOMPRESSED =   re.compile(RE_EXTCOMPRESSED)
REC_FULLIMAGENAME1 =  re.compile(RE_FULLIMAGENAME1)
REC_FULLIMAGENAME2 =  re.compile(RE_FULLIMAGENAME2)
REC_FULLIMAGENAME3 =  re.compile(RE_FULLIMAGENAME3)

def list_of_string(arg):
    "Return True if all the component of the list are of string type."
    return reduce(lambda a, b: a and b, map(lambda s: type(s) == str, arg))

def isExtCompressed(filename):
    "Tells from the filename suffix if the file is externaly compressed"
    if REC_EXTCOMPRESSED.search(filename):
        return True
    else:
        return False

def uncompressedName(filename):
    "Remove the compression extention in filename"
    if isExtCompressed(filename):
        return filename[:filename.rindex(".")]
    else: return filename

def remove_redundancy(seq):
    "Fast way to remove the redundant elements in any sorted sequence."
    _seq_shift = seq[1:]
    _seq_shift.append(None)
    return [E for E, C in zip(seq, map(cmp, seq, _seq_shift)) if C]

def add_reduce(seq):
    "Like Numeric.add.reduce(array)..."
    from operator import add
    return reduce(add, seq)

def importName(moduleName, namedObject):
    """ Import a named object from a module. To do something like:
    from moduleName import namedObject
    but with moduleName and namedObject are runtime computed expressions.
    NOTE: It's more secured solution than just using __import__
    """
    #
    try:
        module = __import__(moduleName, globals(), locals(), namedObject)
    except ImportError:
        "Warning, can't import plugins..."
        return None
    return vars(module)[namedObject]


class XIOError(Exception):
    """This level of exception raises a recoverable error which can be fixed.
    """

class Image:
    """Defines a generic X-ray diffraction image file.
    The constructor tries first to determine the image type (using guessType()),
    and then the internal compression (using guessIntCompression()).

    Then, if a type specific plugin exist for that image type, the header will
    be interpreted and stored in a dictionary with standard keys (as difined by
    the DNA DiffractionImage module).
    The constructor only reads the 9000 first bytes which are supposed to
    contain the header data for all images formats.

    Files externaly compressed with gzip (.gz) or compress (.Z) or bzip2 (.bz2)
    can be read transparently.

    In that case, only the first 9000 bytes are read and uncompressed.

    NOTE: To avoid confusion between "internal" compression type like PCK and
    external compression type like gzip, bzip2, compress, naming as been adapted
    with extCompres... and intCompres...

    #>>> im = Image("e1ip_0_001.mar2000")
    #>>> assert im.type == 'mar'
    #>>> assert im.intCompression == 'pck'
    #>>> assert im.header['ExposureTime'] == 90.0
    #>>> assert im.header['PixelX'] == 0.150

    >>> im = Image("test/ropBr1_1_001.mccd.gz")
    >>> assert im.type == 'marccd'
    >>> assert im.intCompression == None
    >>> assert im.header['ExposureTime'] == 30.0
    >>> assert im.header['BeamX'] == 112.056
    >>> assert im.header['Wavelength'] == 0.91874

    >>> im = Image("test/before_1_001.img.gz")
    >>> assert im.type == 'adsc'
    >>> assert im.intCompression == None
    >>> assert im.header['ExposureTime'] == 1.0
    >>> assert im.header['SerialNumber'] ==  '420'
    >>> assert im.header['PhiWidth'] == 1.0
    >>> assert im.header['Width'] == 2304
    """

    def __init__(self, imageFileName, doInterpret=True):

        self.fileName = imageFileName
        #
        # Add shared library path to sys.path to be used with the "
        # plugin" __import__ based mecanism. __file__ work for
        # python version > 2.2.
        try:
            sys.path.append(os.path.join(os.path.split(__file__)[0],
                                          PLUGIN_DIR_NAME))
            #print "DEBUG: __file__ = ",os.path.split(__file__)[0]
        except:
            pass
            # FIXME Find an alternative method to find the path... for older
            # python version
            # sys.path.append(os.path.join(os.path.split(".")[0])
        try:
            _image = self.open()
            self.rawHead = _image.read(9000)
            _image.close()
        except XIOError, msg:
            print msg
            raise XIOError("Can't read compressed file %s" \
                                             % (imageFileName))
        if len(self.rawHead) < 9000:
            msg = "File %s seems incomplete. Check its size." % imageFileName
            raise XIOError(msg)

        #
        self.type = ""
        self.intCompression = None
        self.interpreter = None
        self.header = {}
        self.detModel = None
        self.beamline = None
        #
        # Methodes that fill the attributes [type, header, compression]
        self.guessType()
        self.guessIntCompression()
        if doInterpret:
            self.headerInterpreter()
        try:
            self.guessDetModel()
        except:
            pass
        try:
            self.guessBeamline()
        except:
            pass

    def open(self):
        "Open with the appropriate function (in read only)."

        # Get the real path (not the eventual link name).
        _realpath = os.path.realpath(self.fileName)
        if isExtCompressed(_realpath):
            if _realpath[-3:].lower() == ".gz":
                try:
                    import gzip
                except ImportError:
                    raise XIOError("Can't read compressed file %s." \
                                   % (self.fileName) + " (No module gzip!)")
                try:
                    return gzip.open(self.fileName)
                except IOError:
                    raise XIOError("Can't open compressed file %s" \
                                            % (self.fileName))
            elif _realpath[-4:].lower() == ".bz2":
                try:
                    import bz2
                except ImportError:
                    raise XIOError("Can't read compressed file %s." \
                                   % (self.fileName) + " (No module bz2!)")
                try:
                    print "trying bz2"
                    return bz2.BZ2File(self.fileName)
                except IOError, _mess:
                    print "IOError", _mess
                    raise XIOError, "argh! Can't open compressed file %s" \
                                            % (self.fileName)
        else:
            return open(self.fileName)

    def guessType(self):
        """Try to guess the image detector type.
           Detector format recognized:
           adsc, mar, marccd, raxis"""

        if self.fileName.endswith(".h5"):
            import h5py
            h5 = h5py.File(self.fileName, "r")
            if "entry/instrument/detector/description" in h5 and "Eiger" in h5["entry/instrument/detector/description"].value:
                self.type = "eiger_hdf5"
                return self.type
            elif "/entry/instrument/detector/detectorSpecific/eiger_fw_version" in h5:
                self.type = "eiger_hdf5" # probably eiger2..
                return self.type
            return

        #
        # Test to identify Mar345 or Mar555 header
        # Look for a swap byte marker
        m1 = struct.unpack('<I', self.rawHead[:4])[0]
        m2 = struct.unpack('>I', self.rawHead[:4])[0]
        t1, t2 = 1234, 4321
        martest1 = (m1 == t1 or m1 == t2 or m2 == t1 or m2 == t2)

        if (martest1 and
            self.rawHead.count("mar research     ") and
            self.rawHead.count("mar555")):
            self.type = "mar555"
            return self.type

        elif (martest1 and self.rawHead.count("mar research     ")):
            self.type = "mar"
            return self.type

        # Test to identify ADSC header
        elif  self.rawHead[:15] == "{\nHEADER_BYTES=" and \
                self.rawHead.count(";\nPIXEL_SIZE="):
            self.type = "adsc"
            return self.type

        # Test to identify MSC ccd header
        elif  self.rawHead[:15] == "{\nHEADER_BYTES=" and \
                self.rawHead.count(";\nCCD_DETECTOR_DESCRIPTION="):
            self.type = "mscccd"
            return self.type

        # Test to identify MarCCD header
        elif self.rawHead[0:3] == "II*" and \
                struct.unpack('<I', self.rawHead[1024:1028])[0] == 2L and \
                self.rawHead[1028:1031] == "MMX" :
            self.type = "marccd"
            return self.type

        # Test to identify imageCIF or CBF header
        elif self.rawHead.count("loop_") >= 3 and \
             self.rawHead.count("data_image_"):
            self.type = "cbf"
            return self.type

        # Test to identify miniCIF (PILATUS)
        elif self.rawHead[0:7] == "###CBF:" and \
                (self.rawHead.count("PILATUS") or self.rawHead.count("Eiger")):
            self.type = "minicbf"
            return self.type

        # Test to identify RAXIS header
        elif self.rawHead[0:5] == "RAXIS":
            self.type = "raxis"
            return self.type

        # Test to identify Oxford Diffraction header
        elif "OD " == self.rawHead[:3]:
            self.type = "oxford"
            return self.type

        # Failled to identify any image type :-(
        else:
            print self.rawHead[:20]
            self.type = "unknown"
            return self.type

    def guessDetModel(self):
        if (not "Width" in self.header) or (not "PixelX" in self.header):
            return
        _size = self.header["Width"]*self.header["PixelX"]
        #print "SIZE= %.2f" % _size
        if self.type == "marccd":
            if 167. > _size > 160.:
                self.detModel = "MarCCD 165"
            elif 228. > _size > 220.:
                self.detModel = "MarCCD 225"
            elif 303. > _size > 297.:
                self.detModel = "MarCCD 300"
            elif 328. > _size > 319.:
                self.detModel = "MarCCD 325"
        elif self.type == "adsc":
            if 190. > _size > 187.:
                self.detModel = "ADSC Q4"
            elif 212. > _size > 207.:
                self.detModel = "ADSC Q210"
            elif 317. > _size > 310.:
                self.detModel = "ADSC 315"

    def guessBeamline(self):
        if hasattr(self.interpreter, "Identifiers"):
            if self.header["SerialNumber"] in self.interpreter.Identifiers:
                self.beamline = self.interpreter.Identifiers[
                                              self.header["SerialNumber"]]
        else:
            return

    def guessIntCompression(self):
        """Try to guess the image detector internal compression type.
        For now the only known used (by MAR and XDS) internal compression
        format is the pck CCP4 packed ."""

        if self.rawHead.count("CCP4 packed image"):
            self.intCompression = "pck"
        elif self.rawHead[0:7] == "###CBF:":
            self.intCompression = "cbf"
        elif "COMPRESSION=" in self.rawHead[:50]:
            self.intCompression = "CRYSALIS"
        return  self.intCompression

    def headerInterpreter(self):
        """Try to interpret the header, and return a dictionary with the header
        information"""
        #print self.type
        interpreterClass = importName("XIO.plugins.%s_interpreter" % \
                                       self.type, "Interpreter")
        if not  interpreterClass:
            interpreterClass = importName("plugins.%s_interpreter" % \
                                       self.type, "Interpreter")
        if not interpreterClass:
            raise XIOError, "Can't import %s interperter" % (self.type)

        # Rules are serial number (or other identifier) based rules
        # To be added
        # Special = interpreter.SpecialRules
        #
        self.interpreter = interpreterClass()
        if self.type == "eiger_hdf5":
            self.RawHeadDict = self.interpreter.getRawHeadDict(self.fileName)
        else:
            self.RawHeadDict = self.interpreter.getRawHeadDict(self.rawHead)
        #VERBOSE = True
        for k in self.interpreter.HTD.keys():
            args, func = self.interpreter.HTD[k]
            #self.header[k] = apply(func, map(self.RawHeadDict.get,args))
            if args[0] in self.RawHeadDict:
                try:
                    self.header[k] = func(*map(self.RawHeadDict.get, args))
                except ValueError:
                    self.header[k] = 0.
                    if VERBOSE:
                        print "WARNING: Can't interpret header KEY %s" % k
        # Check consistancy of beam center coordinates (should be in mm).
        # with pixel size and number...
        # Some time the beam center is expressed in pixels rather than in mm.
        if (self.header["BeamX"] > self.header["Width"]*self.header["PixelX"])\
           and \
           (self.header["BeamX"] > self.header["Width"]*self.header["PixelX"]):
            self.header["BeamX"] = self.header["BeamX"]*self.header["PixelX"]
            self.header["BeamY"] = self.header["BeamY"]*self.header["PixelY"]
        self.header["ImageType"] = self.type

        if self.type == "adsc": # Fix beam center
            """
            Reference: http://strucbio.biologie.uni-konstanz.de/xdswiki/index.php/Generate_XDS.INP
            """
            SNs1 = """\
449
472
474
912
923
933
911
446
916
""".split()
            SNs4 = """\
915
""".split()
            beamx, beamy = self.header["BeamX"], self.header["BeamY"]
            
            if self.header["SerialNumber"] in SNs1:
                self.header["BeamX"], self.header["BeamY"] = beamy, beamx
            elif self.header["SerialNumber"] in SNs4:
                pass
            else: # 3
                self.header["BeamY"] = self.header["Height"]*self.header["PixelY"] - beamy

        return self.header

    def export(self, exportType='xds'):
        """Try to interpret the header, and return a dictionary with the header
          information
          _FIXME_ This need to go to the Collect class
          """
        try:
            #print "Image. TRY import exporter:", exportType.lower()
            exporter = __import__(exportType.lower()+'_export')
        except:
            raise XIOError, "Can't load %s exporter" % (exportType)
        #
        exportDict = {}
        for k in exporter.HTD.keys():
            args, func = exporter.HTD[k]

            exportDict[k] = func(*map(self.header.get, args))
        return exportDict

    def info(self, verbose=0):
        "Print object internal information"
        #import pprint
        print ">> Interpreting header of image:  %s" % self.fileName
        print ">> Image format:      %s" % self.type
        if self.detModel:
            print ">> Detector type:     %s" % self.detModel
        if hasattr(self, "beamline") and self.beamline:
            print ">> Beamline guess:    %s %s %s" % self.beamline
        if self.intCompression:
            print ">> Image Compression:\t%s" % self.intCompression
        if verbose:
            print ">> Header: ",
            pprint.pprint(self.header)
        return self.header

    def info2(self):
        "Print object internal information"
        import pprint
        pprint.pprint(self.header)

    def data_info(self):
        "Print object internal information"
        #import pprint
        print ">> Loading %s" % (self.fileName)
        print ">> %d x %d pixels" % (self.header['Width'],
                                         self.header['Height'])
        print ">> Distance: %.1f mm, Lambda: %.3f A" % \
                           (self.header['Distance'],self.header['Wavelength'])
        try:
            data = self.getData()
            if self.type == 'marccd':
                data = filter(None, data)
                print data[:10]
            if self.type == 'raxis':
                hval = map(lambda val: (0x7fff&val)*8,
                              filter(lambda x: x>0x7fff, data))
                print hval
                print len(hval), max(hval)
            len_data = len(data)
            mean_data = add_reduce(data)/len_data
            print ">> MaxI: %d, AvgI: %.0f" % (max(data), mean_data)
        except XIOError:
            print ">> Don't know how (yet) to read %s compressed raw data." %\
                         self.intCompression

    def getData(self, clipping=False):
        """Read the image bytes. For now only support the 16bits unsigned,
        and uncompressed internaly. Can read compressed file directly
        (like .gz or .Z).
        If clipping=True, set I<0 to O and I>2**16 to 2**16"""

        if not self.interpreter:
            self.headerInterpreter()

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

            #print max(_data)
            #print min(_data)
            return _data

        elif self.intCompression == "CRYSALIS":
            _dataSize = self.header['Width']*self.header['Width']
            _image = self.open()

            # Jump over the header
            head = _image.read(self.header['HeaderSize'])
            # Read the remaining bytes
            OI = int(self.RawHeadDict["OI"])
            OL = int(self.RawHeadDict["OL"])

            endian = "<"
            _data = _image.read(_dataSize)
            _fmt = endian + "B" * _dataSize
            _data = struct.unpack(_fmt, _data)
            file_size = self.header['HeaderSize'] + _dataSize + OI*2 + OL*4
            print "Total file size = %d" % (file_size)
            _overloads_short = 0
            _overloads_long = 0
            if OI:
                #_fmt = self.header['EndianType'] + OI*"H"
                _fmt = endian + OI*"h"
                _overloads_short = _image.read(OI*2)
                _overloads_short = struct.unpack(_fmt, _overloads_short)
                print "short_max=", max(_overloads_short), 2**16
            if OL:
                #_fmt = self.header['EndianType'] + OL*"I"
                _fmt = endian + OL*"i"
                _overloads_long = _image.read(OL*4)
                _overloads_long = struct.unpack(_fmt, _overloads_long)
                print "long_max=", max(_overloads_long), 2**32
            _image.close()
            # unpack
            image = _dataSize*[0,]
            i, j, k = 1, 0, 0

            # first pixel
            if _data[0] == 254:
                image[0] = _overloads_short[j]
                j += 1
            elif _data[0] == 255:
                image[0] = _overloads_long[k]
                k += 1
            else:
                image[0] = _data[0] - 127

            # following pixels
            while i < _dataSize:
                tmp = _data[i]
                if tmp == 254:
                    image[i] = image[i-1] + _overloads_short[j]
                    j += 1
                elif tmp == 255:
                    image[i] = image[i-1] + _overloads_long[k]
                    k += 1
                else:
                    image[i] = image[i-1] + tmp - 127
                i+= 1
            i = 0
            print "MinI: %7d  MaxI: %7d  " % (min(image), max(image)),
            print "AvgI: %.1f" % (add_reduce(image)*1./len(image))
            if clipping:
                while i < _dataSize:
                    if image[i] < 0:
                        image[i] = 0
                    if image[i] > 2**16-1:
                        image[i] = 2**16 -1
                    i += 1
            #_data = [pix-127 for pix in _data if pix < 254]
            #print len(image) - _dataSize
            if _overloads_short:
                print "OI", OI, j#, _overloads_short[:20]
            if _overloads_long:
                print "OL", OL, k, _overloads_long[:20]
            #for i in range(30,500):
            #    print "%8d %8d %8d" % (i, _data[i], image[i])
            #for i in range(-100,-30):
            #    print "%8d %8d %8d" % (i, _data[i], image[i])
            return image

        else:
            raise XIOError, "Sorry, this image is internaly compressed."

    def make_dummy_pilatus_header(self):
        return """
# Detector: NOT PILATUS
# %(date)s
# Pixel_size %(pixel_size)e m x %(pixel_size)e m
# Exposure_time %(exp_time)f s
# Wavelength %(wavelength)f A
# Detector_distance %(distance)f m
# Beam_xy (%(orgx).2f, %(orgy).2f) pixels
# Start_angle %(osc_start)f deg.
# Angle_increment %(osc_width)f deg.
""" % dict(date=datetime.datetime.fromtimestamp(self.header["DateSeconds"]).strftime("%Y-%m-%dT%H:%M:%S.000"),
           pixel_size=self.header["PixelX"]*1e-3,
           exp_time=self.header["ExposureTime"],
           wavelength=self.header["Wavelength"],
           distance=self.header["Distance"]*1e-3,
           orgx=self.header["BeamX"]/self.header["PixelX"],
           orgy=self.header["BeamY"]/self.header["PixelY"],
           osc_start=self.header["PhiStart"],
           osc_width=self.header["PhiWidth"])


class Collect:
    """ A simple class to handle X-ray data collection informations.
    It can be contruct from either a string template or a
    tupple of 4 building elements (dir, prefix, nDigits, suffix).
    dir, prefix and suffix.

    >>> dc = Collect("/data/trp/ref-trypsin_1_001.img")
    >>> assert dc.directory == '/data/trp/'
    >>> assert dc.prefix == 'ref-trypsin_1'
    >>> assert dc.suffix == 'img'
    >>> assert dc.nDigits == 3

    >>> assert dc.mosflmTemplate == '/data/trp/ref-trypsin_1_###.img'
    >>> assert dc.xdsTemplate == '/data/trp/ref-trypsin_1_???.img'
    >>> assert dc.pythonTemplate == '/data/trp/ref-trypsin_1_%03d.img'
    >>> dc.setDirectory('./')
    >>> assert dc.mosflmTemplate == './ref-trypsin_1_###.img'

    The methode getNumber check the imageName and return the imageNumber.
    While checking imageNames, the dirPath and compression are optional.

    >>> dc = Collect("./test/ref-trypsin_1_001.img.gz")
    >>> assert dc.suffix == 'img'
    >>> assert dc.getNumber("/data/trp/ref-trypsin_1_023.img") == 23
    >>> assert dc.getNumber("./ref-trypsin_1_023.img") == 23
    >>> assert dc.getNumber("/data/trp/toto_1_023.img") == None

    The lookup methods are listing the filesystem directory and reporting
    the corresponding Collect image numbers, using the getNumber method.

    >>> assert dc.lookup_imageNumbers() == [1,91]
    >>> assert dc.lookup_imageRanges() == [[1, 1], [91, 91]]
    >>> dc.imageNumbers = [1,2,3,4,76,77,78,79,120,121]
    >>> assert dc.getClosestImage(91) == 79
    >>> mask =[[1,5],[61,65],[121,125]]
    >>> print  dc.lookup_imageRanges(forceCheck=False, mask_range=mask)
    [[1, 4], [121, 121]]

    >>> dc_match = dc.rec_imageNameID.match("/data/trp/ref-trypsin_1_001.img")
    >>> assert dc_match.groups() == ('/data/trp/', '001', None)
    """

    def __init__(self, init):

        _init_list_of_images = None
        self.imageNumbers = [] # set by lookup_imageNumbers()
        self.imageRanges = []  # set by lookup_imageRanges()
        self.image = None

        if type(init) == list and list_of_string(init):
            _init_list_of_images = init
            init = init[0]

        #print "D1: %s" % init
        if type(init) == str:
            self._naming_convension = 1
            M = REC_FULLIMAGENAME1.match(init)

            if not M:
                #print "2nd conv"
                M = REC_FULLIMAGENAME2.match(init)
                self._naming_convension = 2

            if not M:
                #print "3rd conv."
                M = REC_FULLIMAGENAME3.match(init)
                self._naming_convension = 3

            if not M:
                _err_message = "String constructor: %s doesn't match Collect "
                _err_message += "naming convention:\n"
                _err_message += "[dirPath/]+{prefix}_{num}.{suffix}[.compress]"
                raise XIOError, _err_message
                print _err_message

            #print M.groups()
            if self._naming_convension == 2:
                self.directory, self.prefix, num, compres = M.groups()
                self.suffix = ""
            else:
                self.directory, self.prefix, num, \
                                self.suffix, compres = M.groups()

            self.nDigits = len(num)
            # Just a verification that the regexp for the dir works...
            if self.directory:
                assert self.directory == os.path.split(init)[0] + "/"

        elif len(init) == 4 and type(init[0]) == str and type(init[1]) == str \
            and type(init[2]) == int and type(init[3]) == str:

            self.directory = init[0]
            self.prefix = init[1]
            self.nDigits = int(init[2])
            self.suffix = init[3]
        else:
            raise TypeError, "Unexpected Collect contructor argument:", init

        if not self.directory:
            self.directory = "./"

        self.initialImageName = init
        self.setTemplates()

        # self.rec_imageNameID is an identifier for the instanced Collect:
        # Verification is made that the prefix, suffix and number of digits
        # matches. dirPath and extCompression are optional.

        if self._naming_convension == 1:  re_fmt = r"%s_([0-9]{%d,})\.%s"
        elif self._naming_convension == 2: re_fmt = r"%s.([0-9]{%d,})%s"
        elif self._naming_convension == 3: re_fmt = r"%s([0-9]{%d,})\.%s"

        _re_imageNameID = re_fmt % (self.prefix, self.nDigits, self.suffix)
        self.rec_imageNameID = re.compile(RE_DIRPATH +
                                          _re_imageNameID +
                                          RE_EXT_COMPRESSION)
        if _init_list_of_images:
            self.lookup_imageNumbers(_init_list_of_images)

    def interpretImage(self):
        """Try to find and intepret the first image of the collect.
           Store the intepreted Image object as self.image
        """
        if not self.imageNumbers:
            self.lookup_imageNumbers()
        #print self.initialImageName
        self.image = Image(self.initialImageName)
        self.imageType = self.image.type

    def setTemplates(self):
        """Set template format for xds, mosflm and python. Used by constructor,
        and updated by setDirPath()"""

        if self._naming_convension == 1:  sep = "_%s."
        elif self._naming_convension == 3: sep = "%s."
        elif self._naming_convension == 2: sep = ".%s"

        # generic template format without external compression
        self.formatTemplate = self.directory + \
                              self.prefix + sep + \
                              self.suffix

        # python template format
        self.pythonTemplate = self.formatTemplate % \
                                     ("%0" + "%sd" % self.nDigits)

        # MOSFLM Template
        self.mosflmTemplate = self.formatTemplate % ("#" * self.nDigits)

        # XDS Template
        self.xdsTemplate = self.formatTemplate % ("?" * self.nDigits)
        if len(self.xdsTemplate) > 50 and VERBOSE:
            _mess = " Warning NAME_TEMPLATE_OF_DATA_FRAMES has more than 50 ",
            _mess += "characters! XDS will stop. Lengthy path names should "
            _mess += "be abbreviated by a symboliclink for frames directory."
            print _mess

    def setDirectory(self, newpath):
        "Change the absolute dirPath information and update templates."

        if type(newpath) == str or type(newpath) == unicode:
            if newpath[-1] != "/": newpath += "/"
            self.directory = newpath
            self.setTemplates()
        else:
            raise TypeError, "Unexpected argument type in setDirectory:", \
                                                        type(newpath)


    def getNumber(self, imageName):
        "If imageName matchs: return the int(number)."
        m = self.rec_imageNameID.match(imageName)
        if m: return int(m.group(2))

    def lookup_imageNumbers(self, files=None):
        "Return a list of matching image number. Removes duplicate numbers."

        if not files:
            files = os.listdir(self.directory)
        files.sort()

        # Get numbers and remove non matching
        images_num = filter(None, map(self.getNumber, files))
        #print len(images_num), ( images_num[-1] - images_num[0] + 1)

        # If the list contains duplicates: remove duplicate numbers.
        if len(images_num):
            if len(images_num) > ( images_num[-1] - images_num[0] + 1):
                remove_redundancy(images_num)
            self.imageNumbers = images_num
            self.first_imageNumber = images_num[0]
            self.last_imageNumber = images_num[-1]
        else:
            self.imageNumbers = []
            self.first_imageNumber = 0
            self.last_imageNumber = 0

        return images_num


    def _sequence_to_ranges(self, seq):
        _range = []
        prev, i = seq[0]-1, seq[0]
        for n in seq:
            if n != prev+1:
                _range.append([i, prev])
                i = n
            prev  = n
        _range.append([i, prev])
        return _range

    def _ranges_to_sequence(self, ranges):
        _rangePlus1 = lambda R: range(R[0], R[1]+1)
        _add = lambda a, b: a+b
        return reduce(_add, map(_rangePlus1, ranges))


    def lookup_imageRanges(self, forceCheck=False, mask_range=None):
        """Return a range list of consecutive image number.
        For example: [[1,3],[91,93]] or [[1,90]]
        If mask_range: Map a model ranges to the existing image range.
        Example of mask: [[1,10],[41,50],[81,90]]
        """

        # Make sure we have a sequence of imageNumbers to work with
        if not self.imageNumbers or forceCheck:
            self.lookup_imageNumbers()

        # Transform the sequence to ranges
        self.imageRanges = self._sequence_to_ranges(self.imageNumbers)

        # Map the mask if we have one, or return the ranges
        if mask_range:
            mask_seq = self._ranges_to_sequence(mask_range)
            intersect = [m for m in mask_seq if m in self.imageNumbers]
            return self._sequence_to_ranges(intersect)
        else:
            return self.imageRanges

    def isContinuous(self, imagename_list, methode=0, _epsilon=1.5e-1):
        """Return true if the collect is supposed to be a serie of images
           with consecutive phi angles."""
        if not self.imageNumbers:
            self.lookup_imageNumbers()
        if not self.image:
            self.image = Image(self.initialImageName)
        last_image_name = self.pythonTemplate % self.imageNumbers[-1]
        if os.path.exists(last_image_name):
            last_image = last_image_name
        elif os.path.exists(last_image_name+".gz"):
            last_image = last_image_name + ".gz"
        elif os.path.exists(last_image_name+".bz2"):
            last_image = last_image_name + ".bz2"
        # Check only for the initialImageName and last images of the list
        # their phi_osc range is compatible with their numbering.
        if methode == 0:
            phi_range = self.image.header['PhiWidth']
            if phi_range == 0:
                print  "WARNING: Oscillation_Range recorded in image header is null ! Use the -O option to give the true value."
                return True
            phi_start = self.image.header['PhiStart']
            phi_last = Image(last_image).header['PhiStart']
            num_start = self.imageNumbers[0]
            num_last = self.imageNumbers[-1]
            diff = (phi_last - phi_start)/phi_range - (num_last - num_start)
            #print phi_last, phi_start, phi_range, diff
        if diff < _epsilon:
            return True
        else:
            fmt = "Image: %s  phi_start= %8.3f  phi_range= %8.3f"
            print fmt % (imagename_list[0], phi_start, phi_range)
            print fmt % (imagename_list[-1], phi_last, phi_range)
            print "Diff = (phi_last - phi_start)/phi_range -",
            print "(num_last - num_start) = %.4f" % diff
            return False

    def get_range(self, minf=None, maxf=None):
        if self.imageRanges:
            min_c, max_c = self.imageRanges[0][0], self.imageRanges[-1][-1]
            if minf: min_c = max(minf, min_c)
            if maxf: max_c = min(maxf, max_c)
            return [min_c, max_c]
        else:
            return []

    def getClosestImage(self, target, forceCheck=False):
        """Return closest existing image number to target. Quickly."""

        if not self.imageNumbers or forceCheck:
            self.lookup_imageNumbers()

        substract = lambda x: abs(x-int(target))
        diff = map(substract, self.imageNumbers)
        return self.imageNumbers[diff.index(min(diff))]

    def export(self, exportType='xds'):
        "Try to interpret the collect, and return an interpreted dictionary"
        try:
            #print "Collect. TRY import exporter:", exportType.lower()
            exporter = __import__(exportType.lower()+'_export')
        except:
            raise XIOError, "Can't load %s exporter" % (exportType)

        if not self.image.header:
            print "Not sefl.image.header"
            self.interpretImage()

        self.lookup_imageRanges()
        try:
            exportDict = self.image.export(exportType)
        #except ValueError:
        except NameError:
            exportDict = {}
        for k in exporter.CTD.keys():
            args, func = exporter.CTD[k]
            exportDict[k] = func(*map(self.__dict__.get, args))
        exportDict["SPECIFIC_KEYWORDS"] = ""
        det_SN = self.image.header["SerialNumber"]
        try:
            spec_SN = exporter.SPECIFIC_SUPPLEMENTARY_KEYWORDS
            for spec_type in spec_SN.keys():
                if spec_type in det_SN:
                    exportDict["SPECIFIC_KEYWORDS"] = spec_SN[spec_type]
        except:
            pass
        return exportDict

    def export_template(self, exportType='xds'):
        try:
            exporter = __import__(exportType.lower()+'_export')
        except:
            raise XIOError, "Can't load %s exporter" % (exportType)
        exportDict = self.export(exportType)
        #import pprint
        #pprint.pprint(exportDict)
        return exporter.TEMPLATE % exportDict

    def get_export_template(self, exportType='xds'):
        try:
            exporter = __import__(exportType.lower()+'_export')
        except:
            raise XIOError, "Can't load %s exporter" % (exportType)
        return exporter.TEMPLATE

class Interpreter:
    """ A basic interpreter class.
    The interpreter mecanism is based on Translation Dictionaries (TD).
    These TD are containers of the rules used to translate from one dictionary to
    another... Am I clear enought?
    """

    def __init__(self, EndianType='<', ):
        self.EndianType = EndianType

    def interpretHeader(rawData):
        """ Return a header dict """
        # It seems that
        EndianType = ">"
        headerStructureFmt = EndianType + \
                             "".join([fmt for key,fmt in headerStructure])
        headerStructureKeys = [key for key, fmt in headerStructure]

        # uncoding from headerStructureFmt
        read_size = struct.calcsize(headerStructureFmt)
        read_unp = struct.unpack(headerStructureFmt, rawHead[:read_size])

        RawHeadDict = dict(zip(headerStructureKeys, read_unp))
        RawHeadDict.update({'MESSAGE':'', 'EndianType':EndianType})

        return header

    def interpretData(rawData):
        """ Return a data list """
        return data


def test_regexp():
    "Just testing some cases of image names"
    r_1 = REC_FULLIMAGENAME1.match("./ref-trypsin_1_023.img.gz")
    r_2 = REC_FULLIMAGENAME1.match("ref-trypsin_1_023.img")
    r_3 = REC_FULLIMAGENAME1.match("_1_023.img")
    r_4 = REC_FULLIMAGENAME1.match("/d/r/_9923.i.bz2")
    r_5 = REC_FULLIMAGENAME3.match("/d/r/toto9923.i.bz2")
    r_6 = REC_FULLIMAGENAME3.match("/d/r/tot1o9923.i.bz2")
    r_7 = REC_FULLIMAGENAME3.match("/d/r/9923.i.bz2")
    r_8 = REC_FULLIMAGENAME2.match("/d/r/toto.9923")
    r_9 = REC_FULLIMAGENAME2.match("/d/r/toto.9923.bz2")
    r_a = REC_FULLIMAGENAME2.match("/d/r/1.9923.bz2")
    r_b = REC_FULLIMAGENAME2.match("/d/r/654082.0086")
    r_c = REC_FULLIMAGENAME3.match("/TH050916A0001.img")
    print r_1.groups()
    print r_2.groups()
    print r_3.groups()
    print r_4.groups()
    print r_8.groups()
    print r_9.groups()
    print r_a.groups()
    print r_b.groups()
    print r_c.groups()
    assert r_1.groups() == ('./', 'ref-trypsin_1', '023', 'img', '.gz')
    assert r_2.groups() == (None, 'ref-trypsin_1', '023', 'img', None)
    assert r_3.groups() == (None, '_1', '023', 'img', None)
    assert r_4.groups() == ('/d/r/', '', '9923', 'i', '.bz2')
    #print r5.groups()
    #if r6: print r6.groups()
    #if r7: print r7.groups()
    assert r_5.groups() == ('/d/r/', 'toto9', '923', 'i', '.bz2')
    assert r_6.groups() == ('/d/r/', 'tot1o9', '923', 'i', '.bz2')
    assert r_7.groups() == ('/d/r/', '9', '923', 'i', '.bz2')

def test0():
    "Doctest."
    test_regexp()
    import doctest
    doctest.testmod()

def test1(filename):
    "Simple test"
    im = Image(filename)
    im.info()
    dc = Collect(filename)

def test2(filename):
    "More complete test"
    import time
    from DiffractionImage_getHeader import GetHeader
    s1 = time.time()

    im = Image(filename)
    DH1 = im.headerInterpreter()
    s2 = time.time()
    print s2-s1
    DH2 = GetHeader(filename)
    s3 = time.time()
    print s3-s2
    #pprint.pprint(DH2)
    for k in DH2.keys():
        if not DH1[k] == DH2[k]:
            try:
                _diff = abs(DH1[k] - DH2[k])
                assert _diff/DH1[k] < 1.e-7
            except  AssertionError:
                print "Error. Difference %s between value '%s': %s != %s" % \
                                 (_diff/DH1[k],k,DH1[k],DH2[k])
    pprint.pprint(DH1)


def test3(filename):
    im = Image(filename)
    im.info()
    pprint.pprint(im.export('xds'))

if __name__ == "__main__":
    import pprint
    test0()
    #test_regexp()
    if len(sys.argv) > 1:
        filename = sys.argv[1]
        test1(filename)
        test3(filename)
