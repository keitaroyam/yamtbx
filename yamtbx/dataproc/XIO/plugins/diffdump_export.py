# -*- coding: utf-8 -*-

""" XIO plugin for the export parameters to mimic diffdump of the
the DiffractionImage package.
"""

__version__ = "0.1.1"
__author__ = "Pierre Legrand (pierre.legrand@synchrotron-soleil.fr)"
__date__ = "1-12-2009"
__copyright__ = "Copyright (c) 2009 Pierre Legrand"
__license__ = "New BSD, http://www.opensource.org/licenses/bsd-license.php"


DETECTOR_DICT = {
  "format":{
    "mar":       "MAR345",
    "mar555":    "MAR345",
    "marccd":    "MARCCD",
    "adsc":      "SMV",
    "raxis":     "R-AXIS",
    "minicbf":   "MINI-CBF",
    "mscccd":    "SMV",
    "oxford":    "CRYSALIS"
  },
  "manufacturer":{
    "mar":       "MAR",
    "mar555":    "MAR",
    "marccd":    "MAR",
    "adsc":      "ADSC",
    "raxis":     "RIGAKU",
    "minicbf":   "DECTRIS",
    "mscccd":    "RIGAKU",
    "oxford":    "OXFORD-DIFFRACTION",
  },
}

TEMPLATE = """Format : %(format)s
Manufacturer : %(manufacturer)s
Collection date : %(date)s
Exposure time : %(exposure_time).6f s
Detector S/N : %(serial)s
Wavelength : %(wavelength).6f Ang
Beam center : (%(beamx).6f mm,%(beamy).6f mm)
Distance to detector : %(distance).6f mm
Image Size : (%(width)d px, %(height)d px)
Pixel Size : (%(pixelx).6f mm, %(pixely).6f mm)
Oscillation (%(osc_axis)s) : %(osc_start).6f -> %(osc_end).6f deg
Two Theta value: %(twotheta).6f deg"""


#     Header Translator Dictionary.
#     Translate image header entries in a new dictionay
#     newdic['X_RAY_WAVELENGTH'] = float(head['Wavelength'])
#
HTD = {
'serial': (['SerialNumber'], str),
'date': (['DateStr'], str),
'exposure_time': (['ExposureTime'], float),
'wavelength':(['Wavelength'], float),
'distance':(['Distance'], float),
'twotheta':(['TwoTheta'], float),
'osc_start':(['PhiStart'], float),
'osc_end':(['PhiStart', 'PhiWidth'], lambda x,y: x + y),
'osc_range':(['PhiWidth'],float),
'osc_axis':(['OscAxis'], str),
'width':(['Width'], int),
'height':(['Height'], int),
'pixelx':(['PixelX'], float),
'pixely':(['PixelY'], float),
'beamx':(['BeamX'], float),
'beamy':(['BeamY'], float),
}

#     Collect Translator Dictionary.
#     Translate collect object attributes to a new dictionay
#     newdic['SPOT_RANGE'] = list(collect.imageRanges)
#
CTD = {
'format': (['imageType'], lambda x: DETECTOR_DICT["format"][x]),
'manufacturer': (['imageType'], lambda x: DETECTOR_DICT["manufacturer"][x]),
}
