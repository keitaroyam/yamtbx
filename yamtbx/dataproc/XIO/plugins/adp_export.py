from __future__ import unicode_literals

_to_adp_beamlines = {
"mar":       "MAR345",
"marccd":    "MARCCD",
"adsc":      "ADSC",
"raxis":     "RAXIS",
}

#     Header Translator Dictionary.
#     Translate image header entries in a new dictionay
#     newdic['wavelength'] = float(head['Wavelength'])
#
HTD = {
'wavelength':(['Wavelength'], float),
'distance':(['Distance'], float),
'phi_init':(['PhiStart'], float),
'delta_phi':(['PhiWidth'], float),
'beam_y':(['BeamY'], float),
'beam_x':(['BeamX'], float), 
'angle_det':(['TwoTheta'], float),
'resolution_range':(['EdgeResolution'], lambda x: "30. %.2f" % (x))
}
        

#     Collect Translator Dictionary.
#     Translate collect object attributes to a new dictionay
#     newdic['suffix'] = str(collect.suffix)
#
CTD = {
'img_dir':(['directory'], str),
'prefix':(['prefix'], lambda x: str(x)+"_"),
'suffix':(['suffix'], str),
'frame_first':(['first_imageNumber'], int),
'frame_last':(['last_imageNumber'], int),
'beamline':(['imageType'], lambda x: _to_adp_beamlines[x])
}

        
