"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

import os
import json
from yamtbx.dataproc import XIO
from yamtbx.dataproc import cbf
from yamtbx.dataproc.dataset import group_img_files_template
from yamtbx.dataproc.xds import get_xdsinp_keyword
import cStringIO

def sensor_thickness_from_minicbf(img):
    header = cbf.get_pilatus_header(img)
    sensor = filter(lambda x: "sensor, thick" in x, header.splitlines())[0]
    thick, unit = sensor[sensor.index("thickness ")+len("thickness "):].split()
    assert unit == "m"

    return float(thick) * 1.e3
# sensor_thickness_from_minicbf()

def import_geometry(xds_inp=None, dials_json=None):
    assert (xds_inp, dials_json).count(None) == 1

    geom_kwds = set(["DIRECTION_OF_DETECTOR_X-AXIS", "DIRECTION_OF_DETECTOR_Y-AXIS",
                     "DETECTOR_DISTANCE", "ORGX", "ORGY", "ROTATION_AXIS", # "X-RAY_WAVELENGTH",
                     "INCIDENT_BEAM_DIRECTION", "SEGMENT", "DIRECTION_OF_SEGMENT_X-AXIS",
                     "DIRECTION_OF_SEGMENT_Y-AXIS", "SEGMENT_DISTANCE",
                     "SEGMENT_ORGX", "SEGMENT_ORGY"])

    # FIXME in case of multi-segment detector..

    if xds_inp:
        inp = get_xdsinp_keyword(xds_inp)
        inp = filter(lambda x: x[0] in geom_kwds, inp)
        return map(lambda x: "%s= %s"%x, inp)
    elif dials_json:
        import dxtbx.imageset
        from dxtbx.serialize.load import _decode_dict
        from dxtbx.model import BeamFactory
        from dxtbx.model import DetectorFactory
        from dxtbx.model import GoniometerFactory
        from dxtbx.model import ScanFactory
        from dxtbx.serialize.xds import to_xds
        j = json.loads(open(dials_json).read(), object_hook=_decode_dict)
        # dummy
        sweep = dxtbx.imageset.ImageSetFactory.from_template("####",
                                                             image_range=[1,1],check_format=False)[0]
        sweep.set_detector(DetectorFactory.from_dict(j["detector"][0]))
        sweep.set_beam(BeamFactory.from_dict(j["beam"][0]))
        sweep.set_goniometer(GoniometerFactory.from_dict(j["goniometer"][0]))
        sweep.set_scan(ScanFactory.make_scan(image_range=[1,1], exposure_times=[1], oscillation=[1,2], epochs=[0])) # dummy
        inp = get_xdsinp_keyword(inp_str=to_xds(sweep).XDS_INP())
        inp = filter(lambda x: x[0] in geom_kwds, inp)
        return map(lambda x: "%s= %s"%x, inp)

    return []
# import_geometry()

def read_geometry_using_dxtbx(img_file):
    import dxtbx.datablock
    import dxtbx.serialize.xds

    geom_kwds = set(["DIRECTION_OF_DETECTOR_X-AXIS", "DIRECTION_OF_DETECTOR_Y-AXIS",
                     "DETECTOR_DISTANCE", "ORGX", "ORGY", "ROTATION_AXIS", "X-RAY_WAVELENGTH",
                     "DETECTOR", "MINIMUM_VALID_PIXEL_VALUE", "OVERLOAD", "SENSOR_THICKNESS",
                     "NX", "NY", "QX", "QY", "STARTING_ANGLE", "OSCILLATION_RANGE",
                     "FRACTION_OF_POLARIZATION", "POLARIZATION_PLANE_NORMAL", 
                     "INCIDENT_BEAM_DIRECTION", "SEGMENT", "DIRECTION_OF_SEGMENT_X-AXIS",
                     "DIRECTION_OF_SEGMENT_Y-AXIS", "SEGMENT_DISTANCE",
                     "SEGMENT_ORGX", "SEGMENT_ORGY"])

    datablocks = dxtbx.datablock.DataBlockFactory.from_filenames([img_file])
    to_xds = dxtbx.serialize.xds.to_xds(datablocks[0].extract_sweeps()[0])
    inp = get_xdsinp_keyword(inp_str=to_xds.XDS_INP())
    inp = filter(lambda x: x[0] in geom_kwds, inp)
    return to_xds, map(lambda x: " %s= %s"%x, inp)

# read_geometry_using_dxtbx()

def generate_xds_inp(img_files, inp_dir, use_dxtbx=False, anomalous=True,
                     reverse_phi=None, spot_range=None, minimum=False,
                     crystal_symmetry=None, integrate_nimages=None,
                     osc_range=None, orgx=None, orgy=None, rotation_axis=None, distance=None,
                     wavelength=None,
                     minpk=None, exclude_resolution_range=None,
                     fstart=None, fend=None, extra_kwds=None, overrides=None, fix_geometry_when_overridden=False):
    """
    Reference: http://strucbio.biologie.uni-konstanz.de/xdswiki/index.php/Generate_XDS.INP
    """
    if not exclude_resolution_range: exclude_resolution_range = []
    if not extra_kwds: extra_kwds = []
    if not overrides: overrides = []

    is_eiger_hdf5 = (len(img_files) == 1 and "_master.h5" in img_files[0])

    if is_eiger_hdf5:
        template = img_files[0].replace("_master.h5","_??????.h5")
    else:
        groups = group_img_files_template(img_files)
        assert len(groups) == 1
        template, fstart, fend = groups[0] # arguments fstart, fend are overwritten here
    #print inp_dir, img_files[0], template

    tmp = [os.path.dirname(os.path.relpath(img_files[0], inp_dir)),
           os.path.dirname(os.path.abspath(img_files[0]))]
    imdir = min(tmp, key=lambda x:len(x))
    template = os.path.join(imdir, os.path.basename(template))
    #print imdir

    if crystal_symmetry is None:
        sgnum = 0
        cell_str = "50 60 70 90 90 90"
    else:
        sgnum = crystal_symmetry.space_group_info().type().number()
        cell_str = " ".join(map(lambda x: "%.2f"%x, crystal_symmetry.unit_cell().parameters()))

    data_range = "%d %d" % (fstart, fend)
    if spot_range is None:
        spot_range = "%d %d" % (fstart, (fend+fstart)/2)
    elif spot_range == "first":
        spot_range = "%d %d" % (fstart, fstart)
    elif spot_range == "all":
        spot_range = "%d %d" % (fstart, fend)
    elif len(spot_range) == 2:
        spot_range = "%d %d" % spot_range
    else:
        print "Error!"
        return

    friedel = "FALSE" if anomalous else "TRUE"
    is_pilatus_or_eiger = False

    img_files_existed = filter(lambda x: os.path.isfile(x), img_files)
    if not img_files_existed: raise Exception("No actual images found.")

    inp_str = """\
 MAXIMUM_NUMBER_OF_JOBS= 1
 JOB= XYCORR INIT COLSPOT IDXREF DEFPIX INTEGRATE CORRECT
 NAME_TEMPLATE_OF_DATA_FRAMES= %(template)s
 DATA_RANGE= %(data_range)s
 SPOT_RANGE= %(spot_range)s
!BACKGROUND_RANGE=1 10
 FRIEDEL'S_LAW= %(friedel)s
""" % locals()

    if use_dxtbx:
        # This is the current limitation..
        assert all(x is None for x in (osc_range, orgx, orgy, rotation_axis, distance, wavelength))

        toxds, inplst = read_geometry_using_dxtbx(img_files_existed[0])
        inp_str += "\n".join(inplst) + "\n"
        osc_range = toxds.oscillation_range
        nx, ny = toxds.detector_size
        is_pilatus_or_eiger = toxds.get_detector()[0].get_type() == "SENSOR_PAD"

    else:
        im = XIO.Image(img_files_existed[0])

        if osc_range is None: osc_range = im.header["PhiWidth"]

        if rotation_axis is None: # automatic decision
            if "OscAxisVec" in im.header:
                rotation_axis = im.header["OscAxisVec"]
                print "DEBUG::rotation_axis from header:", rotation_axis
            else:
                if im.header["ImageType"] == "raxis": rotation_axis = (0,1,0)
                else: rotation_axis = (1,0,0)

                if reverse_phi is None: # automatic decision
                    REVERSEPHI_SNs=dict(marccd="""\
24
31
38
40
42
106
""".split(), # Known detectors for reversed-phi in SPring-8: 24: BL26B2 Mar225, 31: BL32XU MX225HE, 38: BL44XU MX225HE, 42: BL44XU MX300HE, 40: BL41XU MX225HE, 106: BL32XU MX225HS
                                        adsc="""\
915
458
924
""".split(), # Known detectors for reversed-phi in SPring-8: 915: BL38B1 Q315; APS 19-ID: 458; BM30A: 924
                                        minicbf="""\
PILATUS3 6M, S/N 60-0125
PILATUS3 6M, S/N 60-0132
PILATUS 2M, S/N 24-0109
""".splitlines(), # Known detectors for reversed-phi in SPring-8: BL41XU PILATUS3 6M 60-0125, APS: 19ID PILATUS3 6M 60-0132, MX2 beamline (Brazilian Synchrotron National Laboratory - LNLS)
                                        )
                    if im.header.get("SerialNumber") in REVERSEPHI_SNs.get(im.header["ImageType"], ()):
                        print "DEBUG:: this is reversephi of", rotation_axis
                        reverse_phi = True

                if reverse_phi:
                    rotation_axis = map(lambda x:-1*x, rotation_axis)

        rotation_axis = " ".join(map(lambda x: "%.6f"%x, rotation_axis))

        nx, ny = im.header["Width"], im.header["Height"],
        qx, qy = im.header["PixelX"], im.header["PixelY"]
        if orgx is None: orgx = im.header["BeamX"]/qx
        if orgy is None: orgy = im.header["BeamY"]/qy
        if wavelength is None: wavelength = im.header["Wavelength"]
        if distance is None: distance = im.header["Distance"]

        sensor_thickness = 0 # FIXME

        if im.header["ImageType"] == "marccd":
            detector = "CCDCHESS MINIMUM_VALID_PIXEL_VALUE= 1 OVERLOAD= 65500"
        elif im.header["ImageType"] == "raxis":
            detector = "RAXIS MINIMUM_VALID_PIXEL_VALUE= 0  OVERLOAD= 2000000"
            distance *= -1
        elif im.header["ImageType"] == "minicbf":
            detector = "PILATUS MINIMUM_VALID_PIXEL_VALUE=0 OVERLOAD= 1048576"
            sensor_thickness = sensor_thickness_from_minicbf(img_files[0])
            is_pilatus_or_eiger = True
        elif im.header["ImageType"] == "adsc":
            detector = "ADSC MINIMUM_VALID_PIXEL_VALUE= 1 OVERLOAD= 65000"
        elif im.header["ImageType"] == "mscccd":
            detector = "SATURN MINIMUM_VALID_PIXEL_VALUE= 1 OVERLOAD= 262112" # XXX Should read header!!
            distance *= -1
        elif is_eiger_hdf5:
            detector = "EIGER MINIMUM_VALID_PIXEL_VALUE=0 OVERLOAD= %d" % im.header["Overload"]
            sensor_thickness = im.header["SensorThickness"]
            is_pilatus_or_eiger = True

        inp_str += """\
 ORGX= %(orgx).2f ORGY= %(orgy).2f
 DETECTOR_DISTANCE= %(distance).2f
 OSCILLATION_RANGE= %(osc_range).3f
 X-RAY_WAVELENGTH= %(wavelength).5f
 DETECTOR= %(detector)s
 SENSOR_THICKNESS= %(sensor_thickness).2f
 NX= %(nx)s NY= %(ny)s  QX= %(qx)s  QY= %(qy)s
 ROTATION_AXIS= %(rotation_axis)s
 INCIDENT_BEAM_DIRECTION= 0 0 1
 FRACTION_OF_POLARIZATION= 0.98
 POLARIZATION_PLANE_NORMAL= 0 1 0
""" % locals()

        # XXX Synchrotron can have R-AXIS, and In-house detecotr can have horizontal goniometer!
        if im.header["ImageType"] == "raxis":
            inp_str += """\
 DIRECTION_OF_DETECTOR_X-AXIS= 1 0 0
 DIRECTION_OF_DETECTOR_Y-AXIS= 0 -1 0
 INCIDENT_BEAM_DIRECTION= 0 0 1
!FRACTION_OF_POLARIZATION= 0.98   ! uncomment if synchrotron
 POLARIZATION_PLANE_NORMAL= 1 0 0
"""
        else:
            if im.header["ImageType"] == "mscccd":
                inp_str += """\
 DIRECTION_OF_DETECTOR_X-AXIS= -1 0 0
 DIRECTION_OF_DETECTOR_Y-AXIS=  0 1 0
"""
            else:
                inp_str += """\
 DIRECTION_OF_DETECTOR_X-AXIS= 1 0 0
 DIRECTION_OF_DETECTOR_Y-AXIS= 0 1 0
"""

    if integrate_nimages is None:
        extra_kwds.append(" DELPHI= 5")
    else:
        extra_kwds.append(" DELPHI= %.4f" % osc_range * integrate_nimages)

    if minpk is not None: extra_kwds.append(" MINPK= %.2f" % minpk)
    for r1, r2 in exclude_resolution_range:
        if r1 < r2: r1, r2 = r2, r1
        extra_kwds.append(" EXCLUDE_RESOLUTION_RANGE= %.3f %.3f" % (r1, r2))

    extra_kwds = "\n".join(extra_kwds) + "\n"
    inp_str += extra_kwds

    if not minimum:
        inp_str += """\
 SPACE_GROUP_NUMBER= %(sgnum)d
 UNIT_CELL_CONSTANTS= %(cell)s
 INCLUDE_RESOLUTION_RANGE=50 0

 TRUSTED_REGION=0.00 1.4
 VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS=6000. 30000.
 STRONG_PIXEL=4
 MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT=3
 REFINE(IDXREF)=CELL BEAM ORIENTATION AXIS ! DISTANCE POSITION
 REFINE(INTEGRATE)= DISTANCE POSITION BEAM ORIENTATION ! AXIS CELL
!REFINE(CORRECT)=CELL BEAM ORIENTATION AXIS DISTANCE POSITION
""" % dict(sgnum=sgnum, cell=cell_str)

    if is_pilatus_or_eiger:
        inp_str += """\
 SEPMIN=4 CLUSTER_RADIUS=2
 NUMBER_OF_PROFILE_GRID_POINTS_ALONG_ALPHA/BETA= 13
 NUMBER_OF_PROFILE_GRID_POINTS_ALONG_GAMMA= 13
"""
        if nx == 1475 and ny == 1679: # Pilatus 2M
            inp_str += """\
!EXCLUSION OF VERTICAL DEAD AREAS OF THE PILATUS 2M DETECTOR
 UNTRUSTED_RECTANGLE= 486  496     0 1680
 UNTRUSTED_RECTANGLE= 980  990     0 1680
!EXCLUSION OF HORIZONTAL DEAD AREAS OF THE PILATUS 2M DETECTOR
 UNTRUSTED_RECTANGLE=   0 1476   194  214
 UNTRUSTED_RECTANGLE=   0 1476   406  426
 UNTRUSTED_RECTANGLE=   0 1476   618  638
 UNTRUSTED_RECTANGLE=   0 1476   830  850
 UNTRUSTED_RECTANGLE=   0 1476  1042 1062
 UNTRUSTED_RECTANGLE=   0 1476  1254 1274
 UNTRUSTED_RECTANGLE=   0 1476  1466 1486
"""
        elif nx == 2463 and ny == 2527: # Pilatus 6M
            inp_str += """\
 UNTRUSTED_RECTANGLE= 487  495     0 2528
 UNTRUSTED_RECTANGLE= 981  989     0 2528
 UNTRUSTED_RECTANGLE=1475 1483     0 2528
 UNTRUSTED_RECTANGLE=1969 1977     0 2528
 UNTRUSTED_RECTANGLE=   0 2464   195  213
 UNTRUSTED_RECTANGLE=   0 2464   407  425
 UNTRUSTED_RECTANGLE=   0 2464   619  637
 UNTRUSTED_RECTANGLE=   0 2464   831  849
 UNTRUSTED_RECTANGLE=   0 2464  1043 1061
 UNTRUSTED_RECTANGLE=   0 2464  1255 1273
 UNTRUSTED_RECTANGLE=   0 2464  1467 1485
 UNTRUSTED_RECTANGLE=   0 2464  1679 1697
 UNTRUSTED_RECTANGLE=   0 2464  1891 1909
 UNTRUSTED_RECTANGLE=   0 2464  2103 2121
 UNTRUSTED_RECTANGLE=   0 2464  2315 2333
"""
        elif nx == 3110 and ny == 3269: # Eiger 9M
            inp_str += """\
 UNTRUSTED_RECTANGLE= 1029 1042 0 3269
 UNTRUSTED_RECTANGLE= 2069 2082 0 3269
 UNTRUSTED_RECTANGLE= 0 3110  513  553
 UNTRUSTED_RECTANGLE= 0 3110 1064 1104
 UNTRUSTED_RECTANGLE= 0 3110 1615 1655
 UNTRUSTED_RECTANGLE= 0 3110 2166 2206
 UNTRUSTED_RECTANGLE= 0 3110 2717 2757
"""

    if overrides:
        inp_str += "\n! Overriding parameters:\n"
        inp_str += "\n".join(overrides)+"\n"

        if fix_geometry_when_overridden:
            inp_str += " REFINE(IDXREF)= CELL ORIENTATION ! BEAM AXIS DISTANCE POSITION\n"
            inp_str += " REFINE(INTEGRATE)= CELL ORIENTATION ! DISTANCE POSITION BEAM AXIS\n"

    return inp_str
