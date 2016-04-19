"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

import os
from yamtbx.dataproc import XIO
from yamtbx.dataproc import cbf
from yamtbx.dataproc.dataset import group_img_files_template

def sensor_thickness_from_minicbf(img):
    header = cbf.get_pilatus_header(img)
    sensor = filter(lambda x: "sensor, thick" in x, header.splitlines())[0]
    thick, unit = sensor[sensor.index("thickness ")+len("thickness "):].split()
    assert unit == "m"

    return float(thick) * 1.e3
# sensor_thickness_from_minicbf()

# Reference: http://strucbio.biologie.uni-konstanz.de/xdswiki/index.php/Generate_XDS.INP
def generate_xds_inp(img_files, inp_dir, reverse_phi, anomalous, spot_range=None, minimum=False,
                     crystal_symmetry=None, integrate_nimages=None,
                     osc_range=None, orgx=None, orgy=None, rotation_axis=None, distance=None,
                     wavelength=None,
                     minpk=None, exclude_resolution_range=[],
                     fstart=None, fend=None, extra_kwds=[]):
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

    im = None
    for imgfile in img_files:
        if os.path.isfile(imgfile):
            im = XIO.Image(imgfile)
            break
    if im is None:
        raise Exception("No actual images found.")

    if crystal_symmetry is None:
        sgnum = 0
        cell_str = "50 60 70 90 90 90"
    else:
        sgnum = crystal_symmetry.space_group_info().type().number()
        cell_str = " ".join(map(lambda x: "%.2f"%x, crystal_symmetry.unit_cell().parameters()))

    if osc_range is None: osc_range = im.header["PhiWidth"]

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

    if rotation_axis is None:
        if im.header["ImageType"] == "raxis": rotation_axis = (0,1,0)
        else: rotation_axis = (1,0,0)

    if reverse_phi: rotation_axis = map(lambda x:-1*x, rotation_axis)
    rotation_axis = " ".join(map(lambda x: "%.2f"%x, rotation_axis))

    if integrate_nimages is None:
        delphi = 5
    else:
        delphi = osc_range * integrate_nimages

    nx, ny = im.header["Width"], im.header["Height"],
    qx, qy = im.header["PixelX"], im.header["PixelY"]
    if orgx is None: orgx = im.header["BeamX"]/qx
    if orgy is None: orgy = im.header["BeamY"]/qy
    if wavelength is None: wavelength = im.header["Wavelength"]
    if distance is None: distance = im.header["Distance"]
    friedel = "FALSE" if anomalous else "TRUE"
    sensor_thickness = 0 # FIXME

    if im.header["ImageType"] == "marccd":
        detector = "CCDCHESS MINIMUM_VALID_PIXEL_VALUE= 1 OVERLOAD= 65500"
    elif im.header["ImageType"] == "raxis":
        detector = "RAXIS MINIMUM_VALID_PIXEL_VALUE= 0  OVERLOAD= 2000000"
        distance *= -1
    elif im.header["ImageType"] == "minicbf":
        detector = "PILATUS MINIMUM_VALID_PIXEL_VALUE=0 OVERLOAD= 1048576"
        sensor_thickness = sensor_thickness_from_minicbf(img_files[0])
    elif im.header["ImageType"] == "adsc":
        detector = "ADSC MINIMUM_VALID_PIXEL_VALUE= 1 OVERLOAD= 65000"
    elif im.header["ImageType"] == "mscccd":
        detector = "SATURN MINIMUM_VALID_PIXEL_VALUE= 1 OVERLOAD= 262112" # XXX Should read header!!
        distance *= -1
    elif is_eiger_hdf5:
        detector = "EIGER MINIMUM_VALID_PIXEL_VALUE=0 OVERLOAD= %d" % im.header["Overload"]
        sensor_thickness = im.header["SensorThickness"]

    if minpk is not None: extra_kwds.append(" MINPK= %.2f" % minpk)
    for r1, r2 in exclude_resolution_range:
        if r1 < r2: r1, r2 = r2, r1
        extra_kwds.append(" EXCLUDE_RESOLUTION_RANGE= %.3f %.3f" % (r1, r2))

    extra_kwds = "\n".join(extra_kwds)
    inp_str = """\
 JOB= XYCORR INIT COLSPOT IDXREF DEFPIX INTEGRATE CORRECT
 ORGX= %(orgx).2f ORGY= %(orgy).2f
 DETECTOR_DISTANCE= %(distance).2f
 OSCILLATION_RANGE= %(osc_range).3f
 X-RAY_WAVELENGTH= %(wavelength).5f
 NAME_TEMPLATE_OF_DATA_FRAMES= %(template)s
! REFERENCE_DATA_SET=xxx/XDS_ASCII.HKL ! e.g. to ensure consistent indexing
 DATA_RANGE= %(data_range)s
 SPOT_RANGE= %(spot_range)s
! BACKGROUND_RANGE=1 10 ! rather use defaults (first 5 degree of rotation)
 FRIEDEL'S_LAW= %(friedel)s
 DELPHI= %(delphi).2f
! parameters specifically for this detector and beamline:
 DETECTOR= %(detector)s
 SENSOR_THICKNESS= %(sensor_thickness).2f
! attention CCD detectors: for very high resolution (better than 1A) make sure to specify SILICON
! as about 32* what CORRECT.LP suggests (absorption of phosphor is much higher than that of silicon)
 NX= %(nx)s NY= %(ny)s  QX= %(qx)s  QY= %(qy)s ! to make CORRECT happy if frames are unavailable
 ROTATION_AXIS= %(rotation_axis)s
%(extra_kwds)s
""" % locals()

    # XXX Really, really BAD idea!!
    # Synchrotron can have R-AXIS, and In-house detecotr can have horizontal goniometer..!!
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
        inp_str += """\
 INCIDENT_BEAM_DIRECTION= 0 0 1
 FRACTION_OF_POLARIZATION= 0.98   ! better value is provided by beamline staff!
 POLARIZATION_PLANE_NORMAL= 0 1 0
"""

    if not minimum:
        inp_str += """\
 SPACE_GROUP_NUMBER= %(sgnum)d                   ! 0 if unknown
 UNIT_CELL_CONSTANTS= %(cell)s ! put correct values if known
 INCLUDE_RESOLUTION_RANGE=50 0  ! after CORRECT, insert high resol limit; re-run CORRECT

 TRUSTED_REGION=0.00 1.4  ! partially use corners of detectors; 1.41421=full use
 VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS=6000. 30000. ! often 7000 or 8000 is ok
 STRONG_PIXEL=4           ! COLSPOT: only use strong reflections (default is 3)
 MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT=3 ! default of 6 is sometimes too high
! close spots: reduce SEPMIN and CLUSTER_RADIUS from their defaults of 6 and 3, e.g. to 4 and 2
! for bad or low resolution data remove the "!" in the following line:
 REFINE(IDXREF)=CELL BEAM ORIENTATION AXIS ! DISTANCE POSITION
 REFINE(INTEGRATE)= DISTANCE POSITION BEAM ORIENTATION ! AXIS CELL
! REFINE(CORRECT)=CELL BEAM ORIENTATION AXIS DISTANCE POSITION ! Default is: refine everything

!used by DEFPIX and CORRECT to exclude ice-reflections / ice rings - uncomment if necessary
!EXCLUDE_RESOLUTION_RANGE= 3.93 3.87 !ice-ring at 3.897 Angstrom
!EXCLUDE_RESOLUTION_RANGE= 3.70 3.64 !ice-ring at 3.669 Angstrom
!EXCLUDE_RESOLUTION_RANGE= 3.47 3.41 !ice-ring at 3.441 Angstrom
!EXCLUDE_RESOLUTION_RANGE= 2.70 2.64 !ice-ring at 2.671 Angstrom
!EXCLUDE_RESOLUTION_RANGE= 2.28 2.22 !ice-ring at 2.249 Angstrom
!EXCLUDE_RESOLUTION_RANGE= 2.102 2.042 !ice-ring at 2.072 Angstrom - strong
!EXCLUDE_RESOLUTION_RANGE= 1.978 1.918 !ice-ring at 1.948 Angstrom - weak
!EXCLUDE_RESOLUTION_RANGE= 1.948 1.888 !ice-ring at 1.918 Angstrom - strong
!EXCLUDE_RESOLUTION_RANGE= 1.913 1.853 !ice-ring at 1.883 Angstrom - weak
!EXCLUDE_RESOLUTION_RANGE= 1.751 1.691 !ice-ring at 1.721 Angstrom - weak
""" % dict(sgnum=sgnum, cell=cell_str)
    if im.header["ImageType"] == "minicbf":
        inp_str += """\
 NUMBER_OF_PROFILE_GRID_POINTS_ALONG_ALPHA/BETA= 13 ! Default is 9 - Increasing may improve data
 NUMBER_OF_PROFILE_GRID_POINTS_ALONG_GAMMA= 13      ! accuracy, particularly if finely-sliced on phi,
!                                                   and does not seem to have any downsides.
"""
        if nx == 1475:
            if 1:#! grep -q Flat_field tmp2: #XXX FIXME
                inp_str += """\
! the following specifications are for a detector _without_ proper
! flat_field correction; they cut away one additional pixel adjacent
! to each UNTRUSTED_RECTANGLE
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
            else:
                inp_str += """\
!EXCLUSION OF VERTICAL DEAD AREAS OF THE PILATUS 2M DETECTOR
 UNTRUSTED_RECTANGLE= 487  495     0 1680
 UNTRUSTED_RECTANGLE= 981  989     0 1680
!EXCLUSION OF HORIZONTAL DEAD AREAS OF THE PILATUS 2M DETECTOR
 UNTRUSTED_RECTANGLE=   0 1476   195  213
 UNTRUSTED_RECTANGLE=   0 1476   407  425
 UNTRUSTED_RECTANGLE=   0 1476   619  637
 UNTRUSTED_RECTANGLE=   0 1476   831  849
 UNTRUSTED_RECTANGLE=   0 1476  1043 1061
 UNTRUSTED_RECTANGLE=   0 1476  1255 1273
 UNTRUSTED_RECTANGLE=   0 1476  1467 1485
"""

        elif nx == 2463:
            # Pilatus 6M
            # FIXME: here we could test if a Flat_field correction was applied like we do for 2M
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
    if is_eiger_hdf5 and nx == 3110 and ny == 3269:
            # Eiger 9M
            inp_str += """\
 UNTRUSTED_RECTANGLE= 1029 1042 0 3269
 UNTRUSTED_RECTANGLE= 2069 2082 0 3269
 UNTRUSTED_RECTANGLE= 0 3110  513  553
 UNTRUSTED_RECTANGLE= 0 3110 1064 1104
 UNTRUSTED_RECTANGLE= 0 3110 1615 1655
 UNTRUSTED_RECTANGLE= 0 3110 2166 2206
 UNTRUSTED_RECTANGLE= 0 3110 2717 2757
"""

    return inp_str
