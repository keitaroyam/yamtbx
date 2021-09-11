"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import unicode_literals
# Only for SPring-8 MAR detectors
xds_inp_template = """\
JOB= XYCORR INIT COLSPOT
ORGX= %(orgx).2f ORGY= %(orgy).2f
DETECTOR_DISTANCE= %(distance).3f
OSCILLATION_RANGE= %(osc_range).3f
X-RAY_WAVELENGTH= %(wavelength).5f
NAME_TEMPLATE_OF_DATA_FRAMES= %(template)s
DATA_RANGE= %(framenum)d %(framenum)d
SPOT_RANGE= %(framenum)d %(framenum)d

! COLSPOT paramters
STRONG_PIXEL= %(strong_pixel).2f !4 ! COLSPOT: only use strong reflections (default is 3)
MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT= %(min_pixels)d !3 ! default of 6 is sometimes too high
BACKGROUND_PIXEL= %(background_pixel)f
MAXIMUM_NUMBER_OF_STRONG_PIXELS= %(max_strong_pixels)d
SPOT_MAXIMUM-CENTROID= %(spot_maximum_centroid)f
REFLECTING_RANGE= %(reflecting_range)f

VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS= %(defpix_trusted1).1f %(defpix_trusted2).1f

DETECTOR= CCDCHESS MINIMUM_VALID_PIXEL_VALUE= 1 OVERLOAD= 65500
NX= %(nx)d NY= %(ny)d  QX= %(qx).6f QY= %(qy).6f
ROTATION_AXIS=-1 0 0
DIRECTION_OF_DETECTOR_X-AXIS=1 0 0
DIRECTION_OF_DETECTOR_Y-AXIS=0 1 0
INCIDENT_BEAM_DIRECTION=0 0 1
FRACTION_OF_POLARIZATION=0.98
POLARIZATION_PLANE_NORMAL=0 1 0

"""

from yamtbx.dataproc import XIO
from yamtbx.dataproc import dataset
from yamtbx.dataproc.xds import modify_xdsinp, xparm
from yamtbx.util import call, rotate_file
import os
import sys
import math
import re
import shutil
import tempfile

def resol_to_radius(d, distance, wavelength):
    theta = math.asin(wavelength/2./d)
    return distance * math.tan(2.0*theta)

def coord_to_resol(x, y, header):
    distance, wavelength = header["Distance"], header["Wavelength"]
    qx, qy = header["PixelX"], header["PixelY"]
    orgx, orgy = header["BeamX"]/qx, header["BeamY"]/qy

    l = math.sqrt(((x-orgx)*qx)**2 + ((y-orgy)*qy)**2)
    theta = math.atan(l/distance)/2
    return wavelength / 2. / math.sin(theta)

def res_range_to_trusted_region(dmin, dmax, header):
    ret = [0, 1.4]

    nx, ny = header["Width"], header["Height"],
    distance, wavelength = header["Distance"], header["Wavelength"]
    qx, qy = header["PixelX"], header["PixelY"]

    orgx, orgy = nx/2., ny/2.

    if dmin is not None:
        rad = resol_to_radius(dmin, distance, wavelength)
        ret[1] = max(rad/qx/(nx-orgx), rad/qx/(ny-orgy))
    if dmax is not None:
        rad = resol_to_radius(dmax, distance, wavelength)
        ret[0] = max(rad/qy/(nx-orgx), rad/qy/(ny-orgy)) # should be min?

    return " ".join(["%.2f"%x for x in ret])
#  res_range_to_trusted_region()

def res_range_for_xds(dmin, dmax):
    ret = [100., 0.]
    if dmin is not None:
        ret[1] = dmin
    if dmax is not None:
        ret[0] = dmax
    return " ".join(["%.2f"%x for x in ret])
# res_range_for_xds()

def find_spots(img_file, params):
    """
    Use XDS to locate spots.
    If params.xds.do_defpix is true, DEFPIX will be run. For DEFPIX, dummy XPARM.XDS is needed. Thanks to DEFPIX, we can mask beam stop shadow and remove areas outside the resolution range.
    If false, we need to set TRUSTED_REGION to exclude regions outside resolution range, but it is independent of beam center. Maybe we should remove spots outside the resolution range after XDS run?
    """

    # Test if ramdisk available
    if os.path.isdir("/dev/shm"):
        work_dir = tempfile.mkdtemp(prefix="shika_x_"+os.path.splitext(os.path.basename(img_file))[0], dir="/dev/shm")
    else:
        work_dir = os.path.join(params.work_dir, "xds_"+os.path.splitext(os.path.basename(img_file))[0])

    xdsinp = os.path.join(work_dir, "XDS.INP")
    spot_xds = os.path.join(work_dir, "SPOT.XDS")

    if not os.path.exists(work_dir):
        os.mkdir(work_dir)

    template_str, min_frame, max_frame = dataset.group_img_files_template([img_file])[0]

    im = XIO.Image(img_file)

    # Remove lines if None (means to use default)
    params_maps = [("strong_pixel", "STRONG_PIXEL="),
                   ("minimum_number_of_pixels_in_a_spot", "MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT="),
                   ("background_pixel", "BACKGROUND_PIXEL="),
                   ("maximum_number_of_strong_pixels", "MAXIMUM_NUMBER_OF_STRONG_PIXELS="),
                   ("spot_maximum_centroid", "SPOT_MAXIMUM-CENTROID="),
                   ("reflecting_range", "REFLECTING_RANGE="),
                   ("value_range_for_trusted_detector_pixels", "VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS="),
                   ]
    tmp = xds_inp_template.splitlines()
    for p, x in params_maps:
        if getattr(params.xds, p) is None:
            tmp = [s for s in tmp if not s.startswith(x)]
    inp_template = "\n".join(tmp)

    # Prepare XDS.INP
    inp_str = inp_template%dict(template=os.path.relpath(template_str, work_dir),
                                framenum=min_frame,
                                orgx=im.header["BeamX"]/im.header["PixelX"],
                                orgy=im.header["BeamY"]/im.header["PixelY"],
                                distance=im.header["Distance"],
                                osc_range=im.header["PhiWidth"],
                                wavelength=im.header["Wavelength"],
                                strong_pixel=params.xds.strong_pixel,
                                min_pixels=params.xds.minimum_number_of_pixels_in_a_spot,
                                background_pixel=params.xds.background_pixel,
                                max_strong_pixels=params.xds.maximum_number_of_strong_pixels,
                                spot_maximum_centroid=params.xds.spot_maximum_centroid,
                                reflecting_range=params.xds.reflecting_range,
                                nx=im.header["Width"], ny=im.header["Height"],
                                qx=im.header["PixelX"], qy=im.header["PixelY"],
                                defpix_trusted1=params.xds.value_range_for_trusted_detector_pixels[0],
                                defpix_trusted2=params.xds.value_range_for_trusted_detector_pixels[1]
                                )

    open(xdsinp, "w").write(inp_str)

    if params.xds.do_defpix:
        xp = xparm.XPARM()
        xp.set_info_from_xdsinp(xdsinp)
        open(os.path.join(work_dir, "XPARM.XDS"), "w").write(xp.xparm_str())
        modify_xdsinp(xdsinp, inp_params=[("JOB", "XYCORR INIT DEFPIX"),
                                          ("INCLUDE_RESOLUTION_RANGE", res_range_for_xds(params.distl.res.outer, params.distl.res.inner))
                                          ])
        call("xds", wdir=work_dir, stdout=open(os.path.join(work_dir, "xds.log"), "w"))
        shutil.copy(os.path.join(work_dir, "BKGPIX.cbf"), os.path.join(work_dir, "BKGINIT.cbf"))
        modify_xdsinp(xdsinp, inp_params=[("JOB","COLSPOT")])
    else:
        modify_xdsinp(xdsinp, inp_params=[("TRUSTED_REGION", res_range_to_trusted_region(params.distl.res.outer, params.distl.res.inner, im.header))
                                          ])
        open(os.path.join(work_dir, "xds.log"), "w").write("")

    # Run XDS
    rotate_file(spot_xds)
    call("xds", wdir=work_dir, stdout=open(os.path.join(work_dir, "xds.log"), "a"))

    # Extract results
    spots = [] # (x, y, d, intensity)
    if os.path.isfile(spot_xds):
        for l in open(spot_xds):
            x, y, z, intensity = [float(x) for x in l.strip().split()]
            d = coord_to_resol(x, y, im.header)
            spots.append((x, y, d, intensity))

    # Delete dir
    shutil.rmtree(work_dir)

    return spots
        
    
