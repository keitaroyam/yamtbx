"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals
import os
import shutil
import traceback
import subprocess
import sys
import re
import pickle
import numpy
import glob
import time
import tempfile

from yamtbx.dataproc.xds import xds_inp
from yamtbx.dataproc.xds import modify_xdsinp
from yamtbx.dataproc.xds.xparm import XPARM
from yamtbx.dataproc.xds import integrate_hkl_as_flex
from yamtbx.dataproc.xds import xds_ascii
from yamtbx.dataproc.auto.command_line.run_all_xds_simple import run_xds, try_indexing_hard
from yamtbx.dataproc import crystfel
from yamtbx.dataproc import cbf

import iotbx.phil
import libtbx.phil
from libtbx.utils import Sorry
from libtbx import easy_mp
from libtbx.utils import multi_out
from libtbx.utils import null_out
from cctbx.crystal import reindex
from cctbx import crystal
from cctbx import uctbx


master_params_str = """
topdir = None
 .type = path
lstin = None
 .type = path
 .help = list of input files
nproc = 1
 .type = int
 .help = number of processors
split_num = 100
 .type = int
 .help = split working directory to prevent having too many directories 
anomalous = false
 .type = bool
minpk = 0
 .type = float
 .help = MINPK= for XDS
profile_fitting = true
 .type = bool
 .help = PROFILE_FITTING= for XDS
cell = 10 20 30 90 90 90
 .type = floats(size=6)
sgnum = 0
 .type = int
reintegrate = false
 .type = bool
 .help = Re-integrate based on post-refined parameters (by single image integration)
tryhard = False
 .type = bool
 .help = call try_indexing_hard()
refine_correct = True
 .type = bool
 .help = controls REFINE(CORRECT)=
d_min = 0
 .type = float
d_max = 50
 .type = float
value_range_for_trusted_detector_pixels = 6000. 30000
 .type = floats(size=2)
 .help = VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS=.
strong_pixel = 4
 .type = float
minimum_number_of_pixels_in_a_spot = 3
 .type = int
sepmin = 6
 .type = float
cluster_radius = 3
 .type = float
index_error = 0.05
 .type = float
index_magnitude= 8
 .type = int
index_quality = 0.8
 .type = float
refine_idxref = POSITION *BEAM *AXIS *ORIENTATION *CELL SEGMENT
 .type = choice(multi=True)
refine_integrate = POSITION *BEAM *AXIS *ORIENTATION *CELL SEGMENT
 .type = choice(multi=True)
idxref_d_min = 0
 .type = float

osc_range = None
 .type = float
 .help = Override osc range
orgx = None
 .type = float
orgy = None
 .type = float
distance = None
 .type = float
 .help = Override camera distance in mm
rotation_axis = None
 .type = floats(size=3)
 .help = Override ROTATION_AXIS=
incident_beam_direction = 0,0,1
 .type = floats(size=3)
 .help = Override INCIDENT_BEAM_DIRECTION=
extra_inp_str = None
 .type = str
 .multiple = True

sfx {
 cheetah_mpccd = false
  .type = bool
  .help = SFX setting.
 crystfel_geom = None
  .type = path
 osc_range = 0.1
  .type = float
 clen = None
  .type = float
  .help = Override clen (in m).
}

checkcell {
 check = true
  .type = bool
 tol_length = 0.2
  .type = float
  .help = relative_length_tolerance
 tol_angle = 10
  .type = float
  .help = absolute_angle_tolerance in degree
}

remove_heavy_unuseful_files = false
 .type = bool
keep_file = None
 .type = str
 .multiple = true
pickle_hkl = true
 .type = bool
light_pickle = false
 .type = bool
tmpdir = None
 .type = path
 .help = temporary directory for xds run
"""

bool_to_str = lambda x: ("FALSE", "TRUE")[int(x)]

class ProcFailed(Exception): pass

def sfx_xds_inp(img_in, xdsinp, params):
    import pycbf

    os.symlink(img_in, os.path.join(os.path.dirname(xdsinp), "data_10000.cbf"))

    header = cbf.get_pilatus_header(img_in)
    sensor = [x for x in header.splitlines() if "Wavelength" in x][0]
    wavelength, unit = sensor[sensor.index("Wavelength ")+len("Wavelength "):].split()
    assert unit == "A"

    orgx, orgy = params.orgx, params.orgy
    if orgx is None: orgx = 0
    if orgy is None: orgy = 0

    if params.rotation_axis is not None:
        rotation_axis = " ".join(["%.2f"%x for x in params.rotation_axis])
    else:
        rotation_axis = "1 0 0"

    # Reference: http://strucbio.biologie.uni-konstanz.de/xdswiki/index.php/Generate_XDS.INP
    inpstr = """\
JOB= COLSPOT IDXREF DEFPIX INTEGRATE CORRECT
ORGX= %(orgx).2f ORGY= %(orgy).2f
DETECTOR_DISTANCE= 0
OSCILLATION_RANGE= %(osc_range).5f
X-RAY_WAVELENGTH= %(wavelength)f
NAME_TEMPLATE_OF_DATA_FRAMES=data_?????.cbf
DATA_RANGE=10000 10000

SPACE_GROUP_NUMBER= 0 ! 16
UNIT_CELL_CONSTANTS= 48.3392 77.1339 84.7324 90 90 90
INCLUDE_RESOLUTION_RANGE=50 0  ! after CORRECT, insert high resol limit; re-run CORRECT

FRIEDEL'S_LAW=FALSE     ! This acts only on the CORRECT step

! parameters with changes wrt default values:
TRUSTED_REGION=0.0 2 ! partially use corners of detector (0 1.4143: use all pixels)
VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS=6000. 30000. ! often 7000 or 8000 is ok
STRONG_PIXEL=2           ! COLSPOT: only use strong reflections (default is 3)
MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT=2 ! default of 6 is sometimes too high
! close spots/long cell axis: reduce SEPMIN and CLUSTER_RADIUS from their defaults of 6 and 3
! SEPMIN=4  CLUSTER_RADIUS=2 ! should be default for Pilatus and other detectors with low PSF
! since XDS 01-MAR-2015, POSITION is used instead of DISTANCE. Older versions do it the other way round.
REFINE(IDXREF)=CELL BEAM ORIENTATION  AXIS ! DISTANCE POSITION
REFINE(INTEGRATE)= DISTANCE POSITION BEAM ORIENTATION ! AXIS CELL
! REFINE(CORRECT)=CELL BEAM ORIENTATION AXIS DISTANCE POSITION ! Default is: refine everything

! parameters specifically for this detector and beamline:
DETECTOR= PILATUS MINIMUM_VALID_PIXEL_VALUE=0 OVERLOAD= 1048576  !PILATUS
SENSOR_THICKNESS= 0.05
! attention CCD detectors: for very high resolution (better than 1A) make sure to specify SILICON
! as about 32* what CORRECT.LP suggests (absorption of phosphor is much higher than that of silicon)
NX= 512 NY= 8192  QX= 0.05   QY= 0.05  ! to make CORRECT happy if frames are unavailable

ROTATION_AXIS=%(rotation_axis)
 DIRECTION_OF_DETECTOR_X-AXIS=1 0 0
 DIRECTION_OF_DETECTOR_Y-AXIS=0 1 0
INCIDENT_BEAM_DIRECTION=0 0 1
FRACTION_OF_POLARIZATION=0.98   ! better value is provided by beamline staff!
POLARIZATION_PLANE_NORMAL=0 1 0

""" % dict(osc_range=params.sfx.osc_range,
           wavelength=float(wavelength),
           orgx=orgx, orgy=orgy,
           rotation_axis=rotation_axis)

    from yamtbx.dataproc.crystfel.command_line.geom_for_xds import geom_to_xdsinp_str

    geom = crystfel.geom.Geomfile(params.sfx.crystfel_geom)
    if params.sfx.clen is not None: geom.clen = params.sfx.clen

    inpstr += geom_to_xdsinp_str(geom)

    return inpstr
# sfx_xds_inp()

def remove_heavy_unuseful_files(wdir, keep_files=[]):
    files = [x for x in glob.glob(os.path.join(wdir, "*.cbf")) if "data_" not in os.path.basename(x)]
    files += glob.glob(os.path.join(wdir, "INTEGRATE.HKL"))

    for f in files:
        if not os.path.basename(f) in keep_files:
            os.remove(f)
# remove_heavy_unuseful_files

def check_reidx(params, xs):
    assert params.checkcell.check and params.sgnum > 0
    xsref = crystal.symmetry(params.cell, params.sgnum)
    cosets = reindex.reindexing_operators(xsref, xs,
                                          params.checkcell.tol_length, params.checkcell.tol_angle)

    if cosets.double_cosets is None:
        return None

    return cosets.combined_cb_ops()[0]
# check_reidx()

def check_cell(params, xs):
    assert params.checkcell.check and params.sgnum > 0
    cellref = uctbx.unit_cell(params.cell)
    return cellref.is_similar_to(xs.unit_cell(),
                                     relative_length_tolerance=params.checkcell.tol_length,
                                     absolute_angle_tolerance=params.checkcell.tol_angle)
# check_cell()

def xds_sequence(img_in, topdir, data_root_dir, params):
    relpath = os.path.relpath(os.path.dirname(img_in), data_root_dir)
    workdir = os.path.abspath(os.path.join(topdir, relpath, os.path.splitext(os.path.basename(img_in))[0]))
    print(workdir)
    frame_num = None

    if not os.path.exists(img_in):
        if "<>" in img_in:
            img_in, num = img_in.split("<>")
            frame_num = int(num)
            if not os.path.exists(img_in):
                print("Error: Not found!!", img_in)
                return
            workdir += "_%.6d" % frame_num
            assert img_in.endswith(".h5")
        else:
            for ext in (".bz2", ".gz", ".xz"):
                if os.path.exists(img_in+ext):
                    img_in += ext
                    break

    if params.tmpdir is not None:
        tmpdir = tempfile.mkdtemp(prefix="xds", dir=params.tmpdir)
    else:
        tmpdir = workdir
        
    if not os.path.exists(tmpdir): os.makedirs(tmpdir)

    xparm = os.path.join(tmpdir, "XPARM.XDS")
    xdsinp = os.path.join(tmpdir, "XDS.INP")
    integrate_hkl = os.path.join(tmpdir, "INTEGRATE.HKL")
    decilog = open(os.path.join(tmpdir, "decision.log"), "w")

    try:
        print("Paramters:", file=decilog)
        libtbx.phil.parse(master_params_str).format(params).show(out=decilog, prefix=" ")

        print("\nStarting at %s" % time.strftime("%Y-%m-%d %H:%M:%S"), file=decilog)

        # Prepare XDS.INP
        if params.sfx.cheetah_mpccd:
            xdsinp_str = sfx_xds_inp(img_in, xdsinp, params)
        else:
            if frame_num is not None: # Eiger
                from yamtbx.dataproc import eiger
                img_in_work = os.path.join(os.path.dirname(xdsinp), "data_10000.cbf")
                eiger.extract_to_minicbf(img_in, frame_num, img_in_work)
            else:
                ext = os.path.splitext(img_in)[1]
                img_in_work = os.path.join(os.path.dirname(xdsinp), "data_10000" + ext)
                os.symlink(img_in, img_in_work)
            xdsinp_str = xds_inp.generate_xds_inp(img_files=[img_in_work],
                                                  inp_dir=tmpdir,
                                                  reverse_phi=True, anomalous=params.anomalous,
                                                  spot_range=None, minimum=False,
                                                  crystal_symmetry=None,
                                                  integrate_nimages=None,
                                                  osc_range=params.osc_range,
                                                  orgx=params.orgx, orgy=params.orgy,
                                                  rotation_axis=params.rotation_axis,
                                                  distance=params.distance)

        ofs = open(xdsinp, "w")
        ofs.write(xdsinp_str)
        ofs.close()
        # Start processing
        modify_xdsinp(xdsinp, inp_params=[("JOB", "XYCORR INIT COLSPOT IDXREF"),
                                          ("MAXIMUM_NUMBER_OF_PROCESSORS", "1"),
                                          ("MINPK", "%.2f"%params.minpk),
                                          ("PROFILE_FITTING", bool_to_str(params.profile_fitting)),
                                          ("STRONG_PIXEL", "%.2f"%params.strong_pixel),
                                          ("MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT", "%d"%params.minimum_number_of_pixels_in_a_spot),
                                          ("SEPMIN", "%.2f"%params.sepmin),
                                          ("CLUSTER_RADIUS", "%.2f"%params.cluster_radius),
                                          ("INDEX_ERROR", "%.4f"%params.index_error),
                                          ("INDEX_MAGNITUDE", "%d"%params.index_magnitude),
                                          ("INDEX_QUALITY", "%.2f"%params.index_quality),
                                          ("REFINE(IDXREF)", " ".join([s.upper() for s in params.refine_idxref])),
                                          ("REFINE(INTEGRATE)", " ".join([s.upper() for s in params.refine_integrate])),
                                          ("INCLUDE_RESOLUTION_RANGE", "%.2f %.2f" % (params.d_max, params.idxref_d_min)),
                                          ("VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS", "%.1f %.1f" % tuple(params.value_range_for_trusted_detector_pixels)),
                                          ("INCIDENT_BEAM_DIRECTION", "%.6f %.6f %.6f" % tuple(params.incident_beam_direction))
                                          ])

        if len(params.extra_inp_str) > 0:
            ofs = open(xdsinp, "a")
            ofs.write("\n%s\n" % "\n".join(params.extra_inp_str))
            ofs.close()

        run_xds(wdir=tmpdir, show_progress=False)

        if params.tryhard:
            try_indexing_hard(tmpdir, show_progress=True, decilog=decilog)

            # If Cell hint exists, try to use it..
            if params.sgnum > 0:
                flag_try_cell_hint = False
                if not os.path.isfile(xparm):
                    flag_try_cell_hint = True
                else:
                    xsxds = XPARM(xparm).crystal_symmetry()
                    reidx_op = check_reidx(params, xsxds)
                    if reidx_op is None: flag_try_cell_hint = True

                if flag_try_cell_hint:
                    print(" Worth trying to use prior cell for indexing.", file=decilog)
                    modify_xdsinp(xdsinp, inp_params=[("JOB", "IDXREF"),
                                                      ("UNIT_CELL_CONSTANTS",
                                                       " ".join(["%.3f"%x for x in params.cell])),
                                                      ("SPACE_GROUP_NUMBER", "%d"%params.sgnum),
                                                      ])
                    run_xds(wdir=tmpdir, show_progress=False)
                    modify_xdsinp(xdsinp, inp_params=[("SPACE_GROUP_NUMBER", "0"),
                                                      ])

        if not os.path.isfile(xparm):
            raise ProcFailed("Indexing failed")

        if params.checkcell.check and params.sgnum > 0:
            xsxds = XPARM(xparm).crystal_symmetry()
            reidx_op = check_reidx(params, xsxds)
            if reidx_op is None:
                raise ProcFailed("Incompatible cell. Indexing failed.")

            if not reidx_op.is_identity_op():
                print("Re-idxref to match reference symmetry.", file=decilog)
                xsxds_cb = xsxds.change_basis(reidx_op) # Really OK??
                modify_xdsinp(xdsinp, inp_params=[("JOB", "IDXREF"),
                                                  ("SPACE_GROUP_NUMBER", "%d"%params.sgnum),
                                                  ("UNIT_CELL_CONSTANTS", " ".join(["%.3f"%x for x in xsxds_cb.unit_cell().parameters()]))
                                                  ])
                run_xds(wdir=tmpdir, show_progress=False)

        # Final check
        if params.checkcell.check and params.sgnum > 0:
            if not check_cell(params, xsxds):
                raise ProcFailed(" Incompatible cell. Stop before integration.")

        modify_xdsinp(xdsinp, inp_params=[("INCLUDE_RESOLUTION_RANGE", "%.2f %.2f" % (params.d_max, params.d_min)),
                                          ])

        # To Integration
        modify_xdsinp(xdsinp, inp_params=[("JOB", "DEFPIX INTEGRATE"),])
        run_xds(wdir=tmpdir, show_progress=False)

        if not os.path.isfile(integrate_hkl):
            raise ProcFailed("Integration failed.")

        # Determine unit cell in CORRECT
        if params.refine_correct:
            tmp = [("REFINE(CORRECT)", "ALL"),
                   ("UNIT_CELL_CONSTANTS", " ".join(["%.3f"%x for x in params.cell]))]
        else:
            # XXX What if CELL is refined in INTEGRATE?
            xsxds = XPARM(xparm).crystal_symmetry()
            reidx_op = check_reidx(params, xsxds)
            if reidx_op is None:
                raise ProcFailed(" Incompatible cell. Failed before CORRECT.")

            xsxds_cb = xsxds.change_basis(reidx_op) # Really OK??
            tmp = [("REFINE(CORRECT)", ""),
                   ("UNIT_CELL_CONSTANTS", " ".join(["%.3f"%x for x in xsxds_cb.unit_cell().parameters()]))]

        # PEAK-corrected INTEGRATE.HKL
        ihk = os.path.join(tmpdir, "INTEGRATE.HKL")
        ihk_full = os.path.join(tmpdir, "INTEGRATE_full.HKL")
        ihk_part = os.path.join(tmpdir, "INTEGRATE_part.HKL")
        inhkl = integrate_hkl_as_flex.reader(ihk, [], False)
        inhkl.write_peak_corrected(ihk_part)
        os.rename(ihk, ihk_full)
        
        modify_xdsinp(xdsinp, inp_params=[("JOB", "CORRECT"),
                                          ("DATA_RANGE", "1 20000"),
                                          ("CORRECTIONS", ""),
                                          ("NBATCH", "1"),
                                          ("SPACE_GROUP_NUMBER", "%d"%params.sgnum)] + tmp)

        xac = os.path.join(tmpdir, "XDS_ASCII.HKL")
        xac_full = os.path.join(tmpdir, "XDS_ASCII_full.HKL")
        xac_part = os.path.join(tmpdir, "XDS_ASCII_part.HKL")

        # CORRECT for full
        os.symlink(ihk_full, ihk)
        run_xds(wdir=tmpdir, comm="xds", show_progress=False)
        if os.path.isfile(xac):
            os.rename(xac, xac_full)
            os.rename(os.path.join(tmpdir, "CORRECT.LP"),
                      os.path.join(tmpdir, "CORRECT_full.LP"))
        os.remove(ihk)

        # CORRECT for part
        os.symlink(ihk_part, ihk)
        run_xds(wdir=tmpdir, comm="xds", show_progress=False)
        if os.path.isfile(xac):
            os.rename(xac, xac_part)
            os.rename(os.path.join(tmpdir, "CORRECT.LP"),
                      os.path.join(tmpdir, "CORRECT_part.LP"))
        os.remove(ihk)

        if params.pickle_hkl:
            for f in [x for x in (xac_part, xac_full) if os.path.isfile(x)]:
                print("Pickling %s" % os.path.basename(f), file=decilog)
                x = xds_ascii.XDS_ASCII(f, log_out=decilog)
                if params.light_pickle: x.xd, x.yd, x.zd, x.rlp, x.corr = None, None, None, None, None # To make reading faster
                pickle.dump(x, open(f+".pkl", "wb"), -1)
        if params.pickle_hkl:
            for f in [x for x in (ihk_part, ihk_full) if os.path.isfile(x)]:
                print("Pickling %s" % os.path.basename(f), file=decilog)
                inhkl = integrate_hkl_as_flex.reader(f, read_columns=["IOBS","SIGMA","XCAL","YCAL","RLP","PEAK","MAXC"])
                pickle.dump(inhkl, open(f+".pkl", "wb"), -1)

    except ProcFailed as e:
        print("Processing failed: %s" % e.message, file=decilog)
    except:
        print("Uncaught exception:", file=decilog)
        print(traceback.format_exc(), file=decilog)
    finally:
        if params.remove_heavy_unuseful_files: remove_heavy_unuseful_files(tmpdir, keep_files=params.keep_file)
        print("\nFinished at %s" % time.strftime("%Y-%m-%d %H:%M:%S"), file=decilog)
        decilog.close()
        if tmpdir != workdir: shutil.move(tmpdir, workdir)

# xds_sequence()

def get_file_list(lstin):
    ret = []
    for l in open(lstin):
        l = l.strip()
        if "#" in l: l = l[:l.index("#")]
        if l: ret.append(l)
    return ret
# get_file_list()

def run(params):
    input_files = get_file_list(params.lstin)
    top_dirs = [os.path.join(params.topdir, "split_%.4d" % (i//params.split_num+1)) for i in range(len(input_files))]

    data_root_dir = os.path.dirname(os.path.commonprefix(input_files))

    fun_local = lambda x: xds_sequence(x[0], x[1],  data_root_dir, params)

    #for arg in input_files: fun_local(arg)

    easy_mp.pool_map(fixed_func=fun_local,
                     args=list(zip(input_files, top_dirs)),
                     processes=params.nproc)
# run()

def run_from_args(argv):
    cmdline = iotbx.phil.process_command_line(args=argv,
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args

    for arg in args:
        if os.path.isdir(arg) and params.topdir is None:
            params.topdir = arg
        if os.path.isfile(arg) and params.lstin is None:
            params.lstin = arg

    assert params.lstin is not None

    run(params)
# run_from_args()

if __name__ == "__main__":
    import sys
    run_from_args(sys.argv[1:])
