# LIBTBX_SET_DISPATCHER_NAME kamo.test_installation

"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function

from yamtbx import util
import os
import shutil
import traceback
import tempfile

def tst_jsdir():
    print("Testing location..")

    yamtbx_root = util.yamtbx_module_root()
    if not yamtbx_root:
        print("  Can't locate yamtbx module directory. NG")
        return False

    d3path = os.path.join(yamtbx_root, "dataproc/auto/js/d3-3.5.10")
    if not os.path.isdir(d3path):
        print("  Can't find d3-3.5.10 directory. Please check location of yamtbx. NG")
        return False
    
    print("  %s. OK" % yamtbx_root)
    return True
# tst_jsdir()

def tst_R():
    print("Testing R..")

    rcode, out, err = util.call("Rscript", '-e "print(cos(pi))"')
    if rcode != 0 or out.strip() != '[1] -1':
        print("  Rscript is not avaiable. NG")
        return False
    
    rcode, out, err = util.call("Rscript", '-e "library(rjson)"')
    if rcode != 0:
        print("  rjson is not installed. NG")
        return False
    
    print("  OK")
    return True
# tst_R()

def tst_xds():
    print("Testing XDS..")
    tmpdir = util.get_temp_local_dir("xdstest")
    rcode, out, err = util.call("xds_par", wdir=tmpdir)
    if tmpdir: shutil.rmtree(tmpdir) # just in case; xds shouldn't make any files

    if rcode != 0:
        print("  Not installed. NG")
        return False

    if "license expired" in out:
        print("  license expired. Get latest version. NG")
        return False

    print("  OK")
    return True
# tst_xds()

def tst_h5toxds():
    print("Testing H5ToXds..")
    rcode, out, err = util.call("H5ToXds")
    ignore_msg = "(You can ignore this if you don't process hdf5 files which usually mean Eiger data)"

    if rcode == 127: # 127 is "command not found".
        print("  Not installed. NG %s" % ignore_msg)
        return False

    import numpy
    from yamtbx.dataproc.eiger import make_dummy_h5_for_test
    from yamtbx.dataproc import cbf

    tmpdir = util.get_temp_local_dir("h5test")
    data = numpy.random.randint(0, 65535, size=100).astype(numpy.uint32).reshape((1,10,10))
    master_h5 = make_dummy_h5_for_test(tmpdir, data)
    rcode, out, err = util.call("H5ToXds", "%s 1 1.cbf"%os.path.basename(master_h5), wdir=tmpdir)
    cbfout = os.path.join(tmpdir, "1.cbf")
    if not os.path.isfile(cbfout):
        print("  H5ToXds exists, but not works. Probably Dectris original H5ToXds? Test it with real h5 file. %s" % ignore_msg)
        if out.strip():
            print("  -- stdout:")
            print(out)
        if err.strip():
            print("  -- stderr:")
            print(err)

        shutil.rmtree(tmpdir) 
        return False

    data_read, _, _ = cbf.load_minicbf_as_numpy(cbfout)
    shutil.rmtree(tmpdir)
    
    if numpy.all(data_read.flatten() == data.flatten()):
        print("  OK")
        return True
    else:
        print("  H5ToXds exists, but not correctly works. Probably Dectris original H5ToXds? Test it with real h5 file. %s" % ignore_msg)
        if out.strip():
            print("  -- stdout:")
            print(out)
        if err.strip():
            print("  -- stderr:")
            print(err)

        return False
# tst_xds()

def tst_xdsstat():
    print("Testing XDSSTAT..")
    rcode, out, err = util.call("xdsstat", stdin="\n")
    if rcode != 0:
        print("  Not installed. NG")
        return False

    if "XDSSTAT version" not in err:
        print("  Seems not working. NG")
        return False

    print("  OK")
    return True
# tst_xdsstat()

def tst_ccp4():
    print("Testing ccp4..")
    if "CCP4" not in os.environ or not os.path.isdir(os.environ["CCP4"]):
        print("  Not installed. NG")
        return False
    
    if not os.path.isfile(os.path.join(os.environ["CCP4"], "share/blend/R/blend0.R")):
        print("  BLEND is not available. NG")
        return False

    print("  OK")
    return True
# tst_ccp4()

def tst_dials():
    print("Testing dials (package)..")
    rcode, out, err = util.call("dials.version")

    if rcode != 0:
        print("  Not installed. NG")
        return False

    print("\n".join(["  "+x for x in out.splitlines()]))
    print("  OK")
    return True
# tst_dials()

def tst_dials_module():
    print("Testing dials (module)..")

    try:
        from dials.util.version import dials_version
        print("  %s installed. OK" % dials_version())
        return True
    except ImportError:
        print("  Not installed. NG")
        return False
# tst_dials_module()

def tst_dxtbx_eiger():
    print("Testing eiger hdf5 geometry recognition..")

    try:
        # due to https://github.com/cctbx/cctbx_project/issues/282
        import dxtbx.format.FormatHDFEigerNearlyNexus
    except ImportError:
        print("  FormatHDFEigerNearlyNexus does not exist. This implies h5-recorded geometry can be no longer extracted using dxtbx. I hope it works. OK")
        return True

    try:
        import h5py
        import uuid
        from dxtbx.format.nexus import NXmxReader
        from dxtbx.format.FormatHDFEigerNearlyNexus import EigerNXmxFixer
        from dxtbx.format.nexus import BeamFactory
        from dxtbx.format.nexus import DetectorFactory
        from dxtbx.format.nexus import GoniometerFactory
        import numpy

        tmpfd, tmpf = tempfile.mkstemp()
        os.close(tmpfd)

        h = h5py.File(tmpf, "w")
        h.create_group("/entry/instrument/detector/")
        h["/entry/instrument/detector/description"] = "Dectris Eiger"
        h["/entry/"].attrs["NX_class"] = "NXentry"
        g = h.create_group("entry/instrument/detector/detectorSpecific/detectorModule_000")
        g["countrate_correction_count_cutoff"] = 10
        g = h.create_group("/entry/instrument/beam")
        g["incident_wavelength"] = 1.
        g["incident_wavelength"].attrs["units"] = "angstrom"
        h["/entry/instrument"].attrs["NX_class"] = "NXinstrument"
        h["/entry/instrument/beam"].attrs["NX_class"] = "NXbeam"
        g = h.create_group("/entry/data")
        g.create_dataset("data_000001", (1, 3269, 3110), dtype=numpy.uint16)
        h["/entry/data"].attrs["NX_class"] = "NXdata"
        g = h.create_group("/entry/instrument/detector/geometry/orientation")
        g["value"] = (-1,0,0, 0,-1,0)
        g = h.create_group("/entry/instrument/detector/geometry/translation")
        g["distances"] = (0.11737501, 0.11992501, -0.18)
        h["/entry/instrument/detector"].attrs["NX_class"] = "NXdetector"
        h["/entry/instrument/detector/x_pixel_size"] = h["/entry/instrument/detector/y_pixel_size"] = 7.5e-5
        h["/entry/instrument/detector/sensor_material"] = "Si"
        h["/entry/instrument/detector/sensor_thickness"] = 4.5e-4
        h["/entry/instrument/detector/sensor_thickness"].attrs["units"] = "m"
        h["/entry/instrument/detector/count_time"] = 0.99998

        g = h.create_group("/entry/sample/goniometer")
        g["omega_range_average"] = 1.
        h["/entry/sample"].attrs["NX_class"] = "NXsample"

        h.close()

        temp_file = "tmp_master_%s.nxs" % uuid.uuid1().hex
        fixer = EigerNXmxFixer(tmpf, temp_file)
        reader = NXmxReader(handle=fixer.handle)

        entry = reader.entries[0]
        instrument = entry.instruments[0]
        detector = instrument.detectors[0]
        sample = entry.samples[0]
        beam = sample.beams[0]
        beam_model = BeamFactory(beam).model
        detector_model = DetectorFactory(detector, beam_model).model
        goniometer_model = GoniometerFactory(sample).model

        _test = detector_model[0].get_fast_axis() + detector_model[0].get_slow_axis() + beam_model.get_direction() + goniometer_model.get_rotation_axis()
        os.remove(tmpf)
        
        if _test == (1.0, 0.0, 0.0,  0.0, 1.0, 0.0,  0.0, 0.0, -1.0,  1.0, 0.0, 0.0):
            print("  OK")
            return True
        elif _test == (1.0, 0.0, 0.0,  0.0, -1.0, 0.0,  0.0, 0.0, 1.0,  1.0, 0.0, 0.0):
            print("  OK")
            return True
        else:
            print("  vectors=", _test)
    except:
        print(traceback.format_exc())

    print("  NG! You are using old cctbx. Use latest cctbx or environment of DIALS v1.6 or PHENIX 1.12 or newer.")
    return False
# tst_dxtbx_eiger()

#def tst_h5():
#    print "Testing hdf5..",   

def tst_adxv():
    print("Testing Adxv..")

    rcode, out, err = util.call("adxv", '-help')
    if rcode != 0:
        print("  Adxv is not avaiable. NG")
        return False
    
    print("  OK")
    return True
# tst_R()

def tst_scipy():
    print("Testing SciPy..")

    try: import scipy
    except ImportError:
        print("  Not installed. NG")
        return False

    v = scipy.version.full_version
    vv = list(map(int, v.split(".")))
    if vv < [0, 18, 1]:
        print("  SciPy version installed (%s) is too old. Use 0.18.1 or newer." % v)
        return False

    try: import scipy.optimize
    except ImportError:
        print("  Scipy %s installed, but does not work properly (import scipy.optimize failed). NG" % v)
        return False

    try: scipy.optimize.least_squares
    except AttributeError:
        print("  scipy.optimize.least_squares is not available. Update the version. NG")
        return False

    print("  %s installed. OK" % v)
    return True
# tst_scipy()

def tst_networkx():
    print("Testing networkx..")

    try: import networkx
    except ImportError:
        print("  Not installed. NG")
        return False

    print("  %s installed. OK" % networkx.__version__)
    return True
# tst_networkx()

def tst_numpy():
    print("Testing NumPy..")

    try: import numpy
    except ImportError:
        print("  Not installed. NG")
        return False

    print("  %s installed. OK" % numpy.version.full_version)
    return True
# tst_numpy()

def tst_matplotlib():
    print("Testing Matplotlib..")

    try: import matplotlib
    except ImportError:
        print("  Not installed. NG")
        return False

    print("  %s installed. OK" % matplotlib.__version__)
    return True
# tst_matplotlib()

def tst_wx():
    print("Testing wxPython..")

    try: import wx
    except ImportError:
        print("  Not installed. NG")
        return False

    print("  %s installed. OK" % wx.version())
    return True
# tst_wx()

def show_env():
    import platform
    import sys
    print("Info:")
    print("   Python: %s" % platform.python_version())
    print("     Exec: %s" % sys.executable)
    print(" Platform: %s" % platform.platform())
# show_env()

def run():
    print("Testing installation of KAMO.")
    print("If you have trouble, please report the issue including all outputs:")
    show_env()
    print()

    failed = []

    for f in (tst_jsdir, tst_R, tst_numpy, tst_scipy, tst_networkx, tst_matplotlib, tst_wx,
              tst_xds, tst_xdsstat, tst_ccp4, tst_dials, tst_dials_module, tst_dxtbx_eiger, tst_adxv, tst_h5toxds):
        try:
            ret = f()
        except:
            print(traceback.format_exc())
            ret = False
            
        if not ret: failed.append(f.__name__)

    print()
    if not failed:
        print("All OK!")
        return True
    else:
        print("%d Failures (%s)" % (len(failed), ", ".join(failed)))
        return False
# run()

if __name__ == "__main__":
    run()

