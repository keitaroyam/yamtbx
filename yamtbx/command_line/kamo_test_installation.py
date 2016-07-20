# LIBTBX_SET_DISPATCHER_NAME kamo.test_installation

"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

from yamtbx import util
import os

def tst_jsdir():
    print "Testing location.."

    import libtbx.load_env
    d3path = libtbx.env.find_in_repositories("yamtbx/dataproc/auto/js/d3-3.5.10")
    if not d3path:
        print "  Can't find d3-3.5.10 directory. Please check location of yamtbx. NG"
        return False
    
    print "  %s. OK" % libtbx.env.find_in_repositories("yamtbx")
    return True
# tst_jsdir()

def tst_R():
    print "Testing R.."

    rcode, out, err = util.call("Rscript", '-e "print(cos(pi))"')
    if rcode != 0 or out.strip() != '[1] -1':
        print "  Rscript is not avaiable. NG"
        return False
    
    rcode, out, err = util.call("Rscript", '-e "library(rjson)"')
    if rcode != 0:
        print "  rjson is not installed. NG"
        return False
    
    print "  OK"
    return True
# tst_R()

def tst_xds():
    print "Testing XDS.."
    rcode, out, err = util.call("xds_par")
    if rcode != 0:
        print "  Not installed. NG"
        return False

    if "license expired" in out:
        print "  license expired. Get latest version. NG"
        return False

    print "  OK"
    return True
# tst_xds()

def tst_ccp4():
    print "Testing ccp4.."
    if "CCP4" not in os.environ or not os.path.isdir(os.environ["CCP4"]):
        print "  Not installed. NG"
        return False
    
    if not os.path.isfile(os.path.join(os.environ["CCP4"], "share/blend/R/blend0.R")):
        print "  BLEND is not available. NG"
        return False

    print "  OK"
    return True
# tst_ccp4()

#def tst_h5():
#    print "Testing hdf5..",   

def tst_scipy():
    print "Testing scipy.."

    try: import scipy.optimize
    except ImportError:
        print "  Not installed. NG"
        return False

    try: scipy.optimize.least_squares
    except AttributeError:
        print "  scipy.optimize.least_squares is not available. Update the version. NG"
        return False

    print "  %s installed. OK" % scipy.version.full_version
    return True
# tst_scipy()

def tst_networkx():
    print "Testing networkx.."

    try: import networkx
    except ImportError:
        print "  Not installed. NG"
        return False

    print "  %s installed. OK" % networkx.__version__
    return True
# tst_networkx()

def tst_numpy():
    print "Testing NumPy.."

    try: import numpy
    except ImportError:
        print "  Not installed. NG"
        return False

    print "  %s installed. OK" % numpy.version.full_version
    return True
# tst_numpy()

def tst_matplotlib():
    print "Testing Matplotlib.."

    try: import matplotlib
    except ImportError:
        print "  Not installed. NG"
        return False

    print "  %s installed. OK" % matplotlib.__version__
    return True
# tst_matplotlib()

def tst_wx():
    print "Testing wxPython.."

    try: import wx
    except ImportError:
        print "  Not installed. NG"
        return False

    print "  %s installed. OK" % wx.version()
    return True
# tst_wx()

def show_env():
    import platform
    import sys
    print "Info:"
    print "   Python: %s" % platform.python_version()
    print "     Exec: %s" % sys.executable
    print " Platform: %s" % platform.platform()
# show_env()

def run():
    print "Testing installation of KAMO."
    print "If you have trouble, please report the issue including all outputs:"
    show_env()
    print

    failed = []

    for f in (tst_jsdir, tst_R, tst_xds, tst_ccp4, tst_numpy,
              tst_scipy, tst_networkx, tst_matplotlib, tst_wx):
        ret = f()
        if not ret: failed.append(f.func_name)

    print
    if not failed:
        print "All OK!"
    else:
        print "%d Failures (%s)" % (len(failed), ", ".join(failed))
# run()

if __name__ == "__main__":
    run()

