"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
generated_by_XYCORR = ("XYCORR.LP", "X-CORRECTIONS.cbf", "Y-CORRECTIONS.cbf")

generated_by_INIT = ("INIT.LP", "BKGINIT.cbf", "BLANK.cbf", "GAIN.cbf")

generated_by_COLSPOT = ("COLSPOT.LP", "SPOT.XDS")

generated_by_IDXREF = ("IDXREF.LP", "SPOT.XDS", "XPARM.XDS")

generated_by_DEFPIX = ("DEFPIX.LP", "ABS.cbf", "BKGPIX.cbf")

generated_by_INTEGRATE = ("INTEGRATE.HKL", "INTEGRATE.LP", "FRAME.cbf")

generated_by_CORRECT = ("GXPARM.XDS", "XDS_ASCII.HKL", "CORRECT.LP",
                        "GX-CORRECTIONS.cbf", "GY-CORRECTIONS.cbf",
                        "DX-CORRECTIONS.cbf", "DY-CORRECTIONS.cbf",
                        "MODPIX.cbf", "DECAY.cbf", "ABSORP.cbf",
                        )
generated_by_xdsstat = ("XDSSTAT.LP", "anom.pck", "corr.pck", "misfits.pck",
                        "nobs.pck", "peaks.pck", "rf.pck", "rlps.pck", "scales.pck",
                        )

generated_after_DEFPIX = generated_by_DEFPIX + generated_by_INTEGRATE + generated_by_CORRECT + generated_by_xdsstat

needed_by_COLSPOT = ("BKGINIT.cbf", "BLANK.cbf", "GAIN.cbf", "X-CORRECTIONS.cbf", "Y-CORRECTIONS.cbf")

needed_by_INTEGRATE = ("XPARM.XDS", "BKGPIX.cbf", "BLANK.cbf", "GAIN.cbf", "X-CORRECTIONS.cbf", "Y-CORRECTIONS.cbf")

