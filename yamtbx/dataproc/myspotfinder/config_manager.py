from yamtbx.dataproc.XIO import XIO
from collections import OrderedDict

sp_params_strs =  OrderedDict(((("BL32XU", "EIGER9M", None, None), """\
distl {
  detector_tiles = 1
  peripheral_margin = 0
  minimum_spot_area = 2
  minimum_signal_height = 4.
  minimum_spot_height = None
}
xds {
  strong_pixel = 4
  minimum_number_of_pixels_in_a_spot = 3
  background_pixel = None
}
cheetah {
 ADCthresh = 5
 MinSNR = 8
 MinPixCount = 3
 MaxPixCount = 40
 LocalBGRadius = 2
 MinPeakSeparation = 0
 algorithm = 8
 binning = 1
}
software_binning = False
"""),

                               (("BL32XU", "MX225HS", "2x2", None), """\
distl {
  detector_tiles = 3
  peripheral_margin = 10
  minimum_spot_area = 5
  minimum_signal_height = 2.
  minimum_spot_height = None
}
xds {
  strong_pixel = 4
  minimum_number_of_pixels_in_a_spot = 3
  background_pixel = None
}
software_binning = False
"""),
                              (("BL32XU", "MX225HS", "4x4", None), """\
distl {
  detector_tiles = 3
  peripheral_margin = 10
  minimum_spot_area = 2
  minimum_signal_height = 4.
  minimum_spot_height = None
}
xds {
  strong_pixel = 4
  minimum_number_of_pixels_in_a_spot = 3
  background_pixel = None
}
software_binning = False
"""),
                              (("BL32XU", "MX225HS", "8x8", None), """\
distl {
  detector_tiles = 3
  peripheral_margin = 5
  minimum_spot_area = 1
  minimum_signal_height = 4.
  minimum_spot_height = None
}
xds {
  strong_pixel = 4
  minimum_number_of_pixels_in_a_spot = 3
  background_pixel = None
}
software_binning = False
"""),
                              (("BL32XU", "MX225HE", "2x2", None), """\
distl {
  detector_tiles = 3
  peripheral_margin = 10
  minimum_spot_area = 5
  minimum_signal_height = 2.
  minimum_spot_height = None
}
xds {
  strong_pixel = 4
  minimum_number_of_pixels_in_a_spot = 3
  background_pixel = None
}
software_binning = False
"""),
                              (("BL32XU", "MX225HE", "3x3", None), """\
distl {
  detector_tiles = 3
  peripheral_margin = 10
  minimum_spot_area = 4
  minimum_signal_height = 3.
  minimum_spot_height = None
}
xds {
  strong_pixel = 4
  minimum_number_of_pixels_in_a_spot = 3
  background_pixel = None
}
software_binning = False
"""),
                              (("BL32XU", "MX225HE", "4x4", None), """\
distl {
  detector_tiles = 3
  peripheral_margin = 10
  minimum_spot_area = 2
  minimum_signal_height = 4.
  minimum_spot_height = None
}
xds {
  strong_pixel = 4
  minimum_number_of_pixels_in_a_spot = 3
  background_pixel = None
}
software_binning = False
"""),
                              (("BL32XU", "MX225HE", "8x8", None), """\
distl {
  detector_tiles = 3
  peripheral_margin = 5
  minimum_spot_area = 1
  minimum_signal_height = 4.
  minimum_spot_height = None
}
xds {
  strong_pixel = 4
  minimum_number_of_pixels_in_a_spot = 3
  background_pixel = None
}
software_binning = False
"""),
                              (("BL32XU", "MX225HE", "16x16", None), """\
distl {
  detector_tiles = 3
  peripheral_margin = 5
  minimum_spot_area = 1
  minimum_signal_height = 4.
  minimum_spot_height = None
}
xds {
  strong_pixel = 4
  minimum_number_of_pixels_in_a_spot = 3
  background_pixel = None
}
software_binning = False
"""),
                              (("BL32XU", "Q315r", "2x2", None), """\
distl {
  detector_tiles = 3
  peripheral_margin = 10
  minimum_spot_area = 5
  minimum_signal_height = 2.
  minimum_spot_height = None
}
xds {
  strong_pixel = 4
  minimum_number_of_pixels_in_a_spot = 3
  background_pixel = None
}
software_binning = False
"""),
                              (("BL32XU", "CMOS", "1x1", None), """\
distl {
  detector_tiles = 1
  peripheral_margin = 0
  minimum_spot_area = 3
  minimum_signal_height = 2.
  minimum_spot_height = None
}
xds {
  strong_pixel = 4
  minimum_number_of_pixels_in_a_spot = 3
  background_pixel = None
}
software_binning = False
"""),
                              (("BL41XU", "PILATUS3 6M", None, None), """\
distl {
  detector_tiles = None
  peripheral_margin = 10
  minimum_spot_area = 2
  minimum_signal_height = 4
  minimum_spot_height = None
}
xds {
  strong_pixel = 4
  minimum_number_of_pixels_in_a_spot = 3
  background_pixel = None
}
software_binning = False
""")
                              ))

def get_common_params_str(use_cuda=False, env="oys"):
    if use_cuda:
        dic = dict(use_cuda="True")
    else:
        #if nproc is None: nproc = get_number_of_processors(default=4)
        #if env == "ppu": nproc //= 2
        dic = dict(use_cuda="False")

    return """\
engine = *distl xds
distl {
  res {
    outer = 5.
    inner = 30.
  }
  scanbox_windows = 101 51 51
}
xds {
  do_defpix = True
  value_range_for_trusted_detector_pixels = 9000. 30000
}
cuda_median_background {
  active = %(use_cuda)s
  filter_radius = 10
  filter_repeat = 1
}

#bkg_image = /home/yam/work/smoothing_131114/xds_process_scan/BKGINIT.cbf
#gain_image = /home/yam/work/smoothing_131114/xds_process_scan/GAIN.cbf
#bkg_image = /home/yam/work/smoothing_131114/my_scan172/honki/bkginit_20_20_1.cbf
#gain_image = /home/yam/work/smoothing_131114/my_scan172/honki/test_rev_median5.cbf
#gain_image_nbxy = 3,3
""" % dic

def get_key_by_img(imgfile):
    im = XIO.Image(imgfile)
    if im.header["ImageType"] == "marccd":
        if im.header["SerialNumber"] in ("106", None): # None for 2013B
            if im.header["Height"] == im.header["Width"] == 1440:
                return ("BL32XU", "MX225HS", "4x4", None)
            if im.header["Height"] == im.header["Width"] == 2880:
                return ("BL32XU", "MX225HS", "2x2", None)
            if im.header["Height"] == im.header["Width"] == 720:
                return ("BL32XU", "MX225HS", "8x8", None)

        if im.header["SerialNumber"] == "31":
            if im.header["Height"] == im.header["Width"] == 384:
                return ("BL32XU", "MX225HE", "16x16", None)
            if im.header["Height"] == im.header["Width"] == 768:
                return ("BL32XU", "MX225HE", "8x8", None)
            if im.header["Height"] == im.header["Width"] == 1536:
                return ("BL32XU", "MX225HE", "4x4", None)
            if im.header["Height"] == im.header["Width"] == 2046:
                return ("BL32XU", "MX225HE", "3x3", None)
            if im.header["Height"] == im.header["Width"] == 3072:
                return ("BL32XU", "MX225HE", "2x2", None)
    elif im.header["ImageType"] == "adsc":
        if im.header["Height"]==im.header["Width"]==2352 and int(im.header["PixelX"]*1000)==50:
            return ("BL32XU", "CMOS", "1x1", None) # This may be used at BL26B2.

        if im.header["SerialNumber"] == "915":
            if im.header["Height"] == im.header["Width"] == 3072:
                return ("BL32XU", "Q315r", "2x2", None)

    elif im.header["SerialNumber"] == "PILATUS3 6M, S/N 60-0125":
        return ("BL41XU", "PILATUS3 6M", None, None)

    raise Exception("We do not know such a detector")
# get_key_by_img()

