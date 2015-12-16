"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
import struct

timeconv = lambda t: datetime.datetime.strptime(t, "%m%d%H%M%Y.%S")


class MarCCDHeader:
    def __init__(self):
        pass
    # __init__()
# class MarCCDHeader


class MarCCD:
    def __init__(self, img_in):
        self.img_in = img_in
    # __init__()

    def read_header(self):
        f = open(self.img_in, "rb")
        f.seek(1024)#+256+128+256+128+128+128+128)

        # typedef struct frame_header_type
        # /* File/header format parameters (256 bytes) */
        fmt = "I16s39I"
        bin = f.read(struct.calcsize(fmt))
        data = struct.unpack(fmt, bin)
        keys = "header_type", "header_name", "header_major_version", "header_minor_version", "header_byte_order", "data_byte_order", "header_size", "frame_type", "magic_number", "compression_type", "compression1", "compression2", "compression3", "compression4", "compression5", "compression6", "nheaders", "nfast", "nslow", "depth", "record_length", "signif_bits", "data_type", "saturated_value", "sequence", "nimages", "origin", "orientation", "view_direction", "overflow_location", "over_8_bits", "over_16_bits", "multiplexed", "nfastimages", "nslowimages", "background_applied", "bias_applied", "flatfield_applied", "distortion_applied", "original_header_type", "file_saved",  
        for k,d in zip(keys, data):
            print k, d

        f.seek(256 - struct.calcsize(fmt), 1)

        # /* Data statistics (128) */
        fmt = "3L8I" #+ "%ds" % ((64-40)*4-16)
        bin = f.read(struct.calcsize(fmt))
        data = struct.unpack(fmt, bin)
        keys = "total_counts", "special_counts1", "special_counts2", "min", "max", "mean", "rms", "p10", "p90", "stats_uptodate",  
        for k,d in zip(keys, data):
            print k, d

        f.seek(128 - struct.calcsize(fmt), 1)

        # /* More statistics (256) */
        fmt = "128H"
        assert struct.calcsize(fmt) == 256
        bin = f.read(struct.calcsize(fmt))
        data = struct.unpack(fmt, bin)
        print "percentile", data

        # /* Goniostat parameters (128 bytes) */
        fmt = "28I"
        bin = f.read(struct.calcsize(fmt))
        data = struct.unpack(fmt, bin)
        keys = "xtal_to_detector", "beam_x", "beam_y", "integration_time", "exposure_time", "readout_time", "nreads", "start_twotheta", "start_omega", "start_chi", "start_kappa", "start_phi", "start_delta", "start_gamma", "start_xtal_to_detector", "end_twotheta", "end_omega", "end_chi", "end_kappa", "end_phi", "end_delta", "end_gamma", "end_xtal_to_detector", "rotation_axis", "rotation_range", "detector_rotx", "detector_roty", "detector_rotz"
        for k,d in zip(keys, data):
            print k, d
        
        f.seek(128 - struct.calcsize(fmt), 1)

        # /* Detector parameters (128 bytes) */
        fmt = "5I"
        bin = f.read(struct.calcsize(fmt))
        data = struct.unpack(fmt, bin)
        keys = "detector_type", "pixelsize_x", "pixelsize_y", "mean_bias", "photons_per_100adu"
        for k,d in zip(keys, data):
            print k, d
        
        f.seek(128 - struct.calcsize(fmt), 1)

        # /* X-ray source and optics parameters (128 bytes) */
        fmt = "10I%ds10I%ds%ds" % (4*4, 4*4, (32-28)*4)
        bin = f.read(struct.calcsize(fmt))
        data = struct.unpack(fmt, bin)
        keys = "source_type", "source_dx", "source_dy", "source_wavelength", "source_power", "source_voltage", "source_current", "source_bias", "source_polarization_x", "source_polarization_y", "reserve_source", "optics_type", "optics_dx", "optics_dy", "optics_wavelength", "optics_dispersion", "optics_crossfire_x", "optics_crossfire_y", "optics_angle", "optics_polarization_x", "optics_polarization_y", "reserve_optics", "reserve5"
        for k,d in zip(keys, data):
            print k, d

        # /* File parameters (1024 bytes) */
        fmt = "128s128s64s32s32s32s512s"
        bin = f.read(struct.calcsize(fmt))
        data = struct.unpack(fmt, bin)
        keys = "filetitle", "filepath", "filename", "acquire_timestamp", "header_timestamp", "save_timestamp", "file_comments",  
        for k,d in zip(keys, data):
            print k, d

        f.seek(1024 - struct.calcsize(fmt), 1)

        # /* Dataset parameters (512 bytes) */
        fmt = "512s"
        bin = f.read(struct.calcsize(fmt))
        data = struct.unpack(fmt, bin)
        keys = "dataset_comments",
        for k,d in zip(keys, data):
            print k, d

    # read_header()
# class MarCCD
