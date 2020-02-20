# Class for reading .emd file by Velox.
# This code was written based on FormatSER.py from https://github.com/cctbx/dxtbx/blob/master/format/FormatSER.py
from __future__ import absolute_import, division, print_function

import struct
import h5py
import numpy
import os
import json

from scitbx.array_family import flex
from dxtbx.format.Format import Format
from dxtbx.format.FormatMultiImage import FormatMultiImage

def get_metadata(metadata):
    mds = []
    for i in xrange(metadata.shape[1]):
        metadata_array = metadata[:, i].T
        mdata_string = metadata_array.tostring().decode("utf-8")
        mds.append(json.loads(mdata_string.rstrip('\x00')))

    return mds
# get_metadata()

def analyse_angle(metadata):
    alphas = []
    
    for i, md in enumerate(metadata):
        alpha = numpy.rad2deg(float(md["Stage"]["AlphaTilt"]))
        alphas.append(alpha)

    d_alphas = numpy.diff(alphas)
    q25, q50, q75 = numpy.percentile(d_alphas, [25, 50, 75])
    iqr = q75-q25
    iqrc = 1.5
    lowlim, highlim = q25 - iqrc*iqr, q75 + iqrc*iqr
    d_alphas2 = d_alphas[numpy.where(numpy.logical_and(d_alphas>lowlim, d_alphas<highlim))] # outlier rejected
    d_alpha_z = abs(d_alphas-numpy.mean(d_alphas2))/numpy.std(d_alphas2)

    valid_range = [0, len(metadata)-1]
    for i in xrange(len(metadata)-1):
        if d_alpha_z[i] < 3: break
        valid_range[0] = i+1

    for i in reversed(xrange(len(metadata)-1)):
        if d_alpha_z[i] < 3: break
        valid_range[1] = i

    if valid_range[0] > valid_range[1]:
        valid_range = [0, len(metadata)-1] # reset
        
    mean_alpha_step = (alphas[valid_range[1]] - alphas[valid_range[0]])/(valid_range[1]-valid_range[0])

    return valid_range, mean_alpha_step
# analyse_angle()


class FormatEMD(FormatMultiImage, Format):
    def __init__(self, image_file, **kwargs):

        from dxtbx import IncorrectFormatError

        if not self.understand(image_file):
            raise IncorrectFormatError(self, image_file)
        FormatMultiImage.__init__(self, **kwargs)
        Format.__init__(self, image_file, **kwargs)

    @staticmethod
    def understand(image_file):

        try:
            h = h5py.File(image_file, "r")
        except IOError:
            return False

        if not "/Data/Image" in h:
            return False

        keys = h["/Data/Image"].keys()
        if len(keys) > 1: return False
        d = h["/Data/Image"][keys[0]]
        if "Data" in d and "Metadata" in d and len(d["Data"].shape) == 3:
            return True

        return False

    @staticmethod
    def _read_metadata(image_file):
        h = h5py.File(image_file, "r")
        ret = {}

        image_path = h["/Data/Image"]
        assert len(image_path.keys()) == 1
        k = image_path.keys()[0]

        ret["image_file"] = image_file
        ret["file_handle"] = h
        ret["data_path"] = "/Data/Image/%s/Data" % k
        ret["metadata_path"] = "/Data/Image/%s/Metadata" % k

        metadata = get_metadata(h[ret["metadata_path"]])
        valid_range, mean_alpha_step = analyse_angle(metadata)
        data = h[ret["data_path"]]
        
        ret["n_frames"] = data.shape[2]
        ret["valid_range"] = valid_range
        ret["mean_alpha_step"] = mean_alpha_step
        ret["width"], ret["height"] = data.shape[:2]
        ret["binning"] = int(metadata[0]["BinaryResult"]["ImageSize"]["width"])//ret["width"]
        
        return ret

    def _start(self):
        """Open the image file, read useful metadata into an internal dictionary
        self._header_dictionary"""

        self._header_dictionary = self._read_metadata(self._image_file)

        return

    def _goniometer(self):
        """Dummy goniometer, 'vertical' as the images are viewed. Not completely
        sure about the handedness yet"""

        if self._header_dictionary["mean_alpha_step"] > 0: # XXX is this really OK??
            return self._goniometer_factory.known_axis((0, -1, 0))
        else:
            return self._goniometer_factory.known_axis((0, 1, 0))

    def _detector(self):
        """Dummy detector"""

        image_size = (self._header_dictionary["width"], self._header_dictionary["height"])

        binning = self._header_dictionary["binning"]
        pixel_size = 0.014*binning, 0.014*binning

        distance = 2000
        trusted_range = (-4, 65535)
        beam_centre = [(p * i) / 2 for p, i in zip(pixel_size, image_size)]
        d = self._detector_factory.simple(
            "PAD",
            distance,
            beam_centre,
            "+x",
            "-y",
            pixel_size,
            image_size,
            trusted_range,
        )
        # Not sure what the gain is
        # for p in d: p.set_gain(8)
        return d

    def _beam(self):
        """Dummy unpolarized beam, energy 200 keV"""

        wavelength = 0.02508 # XXX
        return self._beam_factory.make_polarized_beam(
            sample_to_source=(0.0, 0.0, 1.0),
            wavelength=wavelength,
            polarization=(0, 1, 0),
            polarization_fraction=0.5,
        )

    def _scan(self):
        """Dummy scan for this stack"""

        vr = self._header_dictionary["valid_range"]
        image_range = (vr[0]+1, vr[1]+1)
        print("Recommended image_raneg=", image_range)
        image_range = (1, self._header_dictionary["n_frames"])
        exposure_times = 0.0
        nframes = self._header_dictionary["n_frames"] #vr[1]-vr[0]+1
        #nframes = vr[1]-vr[0]+1
        osc_step = abs(self._header_dictionary["mean_alpha_step"])
        oscillation = (osc_step*(vr[0]-1), osc_step)

        # FIXME we do actually have acquisition times, might as well use them
        epochs = [0] * nframes

        #print(len(epochs), self.get_num_images())
        

        return self._scan_factory.make_scan(
            image_range, exposure_times, oscillation, epochs, deg=True
        )


    def get_num_images(self):
        #h = self._header_dictionary["file_handle"]
        #data_path = self._header_dictionary["data_path"]
        #return h[data_path].shape[2]
        #vr = self._header_dictionary["valid_range"]
        return self._header_dictionary["n_frames"] # vr[1] - vr[0] + 1
        #return vr[1] - vr[0] + 1

    # This is still required for dials_regression/test.py
    def get_detectorbase(self):
        pass

    def get_goniometer(self, index=None):
        return Format.get_goniometer(self)

    def get_detector(self, index=None):
        return Format.get_detector(self)

    def get_beam(self, index=None):
        return Format.get_beam(self)

    def get_scan(self, index=None):
        if index == None:
            return Format.get_scan(self)
        else:
            scan = Format.get_scan(self)
            return scan[index]

    def get_image_file(self, index=None):
        return Format.get_image_file(self)

    def get_raw_data(self, index):
        #print(self._header_dictionary["valid_range"][0])
        #index += self._header_dictionary["valid_range"][0]
        h = h5py.File(self._header_dictionary["image_file"], "r")
        data_path = self._header_dictionary["data_path"]
        raw_data = h[data_path][:,:,index].astype(numpy.int32) # flex.int does not suppert int16

        offset_key = "DXTBX_EMD_OFFSET"
        if os.environ.get(offset_key):
            print("DEBUG:: adding %s for %d"%(os.environ[offset_key], index))
            raw_data += int(os.environ[offset_key])

        return flex.int(raw_data)
