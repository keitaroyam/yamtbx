import h5py
import datetime

raw = lambda x: x

def date_seconds(datestr):
    # '2016-02-25T20:00:36.653853'
    a = datetime.datetime.strptime(datestr, "%Y-%m-%dT%H:%M:%S.%f")
    return int(a.strftime("%s"))
# 

class Interpreter:
    HTD = {
    'ExposureTime':(['Exposure_time'], raw),
    'BeamX':(['BeamX', 'PixelX'], lambda x,y:x*y*1.e3),
    'BeamY':(['BeamY', 'PixelY'], lambda x,y:x*y*1.e3),
    'Distance':(['Distance'], lambda x: x*1.e3),
    'Wavelength':(['Wavelength'], raw),
    'PixelX':(['PixelX'], lambda x: x*1.e3),
    'PixelY':(['PixelY'], lambda x: x*1.e3),
    'Width':(['Width'], raw),
    'Height':(['Height'], raw),
    'PhiStart':(['PhiStart'], raw),
    'PhiEnd':(['PhiEnd'], raw),
    'PhiWidth':(['PhiWidth'], raw),
    'SensorThickness':(['SensorThickness'], lambda x:x*1.e3),
    'SensorMaterial':(['SensorMaterial'], str),
    'Nimages_each':(['Nimages_each'], raw),
    'Nimages':(['Nimages'], raw),
    'Overload':(['Overload'], raw),
    #'EdgeResolution':(['Pixel_size','Binary-Size-Second-Dimension','Detector_distance','Wavelength'], \
    #                 get_edge_resolution),
    # Added keys from Graeme's convention.
    #'TwoTheta':(['Detector_2theta'], FLOAT1),   # No example yet...
    'SerialNumber':(['SerialNumber'], str),
    #'HeaderSize':(['HEADER_SIZE'], int),
    #'OscAxis':(['Oscillation_axis'], lambda x: x.split(",")[0].lower().strip()),
    'DateStr':(['DateStr'], str),
    'DateSeconds':(['DateStr'], date_seconds),
    }

    def getRawHeadDict(self, filename):
        self.raw_head_dict = {}

        h5 = h5py.File(filename, "r")

        beam = h5["entry/instrument/beam"]
        detector = h5["entry/instrument/detector"]
        detectorsp = h5["entry/instrument/detector/detectorSpecific"]
        gonio = h5["entry/sample/goniometer"]

        self.raw_head_dict["BeamX"] = detector["beam_center_x"].value # in pixel
        self.raw_head_dict["BeamY"] = detector["beam_center_y"].value # in pixel
        self.raw_head_dict["Distance"] = detector["detector_distance"].value # in m
        self.raw_head_dict["SerialNumber"] = detector["detector_number"].value # S/N
        self.raw_head_dict["SensorThickness"] = detector["sensor_thickness"].value # in m
        self.raw_head_dict["SensorMaterial"] = detector["sensor_material"].value
        self.raw_head_dict["PixelX"] = detector["x_pixel_size"].value
        self.raw_head_dict["PixelY"] = detector["y_pixel_size"].value
        self.raw_head_dict["ExposureTime"] = detector["count_time"].value # in sec
        self.raw_head_dict["Wavelength"] = beam["incident_wavelength"].value # in A
        self.raw_head_dict["Width"] = detectorsp["x_pixels_in_detector"].value #NX
        self.raw_head_dict["Height"] = detectorsp["y_pixels_in_detector"].value
        self.raw_head_dict["Overload"] = detectorsp["countrate_correction_count_cutoff"].value
        self.raw_head_dict["DateStr"] = detectorsp["data_collection_date"].value
        self.raw_head_dict["PhiStart"] = gonio["omega"][0]
        self.raw_head_dict["PhiEnd"] = gonio["omega_end"][-1]
        self.raw_head_dict["PhiWidth"] = gonio["omega_range_average"].value

        self.raw_head_dict["Nimages_each"] = detectorsp["nimages"].value
        self.raw_head_dict["Nimages"] = 0
        for k in sorted(h5["entry/data"].keys()):
            try:
                self.raw_head_dict["Nimages"] += h5["entry/data"][k].shape[0]
            except KeyError:
                break

        return self.raw_head_dict

