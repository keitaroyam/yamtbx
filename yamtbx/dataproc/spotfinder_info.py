"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals
import collections

from cctbx.array_family import flex

all_keys = ('all',
            'spots_non-ice',
            'hi_pass_resolution_spots',
            'ice_free_resolution_spots',
            'lo_pass_resolution_spots',
            'spots_unimodal',
            'spots_low_skew',
            'spots_good_intensity',
            'spots_separated')

def filter_spots_with_ex_resolution_range(spots, resolutions, exclude_resolution_ranges):
    # TODO: more efficient way!
    ret = []
    for spot, res in zip(spots, resolutions):
        found = False
        for rng in exclude_resolution_ranges:
            if min(rng) < res < max(rng):
                found = True
                break
        if not found:
            ret.append(spot)
    return ret
 # filter_spots_with_ex_resolution_range()

class SpotsBase(object):
    def __init__(self):
        self.all_spots = None # flex.distl_spot type
        self.spots = collections.OrderedDict() # {descriptor: flex.int}
        self.intensities = {}
        self.resolutions = {}
        self.total_integrated_signal = {}
        self.mean_integrated_signal= {}
        self.median_integrated_signal= {}
        self.n_spots = {}
    # __init__()

    def keys(self):
        return ("all",) + tuple(self.spots.keys())
    # keys()

    def get_spots(self, key, exclude_resolution_ranges=[]):
        key = self.find_nearest_key(key)
        spots = self.all_spots if key == "all" else [self.all_spots[i] for i in self.spots[key]]

        if exclude_resolution_ranges == []:
            return spots
        else:
            return filter_spots_with_ex_resolution_range(spots, self.resolutions[key], exclude_resolution_ranges)
    # get_spots()

    def find_nearest_key(self, key):
        assert key in all_keys or key == "xds"

        if key in list(self.spots.keys()) or key == "all": # fast ver
            return key

        for key in all_keys[all_keys.index(key)::-1]:
            if key in list(self.spots.keys()):
                return key

        return "all"
    # find_nearest_key()

    def get_total_integrated_signal(self, key, exclude_resolution_ranges=[]):
        key = self.find_nearest_key(key)

        if exclude_resolution_ranges == []:
            return self.total_integrated_signal[key]
        else:
            intensities = filter_spots_with_ex_resolution_range(self.intensities[key],
                                                                self.resolutions[key],
                                                                exclude_resolution_ranges)
            return flex.sum(flex.double(intensities))
    # get_total_integrated_signal()

    def get_mean_integrated_signal(self, key, exclude_resolution_ranges=[]):
        key = self.find_nearest_key(key)
        if exclude_resolution_ranges == []:
            return self.mean_integrated_signal[key]
        else:
            intensities = filter_spots_with_ex_resolution_range(self.intensities[key],
                                                                self.resolutions[key],
                                                                exclude_resolution_ranges)

            if len(intensities) > 0:
                return flex.sum(flex.double(intensities)) / len(intensities)
            else:
                return 0
    # get_mean_integrated_signal()

    def get_median_integrated_signal(self, key, exclude_resolution_ranges=[]):
        key = self.find_nearest_key(key)
        if exclude_resolution_ranges == []:
            if hasattr(self, "median_integrated_signal"): # backward compatibility
                return self.median_integrated_signal[key]
            else:
                print("NOTE: median is not available. using mean value.")
                return self.mean_integrated_signal[key]
        else:
            intensities = filter_spots_with_ex_resolution_range(self.intensities[key],
                                                                self.resolutions[key],
                                                                exclude_resolution_ranges)

            if len(intensities) > 0:
                return flex.median(flex.double(intensities))
            else:
                return 0
    # get_median_integrated_signal()

    def get_n_spots(self, key, exclude_resolution_ranges=[]):
        if exclude_resolution_ranges == []:
            return self.n_spots[self.find_nearest_key(key)]
        else:
            return len(self.get_spots(key, exclude_resolution_ranges))
    # get_n_spots()

    def get_resolutions(self, key):
        return self.resolutions[self.find_nearest_key(key)]
    # get_resolutions()

    def get_summed_peaks(self, key):
        pass

class DummySpot(object): # Dummy spot class for XDS result
    def __init__(self, x, y, d, intensity):
        self.x = x
        self.y = y
        self.resolution = d
        self.intensity = intensity
    # __init__()
    def max_pxl_x(self): return self.y # Swapped!
    def max_pxl_y(self): return self.x
# class DummySpot

class XDSSpots(SpotsBase):

    def __init__(self, spots):
        SpotsBase.__init__(self)
        self.set_spots(spots)

    def set_spots(self, spots):
        self.all_spots = []

        for x, y, d, intensity in spots:
            self.all_spots.append(DummySpot(x, y, d, intensity))

        self.spots["xds"] = list(range(len(self.all_spots)))

        for k in list(self.keys()):
            summed_wts = [spot.intensity for spot in self.get_spots(k)]
            self.intensities[k] = summed_wts
            self.resolutions[k] = [spot.resolution for spot in self.get_spots(k)] # XXX calculate resolution!!

            total_summed = sum(summed_wts)

            if len(summed_wts) > 0:
                self.mean_integrated_signal[k] = total_summed / len(summed_wts)
                self.median_integrated_signal[k] = flex.median(flex.double(summed_wts))
            else:
                self.mean_integrated_signal[k] = 0.
                self.median_integrated_signal[k] = 0.

            self.total_integrated_signal[k] = total_summed
            self.n_spots[k] = len(summed_wts)
    # set_spots()


class DistlSpots(SpotsBase):
    def __init__(self, fstats):
        SpotsBase.__init__(self)
        """
        fstats is the result of practical_heuristics.heuristics_base.oneImage()
        """
        self.set_fstats(fstats)
    # __init__()

    def set_fstats(self, fstats):
        self.all_spots = None
        self.spots = collections.OrderedDict()
        self.total_integrated_signal = {}
        self.mean_integrated_signal= {}
        self.median_integrated_signal= {}
        self.n_spots = {}

        for k in sorted(fstats.nodes.keys()):
            node = fstats.nodes[k] # XXX some data in node.data will be corrupted (e.g. resolution, wts) but x and y coordinates look ok (it works later). why?
            if node.descriptor == "spots_total":
                self.all_spots = node.data
            else:
                self.spots[node.descriptor] = node.data

        # Pre-calculate stats
        for k in list(self.keys()):
            summed_wts = [flex.sum(spot.wts) for spot in self.get_spots(k)]
            self.intensities[k] = summed_wts
            self.resolutions[k] = [spot.resolution for spot in self.get_spots(k)]

            total_summed = flex.sum(flex.double(summed_wts))

            if len(summed_wts) > 0:
                self.mean_integrated_signal[k] = total_summed / len(summed_wts)
                self.median_integrated_signal[k] = flex.median(flex.double(summed_wts))
            else:
                self.mean_integrated_signal[k] = 0.
                self.median_integrated_signal[k] = 0.

            self.total_integrated_signal[k] = total_summed
            self.n_spots[k] = len(summed_wts)
    # set_fstats()

    def get_summed_peaks(self, key):
        return flex.double([flex.sum(spot.wts) for spot in self.get_spots(key)])
    # get_total_integrated_signal()
