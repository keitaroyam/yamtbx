"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

# Original code from iotbx/merging_statistics.py

import sys

from libtbx import slots_getstate_setstate
from iotbx import merging_statistics

class estimate_crude_resolution_cutoffs (slots_getstate_setstate) :
  """
  Uses incredibly simplistic criteria to determine the approximate
  resolution limit of the data based on the merging statistics (using much
  smaller bins than normal).  Not really appropriate for routine use, but
  useful for the pedantic and just-curious.
  """
  cutoffs_attr = ["i_over_sigma_cut", "r_merge_cut", "r_meas_cut",
    "completeness_cut_conservative", "completeness_cut_permissive",
    "cc_one_half_cut",]
  __slots__ = ["n_bins", "i_over_sigma_min", "r_merge_max", "r_meas_max",
    "completeness_min_conservative", "completeness_min_permissive",
    "cc_one_half_min", "d_min_overall",] + cutoffs_attr
  def __init__ (self,
      i_obs,
      crystal_symmetry=None,
      n_bins=None,
      i_over_sigma_min=2.0,
      r_merge_max=0.5,
      r_meas_max=0.5,
      completeness_min_conservative=0.9,
      completeness_min_permissive=0.5,
      cc_one_half_min=0.5) :
    self.n_bins = n_bins
    self.i_over_sigma_min = i_over_sigma_min
    self.r_merge_max = r_merge_max
    self.r_meas_max = r_meas_max
    self.completeness_min_conservative = completeness_min_conservative
    self.completeness_min_permissive = completeness_min_permissive
    self.cc_one_half_min = cc_one_half_min
    for attr in self.cutoffs_attr :
      setattr(self, attr, None)

    # Decide n_bins
    if n_bins is None:
        n_bins = min(200, int(i_obs.complete_set().indices().size()/500.+.5)) # not well tested.
        #print "n_bins=",n_bins

    i_obs.setup_binner(n_bins=n_bins)
    bins = []
    
    for bin in i_obs.binner().range_used():
      unmerged = i_obs.select(i_obs.binner().selection(bin))
      try:
        bin_stats = merging_statistics.merging_stats(unmerged,
                                                     anomalous=False)
        bins.append(bin_stats)
      except RuntimeError: # complains that no reflections left after sigma-filtering.
        continue

    self.d_min_overall = bins[-1].d_min
    d_min_last = float("inf")
    for bin in bins :
      if (self.i_over_sigma_cut is None) :
        if (bin.i_over_sigma_mean < self.i_over_sigma_min) :
          self.i_over_sigma_cut = d_min_last
      if (self.cc_one_half_cut is None) :
        if (bin.cc_one_half < self.cc_one_half_min) :
          self.cc_one_half_cut = d_min_last
      if (self.r_merge_cut is None) :
        if (bin.r_merge > self.r_merge_max) :
          self.r_merge_cut = d_min_last
      if (self.r_meas_cut is None) :
        if (bin.r_meas > self.r_meas_max) :
          self.r_meas_cut = d_min_last
      if (self.completeness_cut_conservative is None) :
        if (bin.completeness < completeness_min_conservative) :
          self.completeness_cut_conservative = d_min_last
      if (self.completeness_cut_permissive is None) :
        if (bin.completeness < completeness_min_permissive) :
          self.completeness_cut_permissive = d_min_last
      d_min_last = bin.d_min

  def show (self, out=sys.stdout, prefix="") :
    def format_d_min (value) :
      if (value is None) :
        return "(use all data)" #% self.d_min_overall
      return "%7.3f" % value
    print >> out, prefix + "Crude resolution cutoff estimates:"
    print >> out, prefix + "  resolution of all data          : %7.3f" % \
      self.d_min_overall
    if (self.cc_one_half_min is not None) :
      print >> out, prefix + "  based on CC(1/2) > %-6g       : %s" % \
        (self.cc_one_half_min, format_d_min(self.cc_one_half_cut))
    if (self.i_over_sigma_min is not None) :
      print >> out, prefix + "  based on mean(I/sigma) > %-6g : %s" % \
        (self.i_over_sigma_min, format_d_min(self.i_over_sigma_cut))
    if (self.r_merge_max is not None) :
      print >> out, prefix + "  based on R-merge < %-6g       : %s" % \
        (self.r_merge_max, format_d_min(self.r_merge_cut))
    if (self.r_meas_max is not None) :
      print >> out, prefix + "  based on R-meas < %-6g        : %s" % \
        (self.r_meas_max, format_d_min(self.r_meas_cut))
    if (self.completeness_min_conservative is not None) :
      print >> out, prefix + "  based on completeness > %-6g  : %s" % \
        (self.completeness_min_conservative,
         format_d_min(self.completeness_cut_conservative))
    if (self.completeness_min_permissive is not None) :
      print >> out, prefix + "  based on completeness > %-6g  : %s" % \
        (self.completeness_min_permissive,
         format_d_min(self.completeness_cut_permissive))

