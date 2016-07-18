#ifndef YAMTBX_DATAPROC_CRYSTFEL_MERGE_EQUIVALENTS_H
#define YAMTBX_DATAPROC_CRYSTFEL_MERGE_EQUIVALENTS_H

#include <cctbx/sgtbx/space_group.h>
#include <scitbx/array_family/ref.h>
#include <boost/unordered_map.hpp>

namespace af = scitbx::af;
namespace cctbx {
  namespace miller {
    size_t hash_value(const index<>& d)
    {
      size_t h = 0;
      boost::hash_combine(h, d[0]);
      boost::hash_combine(h, d[1]);
      boost::hash_combine(h, d[2]);
      return h;
    }
  }
}


template <typename FloatType=double>
struct merging_data{
  merging_data() {}
  merging_data(FloatType iob, FloatType w, FloatType m, int r) 
  : iobs(iob),weight(w),m2(m),red(r) {}
  FloatType iobs, weight, m2;
  int red;
};

template <typename FloatType=double>
struct merge_equivalents_crystfel
{
  af::shared<cctbx::miller::index<> > indices;
  af::shared<FloatType> data;
  af::shared<FloatType> sigmas;
  af::shared<int> redundancies;

  void add_observations(af::const_ref<cctbx::miller::index<> > const& unmerged_indices,
			af::const_ref<FloatType> const& unmerged_data)
  {
    const FloatType w = 1.;
    if (unmerged_indices.size() == 0) return;

    for(int i=0, i_end=unmerged_indices.size(); i<i_end; ++i) {
      const FloatType &iobs = unmerged_data[i];
      const cctbx::miller::index<> &hkl = unmerged_indices[i];
      add_observation(iobs, w, hkl);
    }
  }

  void add_observations(af::const_ref<cctbx::miller::index<> > const& unmerged_indices,
			af::const_ref<FloatType> const& unmerged_data,
			af::const_ref<FloatType> const& unmerged_sigmas)
  {
    if (unmerged_indices.size() == 0) return;

    for(int i=0, i_end=unmerged_indices.size(); i<i_end; ++i) {
      const FloatType &iobs = unmerged_data[i];
      const FloatType w = 1./(unmerged_sigmas[i]*unmerged_sigmas[i]);
      const cctbx::miller::index<> &hkl = unmerged_indices[i];
      add_observation(iobs, w, hkl);
    }
  }

  void merge(const std::string &sigma="population sigma") {
    indices.clear();
    data.clear();
    sigmas.clear();
    redundancies.clear();

    int sigma_flag = 0;
    if (sigma == "population sigma")
      sigma_flag = 1;
    else if (sigma == "experimental sigma")
      sigma_flag = 2;
    else
      throw std::runtime_error("Unknown sigma calculation method");

    for(typename boost::unordered_map<cctbx::miller::index<>, merging_data<FloatType> >::const_iterator p = data_.begin(); p != data_.end(); ++p) {
      indices.push_back(p->first);
      data.push_back(p->second.iobs);
      const int red = p->second.red;
      FloatType sigma = 0;

      if (sigma_flag == 1)
	sigma = std::sqrt(p->second.m2/p->second.weight) / std::sqrt(static_cast<FloatType>(red));
      else if (sigma_flag == 2)
	sigma = std::sqrt(1.0 / p->second.weight);

      sigmas.push_back(sigma);
      redundancies.push_back(red);
    }
  }

private:
  void add_observation(FloatType iobs, FloatType w, const cctbx::miller::index<> &hkl) {
    typename boost::unordered_map<cctbx::miller::index<>, merging_data<FloatType> >::iterator itr = data_.find(hkl);
    if (itr == data_.end()) {
      data_[hkl] = merging_data<FloatType>(iobs, w, w * iobs*iobs, 1);
    } else {
      const FloatType delta = iobs - itr->second.iobs;
      const FloatType tmp = w + itr->second.weight;
      const FloatType R = w * delta / tmp;
      itr->second.iobs += R;
      itr->second.m2 += itr->second.weight * delta * R;
      itr->second.weight = tmp;
      itr->second.red += 1;
    }
  }

protected:
  boost::unordered_map<cctbx::miller::index<>, merging_data<FloatType> > data_;
};

#endif // YAMTBX_DATAPROC_CRYSTFEL_MERGE_EQUIVALENTS_H
