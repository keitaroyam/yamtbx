// (c) RIKEN 2015. All rights reserved. 
// Author: Keitaro Yamashita
//
// This software is released under the new BSD License; see LICENSE.

#include <boost/python.hpp>
#include <boost/python/args.hpp>
#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <scitbx/array_family/flex_types.h>
#include <cctbx/miller.h>
#include <cctbx/miller/match.h>
#include <cctbx/miller/index_span.h>
#include <algorithm>
#include <functional>
#include <map>
#include <vector>

namespace bp = boost::python;
namespace af = scitbx::af;


typedef cctbx::miller::index <> miller_index_type;
typedef af::const_ref <miller_index_type> flex_miller_index;

struct miller_packed_less_than
{
  miller_packed_less_than(const flex_miller_index &indices)
    : span(indices) {
  }
  bool operator()(const miller_index_type &lhs, const miller_index_type &rhs) {
    return span.pack(lhs) < span.pack(rhs);
  }
private:
  cctbx::miller::index_span span;
};

/* Very specific function for add_saturation_info_to_mtz.py
 * "data" are assumed to be ZD or ZOBS (spot centroid on frame number unit)
 * This tries to find the same reflection (indices must be the same, frame number must be similar).
 * PRIOR CONDITION: rhs must be sorted!!
 */
af::flex_bool
make_selection_for_xds_unmerged(const flex_miller_index &lhs_indices,
                                const af::const_ref<double> &lhs_data,
                                const flex_miller_index &rhs_indices,
                                const af::const_ref<double> &rhs_data,
                                double threshold)
{
  CCTBX_ASSERT(lhs_indices.size() == lhs_data.size());
  CCTBX_ASSERT(rhs_indices.size() == rhs_data.size());

  af::flex_bool ret(rhs_data.size());
  for (int i = 0; i < ret.size(); ++i)
    ret[i] = false;

  for (int i = 0; i < lhs_indices.size(); ++i) {
    //    std::cout << i << "/" << lhs_indices.size() << "\n";

    const miller_index_type *p = std::lower_bound(rhs_indices.begin(), rhs_indices.end(), lhs_indices[i],
                                                  miller_packed_less_than(rhs_indices));

    //    std::cout << (*p)[0] << " " << (*p)[1] << " " << (*p)[2]
    //              << " ?= " << lhs_indices[i][0] << " " << lhs_indices[i][1] << " " << lhs_indices[i][2] << "\n";
    int found = 0;
    for (; p != rhs_indices.end(); ++p) {
      if (*p != lhs_indices[i])
        break;

      const int j = p - rhs_indices.begin();

      if (std::abs(rhs_data[j] - lhs_data[i]) < threshold) {
        ret[j] = true;
        ++found;
      }
      // else
      //   std::cout << "Diff= " << std::abs(rhs_data[j] - lhs_data[i]) << " "
      //             << "(" << (*p)[0] << " " << (*p)[1] << " " << (*p)[2] << ") "
      //             << rhs_data[j] << ":" << lhs_data[i] << " j=" << j
      //             << "\n";
    }
    if (found != 1)
      std::cout << found << " Found\n";

    /* // This loop is slow!
    for (int j = 0; j < rhs_indices.size(); ++j) {
      if (lhs_indices[i] != rhs_indices[j])
        continue;
      if (std::abs(rhs_data[i] - lhs_data[j]) < threshold)
        ret[j] = true;
    }
    */
  }
  return ret;
}

struct comp {
  comp(af::shared<cctbx::miller::index<> > const& ind): ind(ind) {}
  bool operator()(size_t i, size_t j) { return _comp(ind[i],ind[j]); }
  cctbx::miller::fast_less_than<> _comp;
  af::shared<cctbx::miller::index<> > const& ind;
};

af::shared<size_t>
sort_permutation_fast_less(af::shared<cctbx::miller::index<> > const& indices)
{
  af::shared<size_t> ret(indices.size());
  for (int i = 0; i < indices.size(); ++i) ret[i] = i;

  std::sort(ret.begin(), ret.end(), comp(indices));
  return ret;
}

boost::python::tuple
my_common_indices(af::shared<cctbx::miller::index<> > const& indices_0,
		  af::shared<cctbx::miller::index<> > const& indices_1)
{
  cctbx::miller::fast_less_than<> lt;

  af::shared<size_t> ret0, ret1;
  ret0.reserve(indices_0.size());
  ret1.reserve(indices_1.size());


  for (int i = 0, j = 0; i < indices_0.size() && j < indices_1.size(); ++j) {
    while (i < indices_0.size() && lt(indices_0[i], indices_1[j])) ++i;
    if (i < indices_0.size() && indices_0[i] == indices_1[j]) {
      ret0.push_back(i);
      ret1.push_back(j);
    }
  }
  return boost::python::make_tuple(ret0, ret1);
}

BOOST_PYTHON_MODULE(yamtbx_utils_ext)
{

  typedef bp::return_value_policy<bp::return_by_value> rbv;

  bp::def("make_selection_for_xds_unmerged", make_selection_for_xds_unmerged);
  bp::def("sort_permutation_fast_less", sort_permutation_fast_less, rbv());
  bp::def("my_common_indices", my_common_indices, rbv());
}
