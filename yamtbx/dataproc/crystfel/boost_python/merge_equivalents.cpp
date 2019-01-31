#include <cctbx/boost_python/flex_fwd.h>
#include <yamtbx/dataproc/crystfel/merge_equivalents.h>

#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

namespace yamtbx { namespace dataproc { namespace crystfel { namespace boost_python {

namespace {

  struct merge_equivalents_crystfel_wrappers
  {
    typedef merge_equivalents_crystfel<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("merge_equivalents_crystfel")
	.def("add_observations", (void(w_t::*)(af::const_ref<cctbx::miller::index<> > const&,
					       af::const_ref<double> const&))&w_t::add_observations,
	     (arg("unmerged_indices"), arg("unmerged_data")))
	.def("add_observations", (void(w_t::*)(af::const_ref<cctbx::miller::index<> > const&,
					       af::const_ref<double> const&,
					       af::const_ref<double> const&))&w_t::add_observations,
	     (arg("unmerged_indices"), arg("unmerged_data"), arg("unmerged_sigmas")))
	.def("merge", &w_t::merge, arg("sigma")="population sigma")
        .add_property("indices", make_getter(&w_t::indices, rbv()))
        .add_property("data", make_getter(&w_t::data, rbv()))
        .add_property("sigmas", make_getter(&w_t::sigmas, rbv()))
        .add_property("redundancies", make_getter(&w_t::redundancies, rbv()))
      ;
    }
  };
} // namespace <anoymous>

  void wrap_merge_equivalents()
  {
    merge_equivalents_crystfel_wrappers::wrap();
  }

}}}} // namespace yamtbx::dataproc::crystfel::boost_python
