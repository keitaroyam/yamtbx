#include <cctbx/boost_python/flex_fwd.h>

//#include <cctbx/miller/sym_equiv.h>
//#include <cctbx/miller/union_of_indices.h>
#include <cctbx/miller/math.h>
#include <boost/python.hpp>
#include <scitbx/boost_python/container_conversions.h>

namespace yamtbx { namespace dataproc { namespace crystfel { namespace boost_python {
  void wrap_merge_equivalents();

namespace {
  void init_module()
  {
    using namespace boost::python;

    wrap_merge_equivalents();
  }
} // namespace <anonymous>
}}}} // namespace yamtbx::dataproc::crystfel::boost_python


BOOST_PYTHON_MODULE(yamtbx_dataproc_crystfel_ext)
{
  yamtbx::dataproc::crystfel::boost_python::init_module();
}
