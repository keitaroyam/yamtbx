// (c) RIKEN 2015. All rights reserved.
// Author: Keitaro Yamashita
//
// This software is released under the new BSD License; see LICENSE.
#include <boost/python.hpp>
#include <boost/python/args.hpp>
#include <string>
#include <pyublas/numpy.hpp>
#include <H5Cpp.h>
#include <cbf_hdf5_filter.h>

namespace bp = boost::python;

herr_t write_byteoffset_data(hid_t group_id, const std::string &name, const pyublas::numpy_vector<int> &data, int width, int height)
{
  hsize_t dim[2] = {height, width};
  H5::DataSpace space2D(2, dim);

  unsigned int cd_values[CBF_H5Z_FILTER_CBF_NELMTS];
  hsize_t chunk2[2] = {dim[0], dim[1]};
  cd_values[CBF_H5Z_FILTER_CBF_COMPRESSION] = CBF_BYTE_OFFSET; // byte_offset
  cd_values[CBF_H5Z_FILTER_CBF_RESERVED]    = 0;
  cd_values[CBF_H5Z_FILTER_CBF_BINARY_ID]   = 1; // ?
  cd_values[CBF_H5Z_FILTER_CBF_PADDING]     = 4095; // ?
  cd_values[CBF_H5Z_FILTER_CBF_ELSIZE]      = 4; // always 32 bit?
  cd_values[CBF_H5Z_FILTER_CBF_ELSIGN]      = 1; // 1 for signed, 0 if unsigned
  cd_values[CBF_H5Z_FILTER_CBF_REAL]        = 0; // 1 if a real array, 0 if an integer array
  cd_values[CBF_H5Z_FILTER_CBF_DIMFAST]     = dim[1];
  cd_values[CBF_H5Z_FILTER_CBF_DIMMID]      = dim[0];
  cd_values[CBF_H5Z_FILTER_CBF_DIMSLOW]     = 1;
  hid_t valprop = H5Pcreate(H5P_DATASET_CREATE);
  herr_t status = H5Pset_chunk(valprop, 2, chunk2);
  status = H5Pset_filter(valprop, CBF_H5Z_FILTER_CBF, H5Z_FLAG_OPTIONAL,
			 CBF_H5Z_FILTER_CBF_NELMTS, cd_values);
  // Don't forget to set HDF5_PLUGIN_PATH!
  if (status < 0) {
    std::cout << "Error in H5Pset_filter\n==============================\n";
    H5Eprint(H5E_DEFAULT, stdout);
    std::cout << "==============================\n";
  }

  hid_t dset = H5Dcreate(group_id, name.c_str(), H5T_STD_I32LE, // always 32 bit?
			 space2D.getId(), H5P_DEFAULT, valprop, H5P_DEFAULT);
  status = H5Dwrite(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, // always 32 bit?
		    &data[0]); // is this really a safe way??

  return status;
}

BOOST_PYTHON_MODULE(yamtbx_byteoffset_h5_ext)
{
  bp::def("write_byteoffset_data", write_byteoffset_data);
}
