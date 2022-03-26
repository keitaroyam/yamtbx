"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import absolute_import, division, print_function, generators
import os
import pycbf
import numpy
from cbflib_adaptbx import cbf_binary_adaptor, CBFWriteAdaptor

def load_cbf_as_numpy(filein, quiet=True):
    assert os.path.isfile(filein)
    if not quiet:
        print("reading", filein, "as cbf")
    h = pycbf.cbf_handle_struct()
    h.read_file(filein.encode("utf-8"), pycbf.MSG_DIGEST)
    ndimfast, ndimslow = h.get_image_size_fs(0)
    arr = numpy.fromstring(h.get_image_fs_as_string(0, 4, 1, ndimfast, ndimslow), dtype=numpy.int32)
    return arr, ndimfast, ndimslow
# load_cbf_as_numpy()

def load_minicbf_as_numpy(filein, quiet=True): # This can also read XDS special cbf
    assert os.path.isfile(filein)
    if not quiet:
        print("reading", filein, "as minicbf")
    h = pycbf.cbf_handle_struct()
    h.read_file(filein.encode("utf-8"), pycbf.MSG_DIGEST)
    h.require_category(b"array_data")
    h.find_column(b"data")
    compression, binary_id, elsize, elsigned, elunsigned, elements, minelement, maxelement, bo, ndimfast, ndimmid, ndimslow, padding = h.get_integerarrayparameters_wdims()
    assert elsize == 4 or elsize == 8
    assert elsigned == 1
    assert ndimslow <= 1
    arr = numpy.fromstring(h.get_integerarray_as_string(), dtype=numpy.int32 if elsize==4 else numpy.int64)
    return arr, ndimfast, ndimmid

# load_minicbf_as_numpy()

def load_cbf_as_flex(filein): # This can also read XDS special cbf
    M = cbf_binary_adaptor(filein)
    data = M.uncompress_implementation("buffer_based").uncompress_data()
    nslow, nfast = M.dim_slow(), M.dim_fast() # can be obtained after getting data
    return data, nfast, nslow
# load_cbf_as_flex()

def load_xds_special(cbfin):
    h = pycbf.cbf_handle_struct()
    h.read_file(cbfin.encode("utf-8"), pycbf.MSG_DIGEST)
    h.require_category(b"array_data")
    h.find_column(b"header_contents")
    header = h.get_value()

    M = cbf_binary_adaptor(cbfin)
    data = M.uncompress_implementation(b"buffer_based").uncompress_data()
    #print "slow, fast=", M.dim_slow(), M.dim_fast() # can be obtained after getting data
    return header, data, M.dim_slow(), M.dim_fast()
# load_xds_special()

def save_numpy_data_as_cbf(data, size1, size2, title, cbfout, pilatus_header=None, header_convention="PILATUS_1.2"):
    h = pycbf.cbf_handle_struct()
    h.new_datablock(title.encode("utf-8"))

    h.require_category(b'array_data')

    if pilatus_header is not None:
        h.require_column(b'header_convention')
        h.set_value(b'"%s"'%header_convention.encode("utf-8"))
        h.require_column(b'header_contents')
        h.set_value(pilatus_header.encode("utf-8"))


    h.require_category(b'array_data')
    h.require_column(b'data')

    assert data.dtype.kind in "iu"
    elsigned = 1 if data.dtype.kind == "i" else 0

    h.set_integerarray_wdims_fs(pycbf.CBF_BYTE_OFFSET, 1, data.tostring(), data.dtype.itemsize,
                                elsigned, len(data), b"little_endian",
                                size1, size2, 1, 0)
    h.write_file(cbfout.encode("utf-8"), pycbf.CBF,
                 pycbf.MIME_HEADERS|pycbf.MSG_DIGEST|pycbf.PAD_4K, pycbf.ENC_NONE)
# save_numpy_data_as_cbf()

def save_flex_int_as_cbf(data, cbfout):
    writer = CBFWriteAdaptor(cbfout)
    writer.write_data(data)
# save_flex_int_as_cbf()

def get_pilatus_header(cbfin):
    h = pycbf.cbf_handle_struct()
    if cbfin.endswith(".bz2"):
        # TODO to speed up, better only bunzip2 the first part of file..
        import tempfile
        import bz2
        junk, tmpf = tempfile.mkstemp()
        os.close(junk)
        open(tmpf, "wb").write(bz2.BZ2File(cbfin).read())
        h.read_file(tmpf.encode("utf-8"), pycbf.MSG_DIGEST)
        os.remove(tmpf)
    else:
        h.read_file(cbfin.encode("utf-8"), pycbf.MSG_DIGEST)
    h.require_category(b"array_data")
    h.find_column(b"header_contents")
    header = h.get_value()
    return header
# get_pilatus_header()
