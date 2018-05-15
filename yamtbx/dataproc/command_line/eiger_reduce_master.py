"""
(c) RIKEN 2017. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

import iotbx.phil
import os
import subprocess
import tempfile
import shutil
import h5py
import bitshuffle.h5

master_params_str = """ 
h5in = None
 .type = path
 .help = "Input file"
h5out = None
 .type = path
 .help = "Output file"
replace = False
 .type = bool
 .help = "Replace original file with updated file"
remove_filters_first = False
 .type = bool
 .help = "Run h5repack -f NONE to remove all filters first."
compress = bslz4 *bslz4_and_gzipshuf
 .type = choice(multi=False)
 .help = "Compress large datasets in master.h5. bslz4_and_gzipshuf is to apply bslz4 for pixel_mask and gzip+shuf for other datasets."
remove_detectorModule_data = flatfield pixel_mask trimbit
 .type = choice(multi=True)
 .help = "Optionally remove data in /entry/instrument/detector/detectorSpecific/detectorModule_* flatfield and pixel_mask are redundant information."
remove_overall = flatfield pixel_mask
 .type = choice(multi=True)
 .help = "Optionally remove flatfield and pixel_mask data in /entry/instrument/detector/detectorSpecific/. NOTE they (especially pixel_mask) may be necessary data."
"""

class ReduceMaster:
    def __init__(self, h5in):
        self.h = h5py.File(h5in, "a")
    # __init__()

    def remove_redundant(self, kinds):
        assert set(kinds).issubset(("flatfield","pixel_mask","trimbit"))

        detSp = self.h["/entry/instrument/detector/detectorSpecific"]
        for k in detSp.keys():
            if k.startswith("detectorModule_"):
                for kk in kinds:
                    if kk not in detSp[k]: continue
                    print "Removing", k+"/"+kk
                    del detSp[k][kk]
    # remove_redundant()

    def remove_datasets(self, kinds):
        assert set(kinds).issubset(("flatfield","pixel_mask"))
        detSp = self.h["/entry/instrument/detector/detectorSpecific"]

        for k in kinds:
            print "Removing %s" % k
            del detSp[k]

    def find_large_dataset_visitor(self, path, obj):
        if type(obj) is h5py.Dataset and obj.size > 100:
            self.large_datasets.append(path)
    # find_large_dataset_visitor()

    def compress_large_datasets(self, compress):
        if not compress: return
        assert compress in ("bslz4", "bslz4_and_gzipshuf")

        self.large_datasets = []
        self.h.visititems(self.find_large_dataset_visitor)

        for path in self.large_datasets:
            print "Compressing %s (size=%d)" % (path, self.h[path].size)
            data = self.h[path][:]
            del self.h[path]
            
            if compress=="bslz4":
                self.h.create_dataset(path, data.shape,
                                      compression=bitshuffle.h5.H5FILTER,
                                      compression_opts=(0, bitshuffle.h5.H5_COMPRESS_LZ4),
                                      chunks=None, dtype=data.dtype, data=data)
            elif compress == "bslz4_and_gzipshuf":
                if "pixel_mask" in path: # bslz4
                    self.h.create_dataset(path, data.shape,
                                          compression=bitshuffle.h5.H5FILTER,
                                          compression_opts=(0, bitshuffle.h5.H5_COMPRESS_LZ4),
                                          chunks=None, dtype=data.dtype, data=data)
                else: # gzip+shuf
                    self.h.create_dataset(path, data.shape,
                                          compression="gzip",shuffle=True,
                                          chunks=None, dtype=data.dtype, data=data)
            else:
                raise "Never reaches here."

    # compress_large_datasets()

    def close(self): self.h.close()
# class ReduceMaster

def run(params):
    print "Parameters:"
    iotbx.phil.parse(master_params_str).format(params).show(prefix="  ")
    print

    master_size_org = os.path.getsize(params.h5in)

    if params.replace:
        shutil.copyfile(params.h5in, params.h5in+".org")
        params.h5out = params.h5in

    tmpfd, tmpout = tempfile.mkstemp(prefix=os.path.basename(params.h5out)+".tmp", dir=os.path.dirname(params.h5out))
    os.close(tmpfd)

    if params.remove_filters_first:
        try:
            p = subprocess.Popen(["h5repack","-f","NONE",params.h5in,tmpout], shell=False)
            p.wait()
        except OSError:
            print "h5repack failed. Is h5repack installed?"
            return
    else:
        shutil.copyfile(params.h5in, tmpout)

    try:
        redmas = ReduceMaster(tmpout)

        if params.remove_detectorModule_data:
            redmas.remove_redundant(params.remove_detectorModule_data)

        if params.remove_overall:
            redmas.remove_datasets(params.remove_overall)

        if params.compress:
            redmas.compress_large_datasets(params.compress)

        redmas.close()
        
        # Run h5repack to clean up the removed space
        try:
            p = subprocess.Popen(["h5repack",tmpout,params.h5out], shell=False)
            p.wait()
        except OSError:
            print "h5repack failed. Is h5repack installed?"
    finally:
        os.remove(tmpout)

    master_size_after = os.path.getsize(params.h5out)

    print
    print "Finished."
    print "  Original file: %s (%.2f MB)" % (params.h5in, master_size_org/1024**2)
    print " Generated file: %s (%.2f MB)" % (params.h5out, master_size_after/1024**2)
# run()

def print_help(command_name):
    print """\
This script (re)modifies master.h5 file of EIGER detectors written by DECTRIS software.
All features are optional:
- Remove filters
- Apply bslz4 filters to large datasets
- Remove unnecessary data to save the disk space
- Remove (maybe necessary but) large data to further save the disk space

You need h5repack program. You also need phenix-1.11 or later if you use phenix.python; dials.python may be used instead.

* Usage:
%(command_name)s yours_master.h5 [h5out=yours_master_reduced.h5] [remove_filters_first=True] [compress=bslz4_and_gzipshuf] [remove_detectorModule_data=flatfield+pixel_mask+trimbit]

* In case you're BL32XU user, collected data before 2017, and want to use Neggia plugin for XDS (and reduce file size anyway):
mv yours_master.h5 yours_master.h5.org
%(command_name)s yours_master.h5.org h5out=yours_master.h5 compress=bslz4 remove_detectorModule_data=flatfield+pixel_mask+trimbit

Default parameters:""" % dict(command_name=command_name)

    iotbx.phil.parse(master_params_str).show(prefix="  ", attributes_level=1)
# print_help()

def run_from_args(argv, command_name="phenix.python this-script"):
    if not argv or "-h" in argv or "--help" in argv:
        print_help(command_name)
        return

    cmdline = iotbx.phil.process_command_line(args=argv,
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args

    for arg in args:
        if not os.path.isfile(arg):
            print "File not found: %s" % arg
            return

        if params.h5in is None and arg.endswith(".h5"):
            params.h5in = arg

    if params.h5in is None:
        print "Please give _master.h5 file."
        return

    if params.h5out is None:
        params.h5out = os.path.splitext(params.h5in)[0] + "_reduce.h5"

    run(params)
# run_from_args()

if __name__ == "__main__":
    import sys

    run_from_args(sys.argv[1:])
