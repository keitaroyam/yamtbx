# EIGER data format and processing

This document explains about EIGER X 9M at BL32XU.
The contents may be updated by the upgrade of the firmware of EIGER.

## EIGER HDF5

Unlike conventional image formats like .img or .cbf, single HDF5 (.h5) file can contain more than one image.
HDF5 has a hierarchical structure like a file system, each of which has arbitrary data including image data and experimental parameters.

HDF5 files written by EIGER consist of two types; master and data. For example, if you collected diffraction data with prefix of sample, you will get

* sample_master.h5 - experimental and instrumental parameters
* sample_data_??????.h5 - diffraction images

Data h5 file contains multiple (by default, up to 100 images) diffraction images.
In case you collected less than 100 images, you will only get sample_data_000001.h5. In case of 200 images or less, you will additionally have sample_data_000002.h5.

sample_master.h5 contains the link to data h5 (using hdf5's ExternalLink function), which allows you to transparently read diffraction images by opening master h5.
Actual data exist in data h5, but you don't have to be aware of this. This is why the data can be processed by just giving master h5 to the program.

### Compression of internal data (filter)

HDF5 has a mechanism called 'filter' that internally compresses data. EIGER data are compressed with LZ4 or by default bitshuffle+LZ4 (bslz4; 30-50% smaller compared to lz4) algorithm.
Normally master.h5, which include a bit huge data e.g. flatfield and pixel mask, is not compressed, but at BL32XU byteshuffle+gzip is used to reduce the size (bslz4 is used since 2017A).

No matter what kind of filters is used, you can read the data in the same way (HDF5 library internally decompress data). However, for non-standard filters including bslz4, you need to install the plugin.


## Installing software

First, you need these:

### HDF5 software
generate_XDS.INP uses h5dump in [HDF5 software](https://www.hdfgroup.org/HDF5/release/obtain5.html).
You can use several package managers (like macports for Mac) to install hdf5.

[Hdfview](https://www.hdfgroup.org/products/java/release/download.html) is also a useful software to visualize the contents of master hdf5 file.

### eiger2cbf (H5ToXds compatible)
eiger2cbf is a tool to convert h5 files to CBF format, which is a standard format for PILATUS.
You can find the link (Converter for Eiger files) in the bottom of [Downloads - Mosflm](http://www.mrc-lmb.cam.ac.uk/harry/imosflm/ver721/downloads.html#Eiger2CBF) to get binaries for Linux and Mac.

eiger2cbf can be used as an alternative of H5ToXds, which is needed by XDS to process h5 files.
H5ToXds itself can be obtained from the bottom of [Dectris website](https://www.dectris.com/EIGER_X_Features.html), but only Linux version is available.
Moreover, H5ToXds does not write header in CBF file.
By these reasons, I recommend to use eiger2cbf instead of H5ToXds.
Eiger2cbf shows message in standard output, which overlaps with XDS output in processing. To avoid this, create shell script named H5ToXds like below.
```
#!/bin/sh
eiger2cbf $@ 2>/dev/null
```

To use eiger2cbf, just give a master.h5; e.g.,

* `eiger2cbf sample_master.h5` .. shows the number of images in the file
* `eiger2cbf sample_master.h5 10 sample_10.cbf` .. Save 10th image (fist image is 1) as sample_10.cbf
* `eiger2cbf sample_master.h5 1:100 sample` .. Save 1st-100th images as sample_??????.cbf

### bitshuffle plugin (advanced)
Although not necessarily indispensable, if you installed [bitshuffle](https://github.com/kiyo-masui/bitshuffle) plugin you can open bslz4-compressed hdf5 files through h5dump, hdfview, and adxv.
Development environment would be needed for installation (Xcode? On Mac [unverified])

```
(sudo) easy_install cython
(sudo) easy_install h5py
cd ~/tmp (arbitrary)
git clone https://github.com/kiyo-masui/bitshuffle.git
cd bitshuffle/bitshuffle
cython ext.pyx
cython h5.pyx
cd ..
(sudo) python setup.py install --h5plugin
```
By default it is installed to /usr/local/hdf5/lib/plugin. To change it, add `--h5plugin-dir=` after `--h5plugin`.
You will need to set environmental variable `HDF5_PLUGIN_PATH` to where you installed (In bash, use `export HDF5_PLUGIN_PATH=`; in csh use `setenv HDF5_PLUGIN_PATH `).

### Neggia plugin (only if you want to use plugin function in XDS)
With plugin, XDS can process h5 files without H5ToXds, which means temporary cbf files are not written. Neggia, the official plugin by DECTRIS, can be obtained from [DECTRIS website](https://www.dectris.com/neggia.html) where you will need registration. Alternatively, you can obtain the source code from [dectris/neggia - Github](https://github.com/dectris/neggia) to build it yourself.


## Viewing image
You can always use any software including adxv by converting images to cbf format using eiger2cbf.

To directly read hdf5, these programs are avialble:

* [ALBULA](https://www.dectris.com/Albula_Overview.html)
* [Adxv](http://www.scripps.edu/tainer/arvai/adxv.html) (need [bitshuffle plugin](#bitshuffle-plugin-advanced); but resolution may be incorrect because experimental parameters are ignored)
* dials.image_viewer (bundled with DIALS)
* yamtbx.adxv_eiger (the default program at BL32XU. launcher for adxv)

yamtbx.adxv_eiger can be installed by following instruction; you can already use it if you installed [KAMO](kamo-en.md#how-can-i-use-kamo-at-home) or yamtbx.

1. Install [PHENIX](http://www.phenix-online.org/)-1.11 or newer
2. Run the following commands (you can clone yamtbx anywhere you like)
```
cd $HOME
git clone https://github.com/keitaroyam/yamtbx.git
cd $PHENIX/modules
ln -s ~/yamtbx/yamtbx .
cd ../build
libtbx.configure yamtbx
```


## Data processing

Two ways to process EIGER data:

1. give h5 files directly (XDS, DIALS, HKL-2000)
2. convert to cbf format first (iMosflm)

### XDS

You can use [generate_XDS.INP](http://strucbio.biologie.uni-konstanz.de/xdswiki/index.php/Generate_XDS.INP) (rev. 0.63 or newer) to prepare XDS.INP (NEED h5dump).
Give a master h5 file to generata\_XDS.INP.
```
generate_XDS.INP ../sample_master.h5
```

If you manually prepare XDS.INP, `NAME_TEMPLATE_OF_DATA_FRAMES=` should be a file name of master h5, but replace 'master' with '??????' like:
```
NAME_TEMPLATE_OF_DATA_FRAMES= ../sample_??????.h5
```

In processing, H5ToXds or plugin is needed.

#### Using H5ToXds
H5ToXds is a helper program for XDS, which generates temporary cbf file. As mentioned above, I recommend to use eiger2cbf as an alternative of H5ToXds.

#### Using plugin
With plugin, XDS can directly process h5 file. The speed-up (about 10%) can be expected especially for the I/O-limited environment.
In XDS.INP, specify plugin (.so file) location by
```
LIB= /where-you-installed/plugin-name.so
```
For example [Neggia plugin](#neggia-plugin-only-if-you-want-to-use-plugin-function-in-xds) can be used.

**CAUTION!** In case of Neggia plugin, you cannot process data collected in 2016 at BL32XU because Neggia does unofficial re-implementation of HDF5 which does not support GZIP+SHUF filters.
In this case, you need to convert master h5 file; if you have yamtbx,
```
mv yours_master.h5 yours_master.h5.org
yamtbx.eiger_reduce_master yours_master.h5.org h5out=yours_master.h5 compress=bslz4
```
If you don't have yamtbx, give [this script](https://github.com/keitaroyam/yamtbx/blob/master/yamtbx/dataproc/command_line/eiger_reduce_master.py) to `phenix.python` (ver. 1.11 or newer).

### DIALS
Ver. dev-652 (1.dev.193) or newer supports EIGER hdf5 file at BL32XU.
Please use the latest version, or at least, the version bundled with CCP4 7.0 update 013 or PHENIX-1.11.
The latest version can be obtained from [DIALS website](https://dials.github.io/installation.html).

Please see the [tutorial](https://dials.github.io/documentation/tutorials/index.html); but first give master h5 file to dials.import:
```
dials.import ../sample_master.h5
```

**CAUTION!** DIALS (before ver. 1.5) had gotten I(+) and I(-) the wrong way round, which resulted in flipped signs of anomalous differences. This bug was fixed in the ver. 1.5 released in April 2017. It doesn't matter if you don't use Bijvoet differences; but you cannot get correct result if you want to see anomalous difference Fourier maps or phase with SAD/MAD/SIRAS/MIRAS methods. After dials.import, you need to edit datablock.json. Find
```json
        "panels": [
          {
            "origin": [
              119.39999904559726,
              119.92500419409488,
              -119.99999428590824
            ],
            "fast_axis": [
              -1.0,
              0.0,
              0.0
            ],
            "name": "/entry/instrument/detector",
            "raw_image_offset": [
              0,
              0
            ],
            "slow_axis": [
              -0.0,
              -1.0,
              0.0
            ],
```
and change `fast_axis` from `-1,0,0` to `1,0,0`, and flip the sign of the first value of `origin` (`origin` depends on data as they express beam center and camera distance). For the example above, it should be fixed like this:
```json
        "panels": [
          {
            "origin": [
              -119.39999904559726,
              119.92500419409488,
              -119.99999428590824
            ],
            "fast_axis": [
              1.0,
              0.0,
              0.0
            ],
            "name": "/entry/instrument/detector",
            "raw_image_offset": [
              0,
              0
            ],
            "slow_axis": [
              -0.0,
              -1.0,
              0.0
            ],
```


### iMosflm

Use [eiger2cbf](#eiger2cbf-h5toxds-compatible) to convert to cbf format.

### HKL-2000

Since ver. 714 (Sep 2016), hdf5 can be directly processed. The previous versions need converted cbf files.

## References
* [EIGER X series (Dectris official)](https://www.dectris.com/EIGER_X_Detectors.html#main_head_navigation)
* [Eiger - XDSwiki](http://strucbio.biologie.uni-konstanz.de/xdswiki/index.php/Eiger)
* [HDF5 - The HDF Group](https://www.hdfgroup.org/HDF5/)
* [NeXus](http://www.nexusformat.org/) In future EIGER hdf5 will be NeXus comaptible (currently not perfect)
