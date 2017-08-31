# Install yamtbx as a part of cctbx

If you just want to use KAMO, please jump to [KAMO documentation](doc/kamo-en.md#installation).

## Dependencies
* [Python 2.7](https://www.python.org/)
* [CCTBX](https://cctbx.github.io/) with [CBFlib](http://www.bernstein-plus-sons.com/software/CBF/)
* [DIALS](https://dials.github.io)
* [wxPython 2.8](http://www.wxpython.org/)
* [Numpy](http://www.numpy.org/)
* [Scipy](http://www.scipy.org/)
* [Matplotlib 1.3](http://matplotlib.org/)
* [Networkx](https://networkx.github.io/)

## Installation using bootstrap.py

CCTBX-provided bootstrap.py can prepare CCTBX+DIALS environment with the dependencies (For details see [Installation for Developers - DIALS](https://dials.github.io/documentation/installation_developer.html)).
Here we use this script and install yamtbx into the built environment.

Set `$PREFIX` to any location you like (e.g. `/usr/local/yamtbx`)

```
mkdir -p $PREFIX
cd $PREFIX
wget https://raw.githubusercontent.com/cctbx/cctbx_project/master/libtbx/auto_build/bootstrap.py
python bootstrap.py --builder=dials
git clone https://github.com/keitaroyam/yamtbx.git yamtbx-github
cd modules/
ln -s ../yamtbx-github/yamtbx .
cd ..
./build/bin/libtbx.configure --enable-openmp-if-possible=True yamtbx
#source setpaths.sh
make # To build C++ codes in yamtbx, you may need to fix SConscript etc.. instruction will be given later (sorry).
```

and you will find yamtbx.\* or kamo\* commands in $PREFIX/build/bin/.
This brings other commands including phenix.\* and dials.\* and you may not want these (you may rather want to use those commands in PHENIX or DIALS environment). Just remove unwanted command scripts from $PREFIX/build/bin/. The alternative way is to make bin2 directory and make links to yamtbx.\* and kamo\* there, and set PATH to bin2/.