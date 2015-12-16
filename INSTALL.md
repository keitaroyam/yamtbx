# Install yamtbx as a part of cctbx
## Dependencies
* [Python 2.7](https://www.python.org/)
* [CCTBX](http://cctbx.sourceforge.net/) with [CBFlib](http://www.bernstein-plus-sons.com/software/CBF/)
* [wxPython 2.8](http://www.wxpython.org/)
* [Numpy](http://www.numpy.org/)
* [Scipy](http://www.scipy.org/)
* [Matplotlib 1.3](http://matplotlib.org/)
* [Networkx](https://networkx.github.io/)

## 1. Prepare python environment
Set `$PREFIX` to any location you like (e.g. `/usr/local/yamtbx`)

### Python 2.x
```
wget http://www.python.org/ftp/python/2.7.5/Python-2.7.5.tar.bz2
tar xvf Python-2.7.5.tar.bz2
cd Python-2.7.5
./configure --prefix=$PREFIX/python-2.7.5
make -j8
make install
```

and prepare easy_install.

```
wget http://peak.telecommunity.com/dist/ez_setup.py
$PREFIX/python-2.7.5/bin/python ez_setup.py
```

### Numpy
You may find later version at numpy [website](http://sourceforge.net/projects/numpy/files/NumPy/)

```
wget "http://downloads.sourceforge.net/project/numpy/NumPy/1.7.1/numpy-1.7.1.tar.gz?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fnumpy%2Ffiles%2FNumPy%2F1.7.1%2F&ts=1380602990&use_mirror=jaist"
tar xvf numpy-1.7.1.tar.gz
cd numpy-1.7.1
$PREFIX/python-2.7.5/bin/python setup.py install
```

### Scipy

You may find later version at scipy [website](http://sourceforge.net/projects/scipy/files/scipy/).
Need to install blas-devel and lapack-devel packages in advance.

```
wget "http://downloads.sourceforge.net/project/scipy/scipy/0.12.0/scipy-0.12.0.tar.gz?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fscipy%2Ffiles%2Fscipy%2F0.12.0%2F&ts=1380603179&use_mirror=jaist"
tar xvf scipy-0.12.0.tar.gz
cd scipy-0.12.0
$PREFIX/python-2.7.5/bin/python setup.py install
```

### Matplotlib

```
$PREFIX/python-2.7.5/bin/easy_install matplotlib
```

### wxPython

```
wget "http://downloads.sourceforge.net/project/wxpython/wxPython/2.8.12.1/wxPython-src-2.8.12.1.tar.bz2?r=http%3A%2F%2Fwww.wxpython.org%2Fdownload.php&ts=1380604137&use_mirror=jaist"
tar xvf wxPython-src-2.8.12.1.tar.bz
cd wxPython-src-2.8.12.1
mkdir -p $PREFIX/wxPython-2.8.12.1
./configure --prefix=$PREFIX/wxPython-2.8.12.1 --enable-shared \
                      --with-gtk \
                      --with-gnomeprint \
                      --with-opengl \
                      --enable-debug \
                      --enable-debug_gdb \
                      --enable-geometry \
                      --enable-graphics_ctx \
                      --enable-sound --with-sdl \
                      --enable-mediactrl \
                      --enable-display \
make -j2; make install
cd contrib; make -j2; make install
cd ../wxPython
$PREFIX/python-2.7.5/bin/python setup.py install WX_CONFIG=$PREFIX/wxPython-2.8.12.1/bin/wx-config
```

## 2. Prepare cctbx and yamtbx

```
mkdir -p $PREFIX/cctbx/bundle
cd $PREFIX/cctbx/bundle
wget http://cci.lbl.gov/cctbx_build/results/last_published/cctbx_bundle.tar.gz
tar xvf cctbx_bundle.tar.gz
cd ..
svn checkout http://svn.code.sf.net/p/cctbx/code/trunk cctbx_svn
git clone https://github.com/keitaroyam/yamtbx.git yamtbx-github
mkdir -p build/src
cd build/src
ln -s ../../cctbx_svn/* .
ln -s ../../bundle/cctbx_sources/{annlib,annlib_adaptbx,boost,cbflib,ccp4io,ccp4io_adaptbx,scons} .
ln -s ../../yamtbx-github/yamtbx .
cd ..
$PREFIX/python-2.7.5/bin/python src/libtbx/configure.py --enable-openmp-if-possible=True cctbx mmtbx iotbx spotfinder dxtbx cbflib yamtbx
source setpaths.sh
make
```

and you may find yamtbx.\* or kamo\* commands in $PREFIX/cctbx/build/bin/.
