######################################################################
Python Computer Graphics Kit v1.2.0
Copyright (C) 2002, Matthias Baas (see license.txt)
######################################################################

The Python Computer Graphics Kit is a collection of Python modules
that contain the basic types and functions to be able to create 3D
computer graphics images. The kit mainly focuses on Pixar's RenderMan
interface, but some modules can also be used for OpenGL programs or
non-RenderMan compliant renderers like POV-Ray, for example.

Installing the package
----------------------

The package uses the Python distutils, so compiling and installing
looks the same on every platform, you simply have to call:

python setup.py install

This will compile the C-modules and install everything in the standard
location. See the distutils documentation that are coming with Python
if you have to install somewhere else.
If you've updated the sources for cgkit I recommend to delete the "build"
directory before compiling so that all the old code is out of the way.

The cgtypes module uses Pyrex (v0.9 or higher required) to generate
the C source file. The generated C file is included in the source
package, so you only need Pyrex if you want to modify the cgtypes module.

Once installed you can check the examples in the examples directory.

IMPORTANT: 

If you've installed a version before v1.0beta3 you have to remove
the directory Lib/site-packages/cgtypes in your Python distribution.
Otherwise Python will access the old Python types instead of the
new C types (which are located directly in Lib/site-packages).

Documentation
-------------

The documentation has to be downloaded separatly or you can browse
it online. Whatever you prefer, you should point your browser to:

http://cgkit.sourceforge.net

