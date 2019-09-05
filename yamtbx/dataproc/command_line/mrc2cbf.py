import mrcfile
import os
from yamtbx.dataproc import cbf

def run(mrcin, prefix=None):
    f = mrcfile.open(mrcin)
    if prefix is None: prefix = os.path.splitext(os.path.basename(mrcin))[0]
    data = f.data
    if data.ndim == 2: data = data.reshape(-1, *data.shape)
    
    print "MRC file loaded"
    print " shape = %s" % (data.shape,)
    print " dtype = %s" % data.dtype
    print
    
    for i in xrange(data.shape[0]):
        size2, size1 = data[i].shape # XXX really?
        cbfout = "%s_%.3d.cbf"%(prefix, i+1)
        d = data[i]
        if d.dtype.kind == "f":
            print "Scaling data x1000"
            d = (data[i]*1000.).astype("int32")
        print "Writing %s" % cbfout
        cbf.save_numpy_data_as_cbf(d.flatten(), size1, size2,
                                   "%s:%d"%(mrcin, i+1),
                                   cbfout, pilatus_header="""
# Detector: %(detname)s
# Pixel_size %(pixelsize)e m x %(pixelsize)e m
# Wavelength %(wavelength)f A
# Detector_distance %(distance)f m
# Beam_xy (%(beamx)f, %(beamy)f) pixels
# Start_angle %(start_angle)%.4f deg.
# Angle_increment %(angle_inc).4f deg.
""" % dict(detname="????",
           pixelsize=75.e-6,
           wavelength=1,
           distance=300.e-3,
           beamx=size1/2., beamy=size2/2.,
           start_angle=0, angle_inc=0))
    

if __name__ == "__main__":
    import sys
    mrcin = sys.argv[1]
    
    run(mrcin)
