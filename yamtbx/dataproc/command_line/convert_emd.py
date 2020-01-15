import os
import datetime
import h5py
import numpy
import json
import iotbx.phil
from yamtbx.dataproc import cbf
from yamtbx.util.xtal import electron_voltage_to_wavelength

master_params_str = """\
clen = .85
 .type = float
 .help = camera distance in mm
beamxy = None
 .type = floats(size=2)
 .help = Beam center in pixel
pixel_size = 14
 .type = float
 .help = "Unbinned pixel size in um"
detector = "BM-Ceta"
 .type = str
 .help = "Detector name in metadata"

tilted_only = False
 .type = bool
 .help = Only convert frames with constant speed
frame_range = None
 .type = ints(size=2)
 .help = Override range 
offset = 0
 .type = int
 .help = offset (pedestal) value to pixel values
"""


def get_metadata(metadata):
    mds = []
    for i in xrange(metadata.shape[1]):
        metadata_array = metadata[:, i].T
        mdata_string = metadata_array.tostring().decode("utf-8")
        mds.append(json.loads(mdata_string.rstrip('\x00')))

    return mds
# get_metadata()

def dump_and_analyse_time_angle(metadata, prefix, detector):
    times, alphas = [], []

    ofs = open(prefix+"_time_alpha.dat", "w")
    ofs.write("frame date alpha\n")
    
    for i, md in enumerate(metadata):
        t = datetime.datetime.fromtimestamp(int(md["CustomProperties"]["Detectors[%s].TimeStamp"%detector]["value"])/1.e6)
        alpha = numpy.rad2deg(float(md["Stage"]["AlphaTilt"]))
        ofs.write('%6d "%s" %s\n' % (i+1, t, alpha))
        times.append(t)
        alphas.append(alpha)

        if i > 0 :
            print "%6d %+f %+f" % (i+1, (times[-1]-times[-2]).total_seconds(), alphas[-1]-alphas[-2])

    ofs.close()

    d_times = numpy.array([(times[i]-times[i-1]).total_seconds() for i in range(1, len(times))])
    d_alphas = numpy.diff(alphas)
    print d_times
    print d_alphas
    if len(d_alphas) == 0: return [], float("nan"), float("nan")
    q25, q50, q75 = numpy.percentile(d_alphas, [25, 50, 75])
    iqr = q75-q25
    iqrc = 1.5
    lowlim, highlim = q25 - iqrc*iqr, q75 + iqrc*iqr
    #print "outliers=", d_alphas[numpy.where(numpy.logical_or(d_alphas<lowlim, d_alphas>highlim))]
    d_alphas2 = d_alphas[numpy.where(numpy.logical_and(d_alphas>lowlim, d_alphas<highlim))] # outlier rejected
    d_alpha_z = abs(d_alphas-numpy.mean(d_alphas2))/numpy.std(d_alphas2)

    valid_range = [0, len(metadata)-1]
    for i in xrange(len(metadata)-1):
        if d_alpha_z[i] < 3: break
        valid_range[0] = i+1

    for i in reversed(xrange(len(metadata)-1)):
        if d_alpha_z[i] < 3: break
        valid_range[1] = i

    if valid_range[0] > valid_range[1]:
        valid_range = [0, len(metadata)-1] # reset
        
    print "valid_range=", valid_range

    print "time_diff= %.6f +/- %.9f seconds" % (numpy.mean(d_times), numpy.std(d_times))
    print "d_alphas= %.6f +/- %.9f degrees" % (numpy.mean(d_alphas2), numpy.std(d_alphas2))

    mean_alpha_step = (alphas[valid_range[1]] - alphas[valid_range[0]])/(valid_range[1]-valid_range[0]+1)
    mean_time_step = (times[valid_range[1]] - times[valid_range[0]]).total_seconds()/(valid_range[1]-valid_range[0]+1)

    return valid_range, mean_alpha_step, mean_time_step
# dump_time_angle()

def write_minicbf(data, metadata, prefix, frame_range, alpha_step, time_step, params):

    detector = str(metadata[0]["BinaryResult"]["Detector"])
    voltage = float(metadata[0]["Optics"]["AccelerationVoltage"])
    wavelen = electron_voltage_to_wavelength(voltage)
    #width = int(metadata[0]["BinaryResult"]["ImageSize"]["width"]) #/ int(metadata[0]["Detectors"]["ImageSize"]["width"])
    #height = int(metadata[0]["BinaryResult"]["ImageSize"]["height"]) # /  binning!!
    width, height = data.shape[:2] # the other way around?
    binning = int(metadata[0]["BinaryResult"]["ImageSize"]["width"])//width
    exp_time = float(metadata[0]["Detectors"]["Detector-0"]["ExposureTime"])
    px = params.pixel_size*1.e-6 * binning

    clen = params.clen
    if params.beamxy is None:
        beamx, beamy = (width/2., height/2.)
    else:
        beamx, beamy = params.beamxy

    print "There are %d frames" % data.shape[2]

    for i in xrange(data.shape[2]):
        if i < frame_range[0] or i > frame_range[1]:
            print "Ignoring frame %d" % (i+1)
            continue
        
        print "Converting", i+1
        frame_date = datetime.datetime.fromtimestamp(int(metadata[i]["CustomProperties"]["Detectors[%s].TimeStamp"%params.detector]["value"])/1.e6)
        datestr = datetime.datetime.strftime(frame_date, "%Y-%m%-dT%H:%M:%S.%f")
        start_angle = str(metadata[i]["Stage"]["AlphaTilt"])
        cbf.save_numpy_data_as_cbf(data[:,:,i].flatten().astype(numpy.int32)+params.offset, data.shape[0], data.shape[1],
                                   "", "%s_%.6d.cbf"%(prefix, i+1),
                                   pilatus_header="""
# Detector: %(detector)s
# %(datestr)s
# Pixel_size %(px).1e m x %(px).1e m
# Exposure_time %(exp_time)s s
# Exposure_period %(time_step)f s
# Wavelength %(wavelen)f A
# Detector_distance %(clen)f m
# Beam_xy (%(beamx).2f, %(beamy).2f) pixels
# Start_angle %(start_angle)s deg.
# Angle_increment %(alpha_step)f deg.
""" % locals(),
                                   header_convention="PILATUS_1.2")
# write_minicbf()

def run(file_in, params):
    prefix = os.path.splitext(os.path.basename(file_in))[0]
    
    h = h5py.File(file_in, "r")
    image_path = h["/Data/Image"]
    assert len(image_path.keys()) == 1
    k = image_path.keys()[0]
    data = image_path[k]["Data"]
    metadata = get_metadata(image_path[k]["Metadata"])

    open("%s_metadata.json"%prefix, "w").write(json.dumps(metadata, indent=2))
    #for i, m in enumerate(metadata): open("%s_metadata_%.6d.json"%(prefix, i+1), "w").write(json.dumps(m, indent=2))
    valid_range, mean_alpha_step, mean_time_step = dump_and_analyse_time_angle(metadata, prefix, params.detector)

    frame_range = [0, len(metadata)-1]
    
    if params.tilted_only:
        frame_range = valid_range
    if params.frame_range:
        frame_range = map(lambda x:x-1, params.frame_range)
        
    write_minicbf(data, metadata, prefix, frame_range, mean_alpha_step, mean_time_step, params)

    
if __name__ == "__main__":
    import sys

    if not sys.argv[1:] or "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
        print """\
This script extracts metadata json file, writes dat file for time and alpha angles for each frame, and converts frames to cbf files.
Only tested with Velox-written files and the BM-Ceta detector.

Usage: convert_emd.py foo.emd [options]

Parameters:"""
        iotbx.phil.parse(master_params_str).show(prefix="  ", attributes_level=1)
        quit()
    
    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str
                                              )
    params = cmdline.work.extract()
    args = cmdline.remaining_args
    file_in = args[0]
    run(file_in, params)
