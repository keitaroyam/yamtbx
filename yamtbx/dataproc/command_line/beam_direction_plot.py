import numpy
import os
from cctbx import sgtbx
from cctbx import crystal
from cctbx.crystal import reindex
from yamtbx.dataproc.xds import xparm
from yamtbx.dataproc.xds import xds_ascii
from yamtbx.util.xtal import abc_convert_real_reciprocal
from yamtbx.util import read_path_list
import iotbx.phil
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt


master_params_str = """
input = None
 .type = path
 .multiple = true

ref_cell = None
 .type = floats(size=6)
ref_symm = None
 .type = str

plot_out = "beam_plot.pdf"
 .type = path
dat_out = None
 .type = path
"""

def get_angles(ib, a_axis, b_axis, c_axis, xs, ref_xs, filename):
    ib = ib / numpy.linalg.norm(ib)

    if ref_xs:
        cosets = reindex.reindexing_operators(ref_xs, xs, 0.2, 20)
        if len(cosets.combined_cb_ops())==0:
            print "Error: cannot find reindex operator"
            return
        
        op = cosets.combined_cb_ops()[0]
        m = numpy.array(op.c_inv().r().as_double()).reshape(3,3).transpose()
        transformed = numpy.dot(m, numpy.array([a_axis,b_axis,c_axis]))
        a, b, c = transformed[0,:], transformed[1,:], transformed[2,:]
        sg = ref_xs.space_group()
        print "Reading %s (%s; reindex: %s)" % (filename, sg.info(), op)
    else:
        a, b, c = a_axis, b_axis, c_axis
        sg = xs.space_group()
        print "Reading %s (%s)" % (filename, sg.info())

    laue = sg.build_derived_reflection_intensity_group(False).build_derived_point_group()

    astar, bstar, cstar = abc_convert_real_reciprocal(a,b,c)
    bstar_e = bstar / numpy.linalg.norm(bstar)
    cstar_e = cstar / numpy.linalg.norm(cstar)

    # incident_beam along a,b,c
    xyz = (numpy.dot(ib, a)/numpy.linalg.norm(a),
           numpy.dot(ib, b)/numpy.linalg.norm(b),
           numpy.dot(ib, c)/numpy.linalg.norm(c))
    angles = []

    #print "Replicating vectors using %2d operators" % len(laue.all_ops())

    # replicate vectors
    for op in laue.all_ops():
        vec_elem = op * xyz
        vec = vec_elem[0] * a + vec_elem[1] * b + vec_elem[2] * c
        veclen = numpy.linalg.norm(vec)

        theta = numpy.arccos(numpy.dot(vec, cstar_e)/veclen)
        if theta > numpy.pi/2: continue
        
        vec_on_ab = vec - cstar_e * veclen*numpy.cos(theta)
        phi = numpy.arccos(numpy.dot(vec_on_ab, a)/numpy.linalg.norm(a)/numpy.linalg.norm(vec_on_ab))

        if numpy.dot(vec_on_ab, bstar_e) < 0: phi = - phi

        angles.append((theta, phi))

    return angles
# get_angles()
    

def from_xds_ascii(xacin, ref_xs):
    xac = xds_ascii.XDS_ASCII(xacin, read_data=False)
    return get_angles(xac.incident_axis,
                      numpy.array(xac.a_axis), numpy.array(xac.b_axis), numpy.array(xac.c_axis),
                      xac.symm,
                      ref_xs,
                      xacin)
# from_xds_ascii()

def from_xparm(xpin, ref_xs):
    xp = xparm.XPARM(xpin)
    return get_angles(xp.incident_beam, 
                      xp.a_axis, xp.b_axis, xp.c_axis,
                      xp.crystal_symmetry(),
                      ref_xs,
                      xpin)
# from_xparm()

def from_crystfel_stream(stin, ref_xs): # XXX Not tested!!
    from yamtbx.dataproc.crystfel import stream
    ret = []
    ib = numpy.array([0., 0., 1.])
    for chunk in stream.stream_iterator(stin, read_reflections=False):
        xs = chunk.indexed_symmetry()
        if xs is None: continue
        ub = chunk.ub_matrix()
        avec, bvec, cvec = numpy.linalg.inv(ub)
        ret.extend(get_angles(ib, avec, bvec, cvec,
                              xs, ref_xs, "%s:%s:%s" % (stin, chunk.filename, chunk.event)))

    return ret
# from_crystfel_stream()
                              
def make_dat(angles, dat_out):
    ofs = open(dat_out, "w")
    ofs.write("theta phi\n")
    for theta, phi in angles: ofs.write("%.4f %.4f\n" % (theta,phi))
    ofs.close()
    print
    print "Data file written. See %s" % dat_out
    print """\
If you want to plot it with R:
R
library(ggplot2)
d <- read.table("%s",h=T)
p <- ggplot(d, aes(x=phi, y=2*sin(theta/2))) +geom_point(alpha=.1,size=2) +coord_polar(start=pi/2,direction=-1) +scale_x_continuous(breaks=seq(-pi,pi,by=20/180*pi),labels=function(x)sprintf("%%.0f",x/pi*180)) 
ggsave("beam.png", p)
""" % dat_out
# make_dat()

def make_plot(angles, plot_out):
    r = map(lambda x: 2*numpy.sin(x[0]/2.), angles)
    phi = map(lambda x: x[1], angles)

    ax = plt.subplot(111, projection='polar')
    ax.scatter(phi, r, color="b", alpha=.1, s=4)
    ax.set_rmax(numpy.sqrt(2.))
    
    theta_labs = numpy.linspace(0, numpy.pi/2, 5, endpoint=True)
    ax.set_yticks(2*numpy.sin(theta_labs/2.))
    ax.set_yticklabels(map(lambda x: "%.0f$^\circ$"%x, numpy.rad2deg(theta_labs)))
    ax.grid(True)

    ax.set_title("Beam direction plot", va='bottom')
    plt.savefig(plot_out, dpi=150)
    plt.gcf().clear() # Required when plot more than once
    print
    print "Plot written. See %s" % plot_out
# make_plot()

def run(params):
    ref_xs = None

    if None not in (params.ref_cell, params.ref_symm):
        ref_xs = crystal.symmetry(params.ref_cell, params.ref_symm)

    if len(params.input) == 1 and params.input[0].endswith(".lst"):
        params.input = read_path_list(params.input[0])

    if len(params.input) == 0: return

    angles = []

    for arg in params.input:
        if ".stream" in arg:            
            angles.extend(from_crystfel_stream(arg, ref_xs))
        elif xds_ascii.is_xds_ascii(arg):
            angles.extend(from_xds_ascii(arg, ref_xs))
        else:
            angles.extend(from_xparm(arg, ref_xs))

    if params.dat_out: make_dat(angles, params.dat_out)
    if params.plot_out: make_plot(angles, params.plot_out)
# run()

def run_from_args(argv):
    if "-h" in argv or "--help" in argv or not argv:
        print """
This program makes a beam-direction plot in crystal frame to see bias of crystal orientations.
The beam direction is replicated using crystal symmetry.

Usages:
  yamtbx.beam_direction_plot formerge.lst # the list of (G)XPARM.XDS or XDS_ASCII.HKL.
  yamtbx.beam_direction_plot formerge.lst ref_cell=10,20,30,90,90,90 ref_symm=p222 # In case you need change of basis

Theta=0 at c*-axis and Phi=0 at a-axis.

Parameters:\
"""
        iotbx.phil.parse(master_params_str).show(prefix="  ", attributes_level=1)
        return

    cmdline = iotbx.phil.process_command_line(args=argv,
                                              master_string=master_params_str)
    params = cmdline.work.extract()

    for arg in cmdline.remaining_args:
        if os.path.isfile(arg):
            params.input.append(arg)
        else:
            print "File not found:",arg
            return

    run(params)

# run_from_args()
if __name__ == "__main__":
    import sys
    run_from_args(sys.argv[1:])
