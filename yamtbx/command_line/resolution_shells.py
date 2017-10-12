# LIBTBX_SET_DISPATCHER_NAME yamtbx.resolution_shells

import numpy

def run(d_max, d_min, nbins, power, quiet=False):
    start = 1./(d_max**power)
    step = ( d_min**(-power) - start ) / float(nbins)
    d_vals = 1./(start + numpy.array(range(nbins+1)) * step)**(1./power)

    if not quiet:
        print "%d resolution shells (%.3f - %.3f A) split by 1/d^%d" % (nbins, d_max, d_min, power)
        print " ".join(map(lambda x: "%.3f"%x, d_vals))
        print
        print "For XSCALE,"
        print " RESOLUTION_SHELLS= %s" % (" ".join(map(lambda x: "%.3f"%x, d_vals[1:])))

    return d_vals
    
# run()

if __name__ == "__main__":
    import sys
    import optparse

    parser = optparse.OptionParser(prog="yamtbx.resolution_shells",
                                   description="Show resolution shells",
                                   usage="usage: %prog [options] d_max d_min")

    parser.add_option("-n","--nshells", action="store", dest="nbins", type=int, default=9,
                      help="Number of shells (default: 9)")
    parser.add_option("-p","--power", action="store", dest="pow", type=int, default=2,
                      help="Split shells by 1/d^power. 2: xds style (default); 3: scalepack style")

    opts, args = parser.parse_args(sys.argv[1:])

    if len(args) != 2:
        parser.print_help()
        quit()

    try:
        d_max, d_min = map(float, args)
    except:
        parser.print_help()
        quit()

    if d_max < d_min: d_max, d_min = d_min, d_max

    run(d_max, d_min, opts.nbins, opts.pow)
