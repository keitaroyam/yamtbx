"""
(c) RIKEN 2017. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import unicode_literals
from yamtbx.dataproc import XIO

def make_geom(header, geom_out):
    s = """\
photon_energy = %(photon).4f
p0/res = %(res).1f ; 1m /x micron
p0/min_fs = 0
p0/max_fs = %(fsmax)d
p0/min_ss = 0
p0/max_ss = %(ssmax)d
p0/corner_x = %(cornerx).2f ; units are pixels
p0/corner_y = %(cornery).2f
p0/fs = -x
p0/ss = -y
p0/clen = %(clen)f
p0/adu_per_photon = 1
p0/max_adu = %(max_adu)d

rigid_group_p0 = p0
rigid_group_collection_connected = p0
rigid_group_collection_independent = p0
""" % dict(fsmax=header["Width"]-1, ssmax=header["Height"]-1,
               cornerx=header["BeamX"]/header["PixelX"], cornery=header["BeamY"]/header["PixelY"],
               clen=header["Distance"]/1000., res=1/(header["PixelX"]*1.e-3), max_adu=header["Overload"], photon=12398.4/header["Wavelength"])

    nh_sensors = (header["Width"]-487)//494+1
    nv_sensors = (header["Height"]-195)//212+1

    for i in range(nh_sensors):
        s += """
badp0v%(i)d/min_fs = %(min_fs)d
badp0v%(i)d/max_fs = %(max_fs)d
badp0v%(i)d/min_ss = %(min_ss)d
badp0v%(i)d/max_ss = %(max_ss)d
""" % dict(i=i, min_fs=494*i+487, max_fs=494*(i+1)-1, min_ss=0, max_ss=header["Height"]-1)

    for i in range(nv_sensors):
        s += """
badp0h%(i)d/min_fs = %(min_fs)d
badp0h%(i)d/max_fs = %(max_fs)d
badp0h%(i)d/min_ss = %(min_ss)d
badp0h%(i)d/max_ss = %(max_ss)d
""" % dict(i=i, min_fs=00, max_fs=header["Width"]-1, min_ss=212*i+195, max_ss=212*(i+1)-1)


    open(geom_out, "w").write(s)
# make_geom()


def run(cbfin):
    header = XIO.Image(cbfin).header
    make_geom(header, "pilatus.geom")
    
# run()

if __name__ == "__main__":
    import sys
    run(sys.argv[1])
