from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import numpy
from yamtbx.dataproc.xds.xparm import XPARM
from yamtbx.dataproc import cbf

"""
From http://xds.mpimf-heidelberg.mpg.de/html_doc/coordinate_systems.html:

ORG(:) = -ORGX*QX*ED(:,1)-ORGY*QY*ED(:,2)+F*ED(:,3)

DIRECTION_OF_SEGMENT_X-AXIS=EDS(1,1) EDS(2,1) EDS(3,1)
DIRECTION_OF_SEGMENT_Y-AXIS=EDS(1,2) EDS(2,2) EDS(3,2)
EDS(:,3)=EDS(:,1) X EDS(:,2) to form a right-handed orthonormal segment system {EDS(:,1), EDS(:,2), EDS(:,3)}.
The origin of the segment system: -ORGXS*QX*EDS(:,1)-ORGYS*QY*EDS(:,2)+FS*EDS(:,3)
where SEGMENT_ORGX= ORGXS, SEGMENT_ORGY= ORGYS, and SEGMENT_DISTANCE= FS.

The representation of the segment system EDS with respect to the laboratory system ED is the matrix product EDSL = MATMUL(ED,EDS).
Thus, a segment pixel at IX,IY has the laboratory coordinates (mm units)
x=QX*(IX-ORGXS)*EDSL(1,1)+QY*(IY-ORGYS)*EDSL(1,2)+FS*EDSL(1,3)+ORG(1)
y=QX*(IX-ORGXS)*EDSL(2,1)+QY*(IY-ORGYS)*EDSL(2,2)+FS*EDSL(2,3)+ORG(2)
z=QX*(IX-ORGXS)*EDSL(3,1)+QY*(IY-ORGYS)*EDSL(3,2)+FS*EDSL(3,3)+ORG(3)
"""

def write_atom(ofs, xyz, iseg):
    s = xyz.shape[:2]
    n = 0
    for i in (0, s[0]-1):
        for j in (0, s[1]-1):
            x,y,z = xyz[i,j,:]
            n+=1
            ofs.write("ATOM      1  C%d  XXX A%4d    %8.3f%8.3f%8.3f  1.00%6.2f           C\n"%(n,iseg,x,y,z,0))

def run(xds_inp):
    xp = XPARM()
    xp.set_info_from_xdsinp_or_inpstr(xdsinp=xds_inp)

    d_map = numpy.zeros((xp.ny, xp.nx)) # resolution mapping

    ed = numpy.zeros((3,3))
    ed[:,0] = xp.X_axis
    ed[:,1] = xp.Y_axis
    ed[:,2] = numpy.cross(xp.X_axis, xp.Y_axis)
    ed /= numpy.linalg.norm(ed, axis=0)


    orgx, orgy = xp.origin
    fs = xp.distance
    qx, qy = xp.qx, xp.qy
    org = -orgx*qx*ed[:,0]-orgy*qy*ed[:,1]+fs*ed[:,2]
    wavelen = xp.wavelength
    s0 = xp.incident_beam/numpy.linalg.norm(xp.incident_beam)/wavelen

    print("qx,qy=", qx, qy)
    print("s0=", s0)
    print("lambda=", wavelen)
    print("ORG=", org)
    print("ED=")
    print(ed)

    ofs_pdb = open("detector_pos.pdb", "w")
    ofs_pml = open("for_pymol.pml", "w")
    
    for iseg, seg in enumerate(xp.segments):
        print("Segment %d"%(iseg+1))
        eds = numpy.zeros((3,3))
        eds[:,0] = seg.eds_x
        eds[:,1] = seg.eds_y
        eds[:,2] = numpy.cross(seg.eds_x, seg.eds_y)
        eds /= numpy.linalg.norm(eds, axis=0)
        edsl = numpy.dot(ed, eds)

        ix, iy = numpy.meshgrid(list(range(seg.x1-1,seg.x2)), list(range(seg.y1-1, seg.y2)))
        tmp = numpy.zeros((ix.shape[0], ix.shape[1], 3))
        for i in range(3):
            tmp[:,:,i] = qx*(ix-seg.orgxs)*edsl[i,0]+qy*(iy-seg.orgys)*edsl[i,1]+seg.fs*edsl[i,2]+org[i]

        write_atom(ofs_pdb, tmp, iseg)
        ofs_pml.write("bond resi %d and name C1, resi %d and name C2\n"%(iseg,iseg))
        ofs_pml.write("bond resi %d and name C2, resi %d and name C4\n"%(iseg,iseg))
        ofs_pml.write("bond resi %d and name C3, resi %d and name C4\n"%(iseg,iseg))
        ofs_pml.write("bond resi %d and name C3, resi %d and name C1\n"%(iseg,iseg))

        tmpdenom = numpy.linalg.norm(tmp,axis=2) * wavelen
        for i in range(3):
            tmp[:,:,i] /= tmpdenom

        s1 = tmp
        s = s1-s0
        d_map[iy,ix] = 1./numpy.linalg.norm(s, axis=2)
        
        """
        for ix in range(seg.x1-1, seg.x2):
            for iy in range(seg.y1-1, seg.y2):
                tmp = (qx*(ix-seg.orgxs)*edsl[0,0]+qy*(iy-seg.orgys)*edsl[0,1]+seg.fs*edsl[0,2]+orgx,
                       qx*(ix-seg.orgxs)*edsl[1,0]+qy*(iy-seg.orgys)*edsl[1,1]+seg.fs*edsl[1,2]+orgy,
                       qx*(ix-seg.orgxs)*edsl[2,0]+qy*(iy-seg.orgys)*edsl[2,1]+seg.fs*edsl[2,2]+fs)
                s1 = tmp/numpy.linalg.norm(tmp)/wavelen
                s = s1-s0
                d_map[iy,ix] = numpy.linalg.norm(s)
        """

    ofs_pml.write("pseudoatom org, pos=(0,0,0)\n")
    ofs_pml.write("pseudoatom s0, pos=(%f,%f,%f)\n"%tuple(s0*100))
    ofs_pml.write("as spheres, org s0\n")
    ofs_pml.write("color red, s0\n")
    ofs_pml.write("distance selection1=org, selection2=s0\n")
    ofs_pml.write("set sphere_scale, 10\n")
        
    data = (d_map*1000).astype(numpy.int32)
    cbf.save_numpy_data_as_cbf(data.flatten(), data.shape[1], data.shape[0], "d_map", "d_map_x1000.cbf")


if __name__ == "__main__":
    import sys
    xds_inp = sys.argv[1]
    run(xds_inp)
