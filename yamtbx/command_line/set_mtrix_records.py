# LIBTBX_SET_DISPATCHER_NAME yamtbx.set_mtrix_records
from __future__ import absolute_import, division, print_function, generators, unicode_literals
import subprocess
import numpy
import re

re_cn = re.compile("^([CD])([0-9]+)$")

def read_cryst1(pdb_in):
    for l in open(pdb_in):
        if l.startswith("CRYST1"):
            try:
                cell = list(map(float, (l[6:15], l[15:24], l[24:33],
                                        l[33:40], l[40:47], l[47:54])))
                return cell
            except:
                print("Invalid CRYST1 card:\n%s"%l.rstrip())
    return None
# read_cryst1()

def mtrix_card(idx, M, transl):
    pdb_mtrix = ""
    for j in range(3):
        pdb_mtrix += "MTRIX%d %3d%10.6f%10.6f%10.6f     %10.5f    %s\n" % (j+1, idx, M[j,0],M[j,1],M[j,2], transl[j], 1 if idx==1 else " ")
    return pdb_mtrix

def read_shifts(file_in):
    for l in open(file_in):
        if l.startswith("pdbout cell"):
            cell = numpy.array(list(map(float, l.split()[2:])))
        elif "pdbin shifts" in l:
            shifts = numpy.array(list(map(float, l.split()[2:])))

    return cell, shifts
# def read_shifts()

def generate_matrices(symbol, rot):
    ret = []
    for i in range(rot):
        deg = 360./rot*i
        t = numpy.deg2rad(deg)
        M = numpy.array([[numpy.cos(t), -numpy.sin(t), 0], [numpy.sin(t), numpy.cos(t), 0], [0,0,1]])
        ret.append(M)
        if symbol == "D":
            M = numpy.dot(numpy.array([[-1.,0,0],[0,1,0],[0,0,-1.]]), M)
            ret.append(M)

    return ret
# generate_matrices()

def get_matrices_using_relion(sym):
    p = subprocess.Popen(["relion_refine", "--sym", sym, "--print_symmetry_ops"],
                         stdout=subprocess.PIPE)
    ret = []
    read_flag = -1
    for l in p.stdout:
        if "R(" in l:
            ret.append(numpy.zeros((3,3)))
            read_flag = 0
        elif 0 <= read_flag < 3:
            ret[-1][read_flag,:] = map(float, l.split())
            read_flag += 1
        elif read_flag >= 3:
            read_flag = -1
    return ret
# get_matrices_using_relion()

if __name__ == "__main__":
    import sys
    import optparse

    parser = optparse.OptionParser(usage="usage: %prog [options] pdbfile symmetry")
    parser.add_option("--box-size", action="store", dest="box_size", type=float, help="box size in angstrom (unit cell length)")
    parser.add_option("--output","-o", action="store", dest="output", help="Output filename")
    parser.add_option("--overwrite", action="store_true", dest="overwrite", help="Overwrite input file.")
    parser.add_option("--shifts", action="store", dest="shifts", help="REFMAC's shifts.txt")
    opts, args = parser.parse_args(sys.argv[1:])

    if len(args) < 2:
        parser.print_help()
        quit()
    
    if opts.overwrite and opts.output:
        print("Error: both --overwrite and --output were specified.")
        quit()

    pdb, symm = args
    r = re_cn.search(symm)
    Ms = []
    if r:
        symbol, rot = r.group(1), int(r.group(2))
        if rot < 1:
            print("Invalid rotation requested: %d" % rot)
            quit()

        print("rotation= %d (%.3f deg)"%(rot, 360./float(rot)))
        Ms = generate_matrices(symbol, rot)
    else:
        print("Trying to obtain matrices using RELION..")
        Ms = get_matrices_using_relion(symm)

    if opts.box_size is None:
        cell = read_cryst1(pdb)
        if cell is None or not numpy.allclose(cell[3:6], 90):
            print("CRYST1 does not look orthorhombic: %s" % cell)
            print("Please give --box-size= in angstrom")
            quit()
    else:
        a = opts.box_size
        cell = [a,a,a,90,90,90]
        
    print("cell parameter= %s"% cell)
    center = numpy.array(cell[:3])/2.

    if opts.shifts:
        org_cell, shift = read_shifts(opts.shifts)
        print("original cell= %s"%org_cell)
        print("shift= %s"%shift)
        assert numpy.allclose(org_cell[3:6], 90)
        center = org_cell[:3]/2 + shift

    print("center= %s"% center)

    to_ncsc = lambda m,t: "ncsc matrix %s %s"%(" ".join(map(str, m.flatten()))," ".join(map(str, t)))
    
    pdb_mtrix = ""
    idx = 0
    ncsc = []

    for i, M in enumerate(Ms):
        transl = numpy.dot(M, -center) + center
        pdb_mtrix += mtrix_card(i+1, M, transl)
        ncsc.append(to_ncsc(M, transl))
        
    print("Refmac instructions:")
    print("\n".join(ncsc))
            
    pdb_lines = open(pdb).readlines()
    i_insert = 0
    for i, l in enumerate(pdb_lines):
        if l.startswith(("ATOM","HETATM")):
            i_insert = i
            break

    if not opts.overwrite and not opts.output:
        quit()
        
    if opts.overwrite: opts.output = pdb
    ofs = open(opts.output, "w")
    for i, l in enumerate(pdb_lines):
        if i_insert == i:
            ofs.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1            \n"%tuple(cell))
            ofs.write(pdb_mtrix)
        if l.startswith(("CRYST1", "MTRIX", "ORIGX", "SCALE")):
            continue
        
        ofs.write(l)
    
