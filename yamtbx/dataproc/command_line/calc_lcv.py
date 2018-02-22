from yamtbx.dataproc import blend_lcv
import numpy
import os

def run(filein):
    cells = []
    files = []
    
    for l in open(filein):
        sp = l.split()
        if len(sp) == 6:
            try:
                sp = map(float, sp)
                cells.append(sp)
                continue
            except:
                pass
            
        elif "INPUT_FILE=" in l:
            tmp = l.split("=")[1].strip()
            if tmp[0] == "*":
                tmp = tmp[1:].strip()
            files.append(tmp)
        elif len(sp) == 1:
            tmp = sp[0].strip()
            files.append(tmp)

    for f in files:
        if not os.path.isabs(f):
            f = os.path.join(os.path.dirname(filein), f)
        
        if not os.path.isfile(f):
            print "not exists:", f
            continue
        for l in open(f):
            if l.startswith("!SPACE_GROUP_NUMBER="):
                sg = l[l.index("=")+1:].strip()
            if l.startswith("!UNIT_CELL_CONSTANTS="):
                cell = map(float, l[l.index("=")+1:].split())
                assert len(cell) == 6
                cells.append(cell)
                break

    print "%d cells loaded." % len(cells)
    cells = numpy.array(cells)
    mean_cell = map(lambda i: cells[:,i].mean(), xrange(6))
    cell_std = map(lambda i: numpy.std(cells[:,i]), xrange(6))
    mean_cell_str = " ".join(map(lambda x:"%.3f"%x, mean_cell))
    cell_std_str = " ".join(map(lambda x:"%.1e"%x, cell_std))
    print "Averaged cell= %s" % mean_cell_str
    print "Cell Std dev.= %s" % cell_std_str
    print ""
    lcv, alcv = blend_lcv.calc_lcv(cells)
    print " LCV= %.2f%%" % lcv
    print "aLCV= %.2f A" % alcv


if __name__ == "__main__":
    import sys
    run(sys.argv[1])
