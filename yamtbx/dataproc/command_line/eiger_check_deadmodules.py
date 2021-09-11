from __future__ import print_function
from __future__ import unicode_literals
import h5py
import numpy
#from yamtbx.dataproc import cbf

def get_module_info(h):
    ret = []
    dsp = h["/entry/instrument/detector/detectorSpecific/"]
    for k in sorted([x for x in list(dsp.keys()) if x.startswith("detectorModule_")]):
        rot = dsp[k]["data_rotation"][()]
        org = dsp[k]["data_origin"][:]
        siz = dsp[k]["data_size"][:]

        assert rot in (0, 180)
        print("#",k, org, siz, rot)
        if rot == 180:
            org = org - siz

        ret.append((k, org, siz))
    return ret

def data_iterator(h):
    for k in sorted(h["/entry/data"].keys()):
        if not h["/entry/data"].get(k): continue
        for i in range(h["/entry/data"][k].shape[0]):
            yield h["/entry/data"][k][i]

def numbers_to_range(n):
    ret = "%d" % n[0]

    for i in range(1, len(n)):
        if n[i]-n[i-1] <= 1:
            if ret[-1] != "-": ret += "-"
            if i==len(n)-1: ret += "%d"%n[i]
        else:
            if ret[-1] == "-":
                ret += "%d"%n[i-1]
                if i<len(n): ret += ",%d"%n[i]
            else: ret += ",%d" % n[i]
    return ret
    
def run(masterf):
    print("# %s" % masterf)
    h = h5py.File(masterf, "r")
    modules = get_module_info(h)
    print("frame module num.bad percent.bad")
    ret = {}
    for fn, data in enumerate(data_iterator(h)):
        #if fn<10000: continue
        #tmp = numpy.zeros(data.shape, dtype=numpy.int32)
        for name, origin, size in modules:
            #tmp[origin[1]:origin[1]+size[1],
            #    origin[0]:origin[0]+size[0]] += j+1

            d = data[origin[1]:origin[1]+size[1],
                     origin[0]:origin[0]+size[0]]
            
            bad_value = 2**(d.dtype.itemsize*8)-1
            num_bad = numpy.sum(d == bad_value)
            frac = num_bad/float(d.size)
            
            if frac>0.5:
                print("%6d %s %6d %5.1f" % (fn+1, name, num_bad, frac*100.))
                ret.setdefault(name, []).append(fn+1)
        
        #cbf.save_numpy_data_as_cbf(tmp.flatten(),
        #                           size1=tmp.shape[1], size2=tmp.shape[0],
        #                           title="test",
        #                           cbfout="junk.cbf")

    if ret:
        print("# Problematic modules:")
        for name in sorted(ret):
            print("#  %s %d frames (%s)" % (name, len(ret[name]), numbers_to_range(ret[name])))

if __name__ == "__main__":
    import sys
    run(sys.argv[1])
