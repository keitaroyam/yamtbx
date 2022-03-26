from __future__ import print_function
from __future__ import unicode_literals
import os
import iotbx.cif
import iotbx.shelx.hklf

def write_mtz(arrays, prefix):
    arrays_for_mtz = []
    for arr in arrays:
        labs = arr.info().labels
        arr = arr.map_to_asu()
        if "_refln_F_squared_meas" in labs:
            print("Iobs found:", labs)
            if not arr.is_unique_set_under_symmetry():
                merge = arr.merge_equivalents(use_internal_variance=False)
                merge.show_summary()
                arr = merge.array()
            arrays_for_mtz.append(("I", "I"+labs[0], arr))

            if arr.anomalous_flag():
                arrays_for_mtz.append(("I", "IMEAN", arr.average_bijvoet_mates()))
        elif "_refln_F_squared_calc" in labs:
            print("Icalc found:", labs)
            arrays_for_mtz.append(("I", "IC", arr))
        else:
            print("skipping", labs)

    if not arrays_for_mtz: return
    mtz_ds = arrays_for_mtz[0][2].as_mtz_dataset(arrays_for_mtz[0][1])
    for arr in arrays_for_mtz[1:]:
        mtz_ds.add_miller_array(arr[2], arr[1])

    print("Writing", prefix+".mtz")
    mtz_ds.mtz_object().write(prefix+".mtz")
# write_mtz()

def run(cif_file, prefix=None):
    if prefix is None: prefix = os.path.splitext(os.path.basename(cif_file))[0]
    
    reader = iotbx.cif.reader(cif_file)
    print(reader.build_miller_arrays())
    print(reader.as_miller_arrays())

    model = reader.model()
    hklf = None
    
    xrss = reader.build_crystal_structures()

    for key in xrss:
        xrs = xrss[key]
        print(key, xrs)
        print("Saving %s_%s.pdb"%(prefix,key))
        open("%s_%s.pdb"%(prefix,key), "w").write(xrs.as_pdb_file(resname="XXX"))

    write_mtz(reader.as_miller_arrays(), prefix)
    
    for k in model:
        if "_shelx_res_file" in model[k]:
            res_str = model[k]["_shelx_res_file"].lstrip("\r\n")
            open("%s_%s.res"%(prefix,k), "w").write(res_str)

            res_lines = res_str.splitlines()
            hklf_line = [x for x in res_lines if x.startswith("HKLF ")]
            if hklf_line:
                hklf = int(hklf_line[0][5])

        if "_shelx_hkl_file" in model[k]:
            hklout = "%s_%s.hkl"%(prefix,k)
            open(hklout, "w").write(model[k]["_shelx_hkl_file"].lstrip("\r\n"))

            if hklf:
                xs = xrss[k].crystal_symmetry()
                arrays = iotbx.shelx.hklf.reader(file_name=hklout).as_miller_arrays(xs)
                if hklf == 4:
                    arrays[0].set_observation_type_xray_intensity()
                    m = arrays[0].merge_equivalents(use_internal_variance=False)
                    m.show_summary()
                    mtz_ds = m.array().as_mtz_dataset("I")
                elif hklf == 3:
                    arrays[0].set_observation_type_xray_amplitude()
                    m = arrays[0].merge_equivalents(use_internal_variance=False)
                    m.show_summary()
                    mtz_ds = m.array().as_mtz_dataset("F")
                else:
                    raise "Unknwon HKLF %d" % hklf

                if m.array().anomalous_flag():
                  mtz_ds.add_miller_array(m.array().average_bijvoet_mates(), "IMEAN")
                
                mtz_ds.mtz_object().write(hklout+".mtz")
                

if __name__ == "__main__":
    import sys
    cif_file = sys.argv[1]
    run(cif_file)
