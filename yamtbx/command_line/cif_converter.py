import os
import iotbx.cif
import iotbx.shelx.hklf

def run(cif_file, prefix=None):
    if prefix is None: prefix = os.path.splitext(os.path.basename(cif_file))[0]
    
    reader = iotbx.cif.reader(cif_file)
    print reader.build_miller_arrays()
    print reader.as_miller_arrays()

    model = reader.model()
    hklf = None
    
    xrss = reader.build_crystal_structures()

    for key in xrss:
        xrs = xrss[key]
        print key, xrs
        print "Saving %s_%s.pdb"%(prefix,key)
        open("%s_%s.pdb"%(prefix,key), "w").write(xrs.as_pdb_file(resname="XXX"))
    
    for k in model:
        if "_shelx_res_file" in model[k]:
            res_str = model[k]["_shelx_res_file"].lstrip("\r\n")
            open("%s_%s.res"%(prefix,k), "w").write(res_str)

            res_lines = res_str.splitlines()
            hklf_line = filter(lambda x: x.startswith("HKLF "), res_lines)
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

                
                mtz_ds.mtz_object().write(hklout+".mtz")
                

if __name__ == "__main__":
    import sys
    cif_file = sys.argv[1]
    run(cif_file)
