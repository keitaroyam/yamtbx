from __future__ import print_function
from __future__ import unicode_literals
import numpy
import cctbx.eltbx.xray_scattering
import cctbx.eltbx.e_scattering

def fetch_equation(table):
    return "+".join(["%f*exp(-%f*s**2)" %(a,b) for a,b in zip(table.array_of_a(), table.array_of_b())]) + "+%f" % table.c()

def run(elements, smin=0, smax=1, sstep=0.01):
    #reg = cctbx.xray.scattering_type_registry()
    #reg.process(scatterers=flex.xray_scatterer([cctbx.xray.scatterer("S")]))
    #reg.assign_from_table("PENG1996")
    #reg.assign_from_table("IT1992")
    #print reg.unique_form_factors_at_d_star_sq(0.05**2)[0]

    print("element s xray electron")
    for el in elements:
        xray = cctbx.eltbx.xray_scattering.it1992(el, True).fetch()
        elec = cctbx.eltbx.e_scattering.ito_vol_c_2011_table_4_3_2_2_entry_as_gaussian(label=el, exact=True)

        print("# Xray for %s    :"%el, fetch_equation(xray))
        print("# electron for %s:"%el, fetch_equation(elec))
    
        for s in numpy.arange(smin, smax, sstep):
            print("%2s %.4f %.4f %.4f" % (el, s, xray.at_d_star(s), elec.at_d_star(s)))

if __name__ == "__main__":
    import sys
    elements = sys.argv[1:]
    run(elements)
