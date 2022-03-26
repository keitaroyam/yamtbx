from __future__ import print_function
from __future__ import unicode_literals
import iotbx.mtz
from cctbx.array_family import flex

def get_I(mtzin):
    mtzobj = iotbx.mtz.object(file_name=mtzin)
    I_arrays = [s for s in mtzobj.as_miller_arrays() if s.info().type_hints_from_file =="intensity"]
    F_arrays = [s for s in mtzobj.as_miller_arrays() if s.info().type_hints_from_file =="amplitude"]

    if len(I_arrays) > 0:
        return I_arrays[0]
    elif len(F_arrays) > 0:
        return F_arrays[0].as_intensity_array()

if __name__ == "__main__":
    import sys

    data = {}

    for f in sys.argv[1:]:
        data[f] = get_I(f)

    for ix in range(len(data)-1):
        for iy in range(ix+1, len(data)):
            x, y = list(data.keys())[ix], list(data.keys())[iy]
            xd, yd = data[x].common_sets(data[y], assert_is_similar_symmetry=False)
            corr = flex.linear_correlation(xd.data(), yd.data())
            assert corr.is_well_defined()
            print(x, "vs", y, " cc=", corr.coefficient())

