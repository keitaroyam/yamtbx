"""
Usage:
yamtbx.python dump_anomalous_differences.py refine_4_tls_001.mtz F-model,I-obs,F-obs-filtered

Then, for example,
R
d<-read.table("ano_refine_4_tls_001.dat",h=T)
cor(d$F.model,d$F.obs.filtered)
plot(d$F.model,d$F.obs.filtered)
sum((d$F.model>0 & d$F.obs.filtered>0) | (d$F.model<0 & d$F.obs.filtered<0))/nrow(d)*100

library(ggplot2)
d$d.cut <- cut(1/d$d**2, 10)
ggplot(d, aes(x=F.model,y=F.obs.filtered))+geom_point()+facet_grid(d.cut~.)
"""
from __future__ import print_function
from __future__ import unicode_literals

import iotbx.file_reader
from cctbx.array_family import flex
import os

def commonalize(Is):
    new_Is = []
    Is0 = Is[0]
    for I in Is[1:]:
        Is0, I = Is0.common_sets(I, assert_is_similar_symmetry=False)
        new_Is.append(I)

    Is = []

    for I in new_Is:
        I = I.common_set(Is0, assert_is_similar_symmetry=False)
        assert len(Is0.data()) == len(I.data())
        Is.append(I)

    return [Is0,] + Is
# commonalize()

def run(hklin, labels):
    labels = labels.split(",")
    data = []
    arrays = [x for x in iotbx.file_reader.any_file(hklin).file_server.miller_arrays if x.anomalous_flag()]

    for lab in labels:
        tmp = [x for x in arrays if lab+"(+)" in x.info().label_string()]

        if len(tmp) == 0:
            print("Such data not found in mtz file:", lab)
            print("All available anomalous labels:", [x.info().labels[0].replace("(+)", "") for x in arrays])
            return
        elif len(tmp) > 1:
            print("The label is ambiguous for label", lab)
            print([x.info().label_string() for x in tmp])
            return

        array = tmp[0]
        if array.is_complex_array():
            array = array.amplitudes().set_observation_type_xray_amplitude()

        anodiff = array.anomalous_differences()

        data.append(anodiff)

    data = commonalize(data)
    indices = data[0].indices()
    diffs = [x.data() for x in data]

    ofs = open("ano_%s.dat" % os.path.splitext(os.path.basename(hklin))[0], "w")
    print("h k l d %s" % " ".join(labels), file=ofs)
    for x in zip(indices, *diffs):
        hkl = x[0]
        d = anodiff.unit_cell().d(hkl)
        h,k,l = hkl
        print(h, k, l, d, " ".join(["%.3f"%y for y in x[1:]]), file=ofs)

if __name__ == "__main__":
    import sys
    hklin, labels = sys.argv[1:]
    run(hklin, labels)
