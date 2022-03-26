"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals
from yamtbx.dataproc.xds import xscalelp
import os

def run(lpin):
    prefix = os.path.basename(lpin)

    data_read = xscalelp.get_read_data(lpin) # i, <I>, Naccepted, Nrejected, filename
    k_b = xscalelp.get_k_b(lpin) # K, B, filename
    ISas = xscalelp.get_ISa(lpin) # a, b, ISa, ISa0, filename
    rfactors = xscalelp.get_rfactors_for_each(lpin) # {filename: list of [dmin, Robs, Rexpt, Compared]}

    # Transform them to dict(filename=data)
    data_read = dict([(x[4], x[:4]) for x in data_read])
    k_b = dict([(x[2], x[:2]) for x in k_b])
    ISas = dict([(x[4], x[:4]) for x in ISas])

    ofs = open(prefix+"_table.dat", "w")
    print("#", os.path.abspath(lpin), file=ofs)
    print("file I.mean n.accept n.reject K B em.a em.b ISa ISa0 R", file=ofs)
    nan = float("nan")
    maxlen = str(max([len(x) for x in data_read]))

    for f in data_read:
        Imean, nacc, nrej = data_read[f][1:]
        K, B = k_b.get(f, [nan,nan])
        a, b, ISa, ISa0 = ISas.get(f, [nan]*4)
        r_total = rfactors[f][-1][1] if f in rfactors else None
        data = f, Imean, nacc, nrej, K, B, a, b, ISa, ISa0, r_total
        print(("%-"+maxlen+"s %.4e %6d %6d %7.3f %7.3f %.3e %.3e %5.2f %5.2f %7.1f") % data, file=ofs)

    ofs = open(prefix+"_rfactors.dat", "w")
    print("#", os.path.abspath(lpin), file=ofs)
    print("file dmin R.obs R.exp comp", file=ofs)
    for f in rfactors:
        for dmin, robs, rexp, comp in rfactors[f]:
            if dmin is not None:
                print(("%-"+maxlen+"s %5.2f %7.1f %7.1f %5d") % (f, dmin, robs, rexp, comp), file=ofs)
# run()

if __name__ == "__main__":
    import sys
    run(sys.argv[1]) # give pointless log
