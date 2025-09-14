from __future__ import print_function
from __future__ import unicode_literals
import os
import re

class XtriageLogfile(object):
    def __init__(self, logfile):
        self.logfile = logfile
        self.wilson_b = float("nan")
        self.anisotropy = float("nan")

        self.read_logfile()
    # __init__()

    def read_logfile(self):
        if not os.path.isfile(self.logfile):
            print("Can't open file:", self.logfile)
            return

        reading = ""
        Bs = []
        for l in open(self.logfile):
            if "ML estimate of overall B value" in l:
                reading = "wilsonb"
            elif reading == "wilsonb":
                self.wilson_b = float(l.split()[0])
                reading = ""
            elif "Eigen analyses of B-cart:" in l:
                reading = "B_cart"
            elif l.strip() == "":
                reading = ""
            elif reading == "B_cart":
                r = re.search(r"^  \| ([123]) *\| *([-0-9\.]*) *\| ",l)
                if r:
                    Bs.append(float(r.group(2)))
                    if int(r.group(1)) == 3:
                        reading = ""

        if Bs:
            self.anisotropy = max(Bs) - min(Bs)
    # read_logfile()
# class XtriageLogfile

if __name__ == "__main__":
    import sys

    xlf = XtriageLogfile(sys.argv[1])
    print(xlf.wilson_b, xlf.anisotropy)
