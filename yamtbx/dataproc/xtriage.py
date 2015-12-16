import os

class XtriageLogfile:
    def __init__(self, logfile):
        self.logfile = logfile
        self.wilson_b = float("nan")
        self.anisotropy = float("nan")

        self.read_logfile()
    # __init__()

    def read_logfile(self):
        if not os.path.isfile(self.logfile):
            print "Can't open file:", self.logfile
            return

        reading = ""
        for l in open(self.logfile):
            if l.startswith("ML estimate of overall B value of"):
                reading = "wilsonb"
            elif reading == "wilsonb":
                self.wilson_b = float(l.split()[0])
                reading = ""
            elif l.startswith("Anisotropy"):
                self.anisotropy = float(l.split()[-1])
    # read_logfile()
# class XtriageLogfile

if __name__ == "__main__":
    import sys

    xlf = XtriageLogfile(sys.argv[1])
    print xlf.wilson_b, xlf.anisotropy
