from __future__ import print_function
import numpy

class InitLp(object):
    def __init__(self, initlp):
        self.initlp = initlp
        self.parse()

    def parse(self):
        self.bkginit_table = []
        flag = None
        for l in open(self.initlp):
            if l.startswith(" FRAME #    GAIN     MGAIN"):
                flag = "gain"
            elif l.startswith(" FRAME #    SCALE     COUNTS     MBKG    MSCALE      MSPOT"):
                flag = "bkginit"
            elif l.startswith("  FRAME#     SCALE     BACKGND       GAIN     #PIXEL"):
                # From ver Version June 30, 2024
                # number of items has changed, but we only need first and third items..
                flag = "bkginit2"
            elif flag == "gain":
                try:
                    frame, gain, mgain = l.split()
                except:
                    continue
            elif flag == "bkginit":
                try:
                    frame, scale, counts, mbkg, mscale, mspot = l.split()
                    frame, mbkg, mscale, mspot = list(map(int, (frame, mbkg, mscale, mspot)))
                    scale, counts = list(map(float, (scale, counts)))
                    self.bkginit_table.append([frame, scale, counts, mbkg, mscale, mspot])
                except:
                    continue
            elif flag == "bkginit2":
                try:
                    frame, scale, counts, gain, npix = l.split()
                    self.bkginit_table.append([int(frame), float(scale), float(counts), float(gain), int(npix)])
                except:
                    continue            
            elif "NUMBER OF GAIN VALUES IN TABLE" in l:
                flag = None
            elif "NUMBER OF PIXELS USED FOR DATA COLLECTION" in l:
                flag = None
    # parse()        

    def check_bad_first_frames(self, iqrc=1.5):
        # /isilon/users/rikenbl/rikenbl/Staff_26B2/ito/HSA/20190730/_kamoproc/cps1998/w02/data/multi_1-100
        counts = numpy.array([x[2] for x in self.bkginit_table])
        frames = [x[0] for x in self.bkginit_table]
        q25, q50, q75 = numpy.percentile(counts, [25, 50, 75])
        iqr = q75 - q25
        lowlim, highlim = q25 - iqrc*iqr, q75 + iqrc*iqr
        bad_sel = (counts < lowlim)
        bad_first_frames = []
        if bad_sel[0]:
            print("Bad first frames for", self.initlp)
            last_i = -1
            for i in numpy.where(counts < lowlim)[0]:
                if i - last_i > 1: break
                last_i = i
                print("", frames[i], counts[i])
                bad_first_frames.append(frames[i])

        return bad_first_frames
    # check_bad_first_frames()

if __name__ == "__main__":
    import sys
    lp = InitLp(sys.argv[1])
    lp.check_bad_first_frames()
