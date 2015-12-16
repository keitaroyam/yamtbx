"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
class XdsstatLp:
    def __init__(self, lpin):
        if lpin is not None:
            self.parse(lpin)
    # __init__()

    def parse(self, lpin):
        self.by_frame = {}
        self.by_framediff = {}

        reading = ""

        for l in open(lpin):
            if "Frame #refs #misfits  Iobs     sigma   Iobs/sigma   Peak      Corr     Rmeas  #rmeas #unique" in l:
                reading = "by_frame"
            elif reading == "by_frame":
                if "For inline graphs use a Java browser" in l:
                    reading = ""
                    continue

                if len(l) < 94: continue

                # frame, nrefs, nmisfits, iobs, sigma, ios, peak, corr, rmeas, nrmeas, nuniq, (dummy)
                sp = l[:6], l[6:13], l[13:18], l[18:28], l[28:38], l[38:48], l[48:58], l[58:68], l[68:78], l[78:85], l[85:92]
                self.by_frame.setdefault("frame", []).append(int(sp[0]))
                self.by_frame.setdefault("nrefs", []).append(int(sp[1]))
                self.by_frame.setdefault("nmisfits", []).append(int(sp[2]))
                self.by_frame.setdefault("iobs", []).append(float(sp[3]))
                self.by_frame.setdefault("sigma", []).append(float(sp[4]))
                self.by_frame.setdefault("ios", []).append(float(sp[5]))
                self.by_frame.setdefault("peak", []).append(float(sp[6]))
                self.by_frame.setdefault("corr", []).append(float(sp[7]))
                self.by_frame.setdefault("rmeas", []).append(float(sp[8]))
                self.by_frame.setdefault("nrmeas", []).append(int(sp[9]))
                self.by_frame.setdefault("nuniq", []).append(int(sp[10]))
            elif "Framediff #refs R_d n-notfriedel Rd-notfriedel n-friedel Rd-friedel" in l:
                reading = "by_framediff"
            elif reading == "by_framediff":
                if "For inline graphs use a Java browser" in l:
                    reading = ""
                    continue

                if len(l) < 71: continue

                # Framediff #refs R_d n-notfriedel Rd-notfriedel n-friedel Rd-friedel, (dummy)
                sp = l[:6], l[6:14], l[14:24], l[24:32], l[32:42], l[42:50], l[50:60]

                self.by_framediff.setdefault("framediff", []).append(int(sp[0]))
                self.by_framediff.setdefault("nrefs", []).append(int(sp[1]))
                self.by_framediff.setdefault("rd", []).append(float(sp[2]))
    # parse()
# class XdsstatLp
