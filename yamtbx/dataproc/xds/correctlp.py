"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from cctbx import sgtbx
from cctbx import uctbx
from yamtbx.util import safe_float

def get_ISa(lp, check_valid=False):
    read_flag = False
    for l in open(lp):
        if "a        b          ISa" in l:
            read_flag = True
        elif read_flag:
            tmp = l.split()
            if len(tmp) != 3:
                return float("nan")
            a, b, ISa = tmp
            if not check_valid or not (a == "4.000E+00" and b == "1.000E-04"):
                return float(ISa)
    return float("nan")
# get_ISa()

def get_P1_cell(lp, force_obtuse_angle=False):
    read_flag = False
    for l in open(lp):
        if l.startswith(" CHARACTER  LATTICE     OF FIT      a      b      c   alpha  beta gamma"):
            read_flag = True
        elif read_flag and l[14:16] == "aP":
            cell = map(float, (l[32:39], l[39:46], l[46:53], l[53:59], l[59:65], l[65:71]))
            if force_obtuse_angle:
                tmp = map(lambda x: (x[0]+3,abs(90.-x[1])), enumerate(cell[3:])) # Index and difference from 90 deg
                tmp.sort(key=lambda x: x[1], reverse=True)
                if cell[tmp[0][0]] < 90:
                    tmp = map(lambda x: (x[0]+3,90.-x[1]), enumerate(cell[3:])) # Index and 90-val.
                    tmp.sort(key=lambda x: x[1], reverse=True)
                    for i,v in tmp[:2]: cell[i] = 180.-cell[i]
            return uctbx.unit_cell(cell)
    return None
# get_P1_cell()

def table_split(l):
    # TODO need to change behavior by version?
    assert l[50] == l[61] == l[71] == l[98] == "%"
    return map(lambda x:x.strip(),
               (l[:9], # RESOLUTION LIMIT
                l[9:21], l[21:29], l[29:39], # NUMBER OF REFLECTIONS: OBSERVED  UNIQUE  POSSIBLE
                l[39:50], # COMPLETENESS
                l[51:61], l[62:71], # R-FACTOR: observed  expected
                l[72:81], # COMPARED
                l[81:89], # I/SIGMA
                l[89:98], # R-meas 
                l[99:107], # CC(1/2)
                l[108:114], # AnomalCorr
                l[115:123], # SigAno 
                l[123:] # Nano
               ))

# table_split

def errortable_split(l):
    # TODO need to change behavior by version?
    return map(lambda x:x.strip(),
               (l[:9], l[9:17], l[17:26], l[26:33], l[33:43], l[43:53], l[53:61], l[61:69], l[69:77]))
# errortable_split


class CorrectLp:
    def __init__(self, lpin):
        self.anomalous_flag = None

        if lpin is not None:
            self.parse(lpin)
    # __init__()

    def get_ISa(self): return self.a_b_ISa[2]
    def space_group_str(self): return "?" if self.space_group is None else self.space_group.info()
    def parse(self, lpin):
        reading = ""
        self.table = {}
        self.error_table = {}
        self.space_group, self.unit_cell, self.unit_cell_esd = None, None, None
        self.a_b_ISa = [float("nan"), float("nan"),float("nan")]
        self.snippets = {}

        lines = open(lpin).readlines()
        for i, l in enumerate(lines):
            if l.startswith(" FRIEDEL'S_LAW="):
                self.anomalous_flag = "FALSE" in l

            # ************ SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET ************
            elif "SELECTED SPACE GROUP AND UNIT CELL FOR THIS DATA SET" in l:
                reading = "selected_symmetry"
            elif reading == "selected_symmetry":
                # Info before refinement
                if "SPACE_GROUP_NUMBER=" in l:
                    self.space_group = sgtbx.space_group_info(int(l[l.index("=")+1:])).group()
                elif "UNIT_CELL_CONSTANTS=" in l:
                    # Will be overridden if refined
                    tmp = l[l.index("=")+1:]
                    self.unit_cell = map(float, (tmp[0:9],tmp[9:18],tmp[18:27],tmp[27:35],tmp[35:43],tmp[43:51]))
                    
                if "********" in l:
                    reading = ""

            # ******************************************************************************
            #  REFINEMENT OF DIFFRACTION PARAMETERS USING ALL IMAGES
            # ******************************************************************************
            elif l.startswith(" UNIT CELL PARAMETERS"):
                cell_params = map(float, (l[21:32], l[32:42], l[42:52], l[52:60], l[60:68], l[68:76]))
                self.unit_cell = cell_params
            elif l.startswith(" E.S.D. OF CELL PARAMETERS"):
                esd_params = map(float, (l[26:35], l[35:43], l[43:51], l[51:59], l[59:67], l[67:75]))
                self.unit_cell_esd = esd_params
            elif l.startswith(" SPACE GROUP NUMBER"):
                self.space_group = sgtbx.space_group_info(int(l[20:])).group()

            # ******************************************************************************
            #    CORRECTION PARAMETERS FOR THE STANDARD ERROR OF REFLECTION INTENSITIES
            # ******************************************************************************
            elif "a        b          ISa" in l:
                reading = "ISa"
                self.snippets["ISa"] = l
            elif reading == "ISa":
                a, b, ISa = map(float, l.split())
                reading = ""
                self.a_b_ISa = a, b, ISa
                self.snippets["ISa"] += l

            # ******************************************************************************
            #      STANDARD ERROR OF REFLECTION INTENSITIES AS FUNCTION OF RESOLUTION
            #      FOR DATA SET  XDS_ASCII.HKL
            # ******************************************************************************
            elif "RESOLUTION RANGE  I/Sigma  Chi^2  R-FACTOR  R-FACTOR  NUMBER ACCEPTED REJECTED" in l:
                reading = "error_table"
            elif reading == "error_table":
                if l.startswith(" ***"):
                    reading = ""
                    continue
                if len(l.split()) < 3: continue
                sp = errortable_split(l)
                # Note that the last line is about overall data
                self.error_table.setdefault("dmax", []).append(float(sp[0]))
                self.error_table.setdefault("dmin", []).append(float(sp[1]))
                self.error_table.setdefault("ios", []).append(safe_float(sp[2]))
                self.error_table.setdefault("chisq", []).append(safe_float(sp[3]))
                self.error_table.setdefault("r_merge", []).append(safe_float(sp[4]))
                self.error_table.setdefault("number", []).append(int(sp[6]))
                self.error_table.setdefault("nacc", []).append(int(sp[7]))
                self.error_table.setdefault("nrej", []).append(int(sp[8]))

            #if "SUMMARY OF DATA SET STATISTICS FOR IMAGE":

            # ******************************************************************************
            #  STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL" (DATA_RANGE=       x      x)
            # FILE TYPE:         XDS_ASCII      MERGE=FALSE          FRIEDEL'S_LAW=x
            # ******************************************************************************
            elif 'STATISTICS OF SAVED DATA SET "XDS_ASCII.HKL"' in l:
                reading = "stats_all"
                continue
            elif "SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION OF RESOLUTION" in l:
                self.snippets["table1"] = l
                if reading == "stats_all":
                    key = "all"
                    self.table[key] = {}
                    for ll in lines[i+1:i+4]: self.snippets["table1"] += ll
                    for ll in lines[i+4:i+4+10]:
                        self.snippets["table1"] += ll
                        sp = table_split(ll)
                        assert len(sp) == 14
                        self.table[key].setdefault("dmin", []).append(float(sp[0]) if sp[0]!="total" else None)
                        self.table[key].setdefault("redundancy", []).append(float(sp[1])/float(sp[2]) if float(sp[2]) > 0 else 0)
                        self.table[key].setdefault("cmpl", []).append(float(sp[4][:-1]))
                        self.table[key].setdefault("r_merge", []).append(float(sp[5][:-1]))
                        self.table[key].setdefault("i_over_sigma", []).append(float(sp[8]))
                        self.table[key].setdefault("r_meas", []).append(safe_float(sp[9][:-1]))
                        self.table[key].setdefault("cc_half", []).append(float(sp[10].replace("*","")))
                        self.table[key].setdefault("cc_ano", []).append(float(sp[11].replace("*","")))
                        self.table[key].setdefault("sig_ano", []).append(float(sp[12]))
                        if sp[0]=="total": # in case less than 9 shells
                            break
                    #print "".join(lines[i+4:i+4+10])
                    
    # parse()

    def resolution_based_on_ios_of_error_table(self, min_ios):
        # If all ios is less than min_ios, can't decide resolution.
        if "ios" not in self.error_table or max(self.error_table["ios"]) < min_ios:
            return float("nan")

        flag_last = False
        for dmin, ios in zip(self.error_table["dmin"], self.error_table["ios"]):
            if flag_last:
                return dmin
            if ios <= min_ios:
                flag_last = True
                continue

        return self.error_table["dmin"][-1] # the last
    # resolution_based_on_ios_of_error_table()

    def is_ISa_valid(self):
        if "ISa" not in self.snippets: return False

        sp = self.snippets["ISa"].splitlines()[-1].split()
        a, b, ISa = sp
        return not (a == "4.000E+00" and b == "1.000E-04") # in this case error model adjustment is failed.
    # is_ISa_valid()

if __name__ == "__main__":
    import sys
    lp = CorrectLp(sys.argv[1])
    print lp.table
    print lp.space_group.info(), lp.unit_cell
