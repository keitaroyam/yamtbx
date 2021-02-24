def parse_dat(f_in, ret):
    ifs = open(f_in)
    first = ifs.readline()
    if first.startswith("  1/d centre"):
        key = ""
        if "Rsplit/%" in first:
            key = "table_rsplit"
        elif first.split()[2] == "CC":
            key = "table_cc"
        else:
            return None

        if key in ret: print "Warning: duplidated %s" %key
        ret[key] = {}
        for l in ifs:
            _, fom, nref, d, dmax, dmin = l.split()
            dmax, dmin = 10./float(dmax), 10./float(dmin)
            ret[key].setdefault("dmax", []).append("%7.3f"%dmax)
            ret[key].setdefault("dmin", []).append("%7.3f"%dmin)
            if key == "table_rsplit":
                fom = "%.4f"%(float(fom)/100.)
            elif key == "table_cc":
                fom = "%.4f"%(float(fom))
            ret[key].setdefault("value", []).append(fom)

    elif first.startswith("Center 1/nm"):
        if "table" in ret: print "Warning: duplicated table"
        ret["table"] = {}
        nposs_all = 0
        nref_all = 0
        nmeas_all = 0
        for l in ifs:
            _, nref, nposs, cmpl, nmeas, red, snr, _, _, d, dmax, dmin = l.split()
            dmax, dmin = 10./float(dmax), 10./float(dmin)
            ret["table"].setdefault("dmax", []).append("%7.3f"%dmax)
            ret["table"].setdefault("dmin", []).append("%7.3f"%dmin)
            ret["table"].setdefault("snr", []).append(snr)
            ret["table"].setdefault("cmpl", []).append(cmpl)
            ret["table"].setdefault("redun", []).append(red)
            ret["table"].setdefault("nuniq", []).append(nref)
            ret["table"].setdefault("nposs", []).append(nposs)
            nposs_all += int(nposs)
            nref_all += int(nref)
            nmeas_all += int(nmeas)

        ret["cmpl"] = "%6.2f"%(100.*nref_all/nposs_all)
        ret["redun"] = "%6.1f"%(nmeas_all/nref_all)
# parse_dat()

def parse_file(logfile, ret):
    parse_dat(logfile, ret)

    for l in open(logfile):
        if "Overall Rsplit = " in l:
            if "rsplit" in ret: print "Warning: duplicated Rsplit"
            ret["rsplit"] = "%.3f"%(float(l[l.rindex("=")+1:l.rindex("%")])/100.)
        elif "Overall CC = " in l:
            if "cc12" in ret: print "Warning: duplicated CC"
            ret["cc12"] = "%.4f"%(float(l[l.rindex("=")+1:].strip()))
        elif "overall <snr> = " in l:
            if "snr" in ret: print "Warning: duplicated SNR"
            ret["snr"] = "%.2f"%(float(l[l.rindex("=")+1:].strip()))
        elif "measurements in total" in l:
            if "nmeas" in ret: print "Warning: duplicated measurements"
            ret["nmes"] = l.split()[0]
        elif "reflections in total" in l:
            if "nuniq" in ret: print "Warning: duplicated reflections"
            ret["nuniq"] = l.split()[0]
        elif "Accepted resolution range" in l:
            ret["dmin"] = l[l.rindex("to")+3:l.rindex("Ang")].strip()
            ret["dmax"] = l[l.rindex("(")+1:l.rindex(" to")].strip()
# parse_file()

def run(log_files):
    ret = {}
    for f in log_files:
        parse_file(f, ret)

    s = """\
_reflns.entry_id                         UNNAMED
_reflns.d_resolution_low                 %(dmax)s
_reflns.d_resolution_high                %(dmin)s
_reflns.number_all                       ?
_reflns.number_obs                       %(nuniq)s
_reflns.observed_criterion               ?
_reflns.observed_criterion_F_max         ?
_reflns.observed_criterion_F_min         ?
_reflns.observed_criterion_I_max         ?
_reflns.observed_criterion_I_min         ?
_reflns.observed_criterion_sigma_F       ?
_reflns.observed_criterion_sigma_I       ?
_reflns.percent_possible_obs             %(cmpl)s
_reflns.pdbx_redundancy                  %(redun)s
_reflns.pdbx_netI_over_sigmaI            %(snr)s
_reflns.pdbx_number_measured_all         ?
_reflns.pdbx_diffrn_id                   1
_reflns.pdbx_ordinal                     1
_reflns.pdbx_CC_half                     %(cc12)s
_reflns.pdbx_R_split                     %(rsplit)s
""" % ret

    s += """\
#
loop_
_reflns_shell.d_res_low
_reflns_shell.d_res_high
_reflns_shell.meanI_over_sigI_obs
_reflns_shell.number_unique_obs
_reflns_shell.percent_possible_all
_reflns_shell.pdbx_redundancy
_reflns_shell.pdbx_diffrn_id
_reflns_shell.pdbx_CC_half
_reflns_shell.pdbx_R_split
"""
    for i in xrange(len(ret["table"]["dmax"])):
        tmp = dict(dmin=ret["table"]["dmin"][i],
                   dmax=ret["table"]["dmax"][i],
                   snr=ret["table"]["snr"][i],
                   nuniq=ret["table"]["nuniq"][i],
                   cmpl=ret["table"]["cmpl"][i],
                   redun=ret["table"]["redun"][i],
                   cc12=ret["table_cc"]["value"][i],
                   rsplit=ret["table_rsplit"]["value"][i])
        s += "%(dmax)s %(dmin)s %(snr)s %(nuniq)s %(cmpl)s %(redun)s 1 %(cc12)s %(rsplit)s\n" % tmp

    print s

if __name__ == "__main__":
    import sys
    log_files = sys.argv[1:]
    run(log_files)
