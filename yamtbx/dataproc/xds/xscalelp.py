"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
import re
import os
import collections
from yamtbx.dataproc.xds import correctlp
from yamtbx.dataproc import xds
from yamtbx.dataproc import cbf

re_data_info = re.compile("([0-9]+) *([-\.0-9E\+]+) * ([0-9]+) *([0-9]+) * ([^ ]*)")

def get_pairwise_correlations(lpin):
    read_flag = False
    ret = []
    for l in open(lpin):
        if "  #i   #j     REFLECTIONS     BETWEEN" in l:
            read_flag = True
        elif read_flag and l.strip() != "":
            sp = l.split()
            if len(sp) != 6 or re.search("[A-Za-z]", l):
                break

            i, j, common_refs, corr, ratio, bfac = sp
            ret.append((int(i), int(j), int(common_refs), float(corr), float(ratio), float(bfac)))
    return ret
# get_pairwise_correlations()


def construct_data_graph(lpin, min_common_refs=10):
    import networkx as nx

    G = nx.Graph()
    data = get_pairwise_correlations(lpin)
    for i, j, common_refs, corr, ratio, bfac in data:
        if common_refs > min_common_refs:
            G.add_edge(i-1, j-1)

    return G
# construct_data_graph

def read_no_common_ref_datasets(lpin):
    ret = []
    for l in open(lpin):
        if "no common reflections with data set" in l:
            idx = int(l.split()[-1])
            ret.append(idx-1)
    return ret
# read_no_common_ref_datasets()

def get_read_data(lpin):
    """
 DATA    MEAN       REFLECTIONS        INPUT FILE NAME
 SET# INTENSITY  ACCEPTED REJECTED
   1  0.1283E+04     1062      0  ../../XDS_ASCII.HKL
   """

    read_flag = False
    ret = []
    for l in open(lpin):
        if l.startswith(" SET# INTENSITY  ACCEPTED REJECTED"):
            read_flag = True
        elif read_flag:
            if "************" in l: break

            r = re_data_info.search(l)
            if r:
                idx, mean_i, n_accepted, n_rejected, filename = r.groups()
                ret.append([int(idx), float(mean_i), int(n_accepted), int(n_rejected), filename.strip()])
    return ret
# get_read_data()

def get_k_b(lpin):
    """
     K        B           DATA SET NAME
    1.0000   0.000    XDS_ASCII_fullres.HKL
    """

    read_flag = False
    ret = [] # list of [K, B, filename]
    for l in open(lpin):
        if "K        B           DATA SET NAME" in l:
            read_flag = True
        elif read_flag:
            if "********************" in l: break
            sp = l.strip().split()
            if len(sp) == 3:
                ret.append([float(sp[0]), float(sp[1]), sp[2]])
    return ret
# get_k_b()

def get_ISa(lpin):
    """
     a        b          ISa    ISa0   INPUT DATA SET
 1.986E+00  1.951E-01    1.61   50.00 /isilon/users/target/target/Iwata/_proc_ox2r/150415-hirata/1010/06/DS/multi011_1-5/XDS_ASCII_fullres.HKL
    """

    read_flag = False
    ret = [] # list of [a, b, ISa, ISa0, filename]
    for l in open(lpin):
        if "a        b          ISa    ISa0   INPUT DATA SET" in l:
            read_flag = True
        elif read_flag:
            if "********************" in l: break

            sp = l.strip().split()
            if len(sp) == 5:
                a, b, ISa, ISa0 = map(float, sp[:4])
                ret.append([a,b,ISa,ISa0, sp[4]])
    return ret
# get_ISa()

def get_rfactors_for_each(lpin):
    """
  R-FACTORS FOR INTENSITIES OF DATA SET /isilon/users/target/target/Iwata/_proc_ox2r/150415-hirata/1010/06/DS/multi011_1-5/XDS_ASCII_fullres.HKL

 RESOLUTION   R-FACTOR   R-FACTOR   COMPARED
   LIMIT      observed   expected

     5.84        60.4%      50.1%       174
     4.13        58.1%      51.5%       310
     3.38        60.0%      54.6%       410
     2.92        90.3%      76.1%       483
     2.62       130.4%     100.3%       523
     2.39       241.1%     180.5%       612
     2.21       353.9%     277.9%       634
     2.07       541.1%     444.0%       673
     1.95       -99.9%     -99.9%       535
    total        84.5%      71.2%      4354
    """

    read_flag = False
    filename = None
    ret = collections.OrderedDict() # {filename: list of [dmin, Robs, Rexpt, Compared]}
    for l in open(lpin):
        if "R-FACTORS FOR INTENSITIES OF DATA SET" in l:
            filename = l.strip().split()[-1]
        elif "LIMIT      observed   expected" in l:
            read_flag = True
        elif read_flag:
            sp = l.strip().replace("%","").split()
            if len(sp) == 4:
                dmin, robs, rexp, compared = sp
                if dmin != "total": dmin = float(dmin)
                else: dmin, read_flag = None, False

                robs, rexp = map(float, (robs, rexp))
                compared = int(compared)
                ret.setdefault(filename, []).append([dmin, robs, rexp, compared])
    return ret
# get_rfactors_for_each()

def snip_symm_and_cell(lpin):
    s = ""
    read_flag = False
    for l in open(lpin):
        if "THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:" in l:
            read_flag = True
        elif read_flag:
            s += l
            if "UNIT_CELL_CONSTANTS=" in l: break

    return s
# snip_symm_and_cell()

def snip_stats_table(lpin):
    s = ""
    read_flag = False
    for l in open(lpin):
        if "SUBSET OF INTENSITY DATA WITH SIGNAL/NOISE >= -3.0 AS FUNCTION OF RESOLUTION" in l:
            read_flag = True

        if read_flag:
            s += l

        if "total" in l:
            break
    return s
# snip_stats_table()

def read_stats_table(lpin):
    lines = snip_stats_table(lpin).splitlines()
    if len(lines)<2 or lines[1].strip() != "RESOLUTION     NUMBER OF REFLECTIONS    COMPLETENESS R-FACTOR  R-FACTOR COMPARED I/SIGMA   R-meas  CC(1/2)  Anomal  SigAno   Nano":
        return None

    table = {}

    for ll in lines[4:]:
        sp = correctlp.table_split(ll)
        assert len(sp) == 14
        table.setdefault("dmin", []).append(float(sp[0]) if sp[0]!="total" else None)
        table.setdefault("nuniq", []).append(int(sp[2]))
        table.setdefault("redundancy", []).append(float(sp[1])/float(sp[2]) if float(sp[2]) > 0 else 0)
        table.setdefault("cmpl", []).append(float(sp[4][:-1]))
        table.setdefault("r_merge", []).append(float(sp[5][:-1]))
        table.setdefault("i_over_sigma", []).append(float(sp[8]))
        table.setdefault("r_meas", []).append(float(sp[9][:-1]))
        table.setdefault("cc_half", []).append(float(sp[10].replace("*","")))
        table.setdefault("cc_ano", []).append(float(sp[11].replace("*","")))
        table.setdefault("sig_ano", []).append(float(sp[12]))

    for i in xrange(len(table["dmin"])):
        lrange = float("Inf") if i in (0,len(table["dmin"])-1) else table["dmin"][i-1]
        rrange = table["dmin"][i] if i != len(table["dmin"])-1 else table["dmin"][i-1]
        table.setdefault("d_range", []).append((lrange, rrange))

    return table
# read_stats_table()

def snip_control_cards(lpin):
    ret = ""
    read_flag = False
    for l in open(lpin):
        if "CONTROL CARDS" in l:
            read_flag = True
        elif read_flag:
            if "*****" in l or l.strip()=="": continue
            if "THE DATA COLLECTION STATISTICS REPORTED BELOW ASSUMES:" in l: break
            ret += l
    return ret
# snip_control_cards()

def read_control_cards(lpin):
    res = []
    # XXX Need special care for xscale specific manner (order matters!)

    for l in snip_control_cards(lpin).splitlines():
        if "!" in l:  l = l[:l.find("!")] # Remove comment
        r = xds.re_xds_kwd.findall(l)
        res.extend(r)

    return res
# read_control_cards()

def cbf_to_dat(lpin):
    lpdir = os.path.dirname(lpin)
    filename = None
    params = {}

    ofs_abs = open(os.path.join(lpdir, "corfac_absorp.dat"), "w")
    ofs_abs.write("file ix xmin xmax ipos posx posy fac\n")
    ofs_mod = open(os.path.join(lpdir, "corfac_modpix.dat"), "w")
    ofs_mod.write("file ix xmin xmax iy ymin ymax fac\n")
    ofs_dec = open(os.path.join(lpdir, "corfac_decay.dat"), "w")
    ofs_dec.write("file ix xmin xmax iy ymin ymax fac\n")

    for l in open(lpin):
        if l.startswith(" CORRECTION FACTORS for visual inspection"):
            filename = l.split()[-1]
        elif l.startswith(" NUMBER OF REFLECTIONS USED FOR DETERMINING"):
            nref = int(l.split()[-1]) # not used now

            if "_***" in filename:
                # reset
                filename = None
                params = {}
                continue

            # open cbf
            data, ndimfast, ndimmid = cbf.load_minicbf_as_numpy(os.path.join(lpdir, filename))
            data = data.reshape(ndimmid, ndimfast)

            filenum = filename[filename.index("_")+1:-4]
            xmin, xmax = float(params["XMIN"][0]), float(params["XMAX"][0])
            nxbin = int(params["NXBIN"][0])
            xstep = (xmax-xmin)/nxbin

            if "ABSORP_" in filename:
                for ix in xrange(nxbin):
                    x1, x2 = xmin + ix*xstep, xmin + (ix+1)*xstep
                    for ipos in xrange(len(params["DETECTOR_SURFACE_POSITION"])):
                        pos = params["DETECTOR_SURFACE_POSITION"][ipos]
                        fac = data[ipos, ix]
                        ofs_abs.write("%s %2d %.2f %.2f %2d %s %5d\n" % (filenum, ix, x1, x2, ipos, pos, fac))
            elif "DECAY_" in filename or "MODPIX_" in filename:
                ymin, ymax = float(params["YMIN"][0]), float(params["YMAX"][0])
                nybin = int(params["NYBIN"][0])
                ystep = (ymax-ymin)/nybin

                for ix in xrange(nxbin):
                    x1, x2 = xmin + ix*xstep, xmin + (ix+1)*xstep
                    for iy in xrange(nybin):
                        y1, y2 = ymin + iy*ystep, ymin + (iy+1)*ystep
                        fac = data[iy, ix]
                        ofs = ofs_mod if "MODPIX_" in filename else ofs_dec
                        ofs.write("%s %2d %.2f %.2f %2d %.5f %.5f %5d\n" % (filenum, ix, x1, x2, iy, y1, y2, fac))
            else:
                print "What is this file!?", filename

            # reset
            filename = None
            params = {}
        elif filename is not None:
            r = xds.re_xds_kwd.findall(l)
            for k, v in r:
                params.setdefault(k, []).append(v)

