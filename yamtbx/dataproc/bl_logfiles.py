"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals
import re
import os
import datetime
import math
import traceback
from collections import OrderedDict

from yamtbx.dataproc.dataset import template_to_filenames, re_pref_num_ext
from functools import reduce

class ScanInfo(object):
    def __init__(self):
        self.vpoints, self.vstep = 0, 0
        self.hpoints, self.hstep = 0, 0
        self.date = datetime.date.fromtimestamp(0) # 1970-1-1
        self.filename_template = None
        self.wavelength = float("nan")
        self.scan_direction, self.scan_path = None, None
        self.attenuator = ("??", float("nan"))
        self.beam_hsize, self.beam_vsize = -1., -1.
        self.osc_start, self.osc_step = float("nan"), float("nan")
        self.fixed_spindle = float("nan")
        self.frame_rate = float("nan")
        self.exp_time = -1.
        self.distance = -1.
        self.has_extra_images = True # whether one extra image in each row is created in shutterless mode.
        self.need_treatment_on_image_numbers = True # old bss does not write "dummy" lines in shutterless mode and therefore number in the logfile has offset from filename.
        self.filename_coords = [] # becomes dict later?
        self.filename_idxes = [] # for speed?
        self.filename_tags = [] # for SACLA
    # __init__ ()

    def is_shutterless(self):
        return self.fixed_spindle == self.fixed_spindle and self.frame_rate == self.frame_rate
    # is_shutterless()

    def match_file_with_template(self, filename):
        tmpl = os.path.basename(self.filename_template)
        re_tmpl = re.compile("^" + re.escape(tmpl).replace("\\?","[0-9]") + "$")
        return re_tmpl.search(os.path.basename(filename)) is not None
    # match_file_with_template()

    def get_file_number_based_on_template(self, filename):
        tmpl = os.path.basename(self.filename_template)
        re_tmpl = re.compile("^" + re.escape(tmpl).replace("\\?","[0-9]").replace("[0-9]","([0-9])").replace(")(","") + "$")
        return re_tmpl.search(os.path.basename(filename))
    # get_file_number_based_on_template()

    def get_prefix(self): return self.filename_template[:self.filename_template.index("?")]

    def scan_completed(self):
        return self.vpoints*self.hpoints == len(self.filename_coords)

# class ScanInfo

class BssDiffscanLog(object):
    def __init__(self, scanlog):
        self.scanlog = scanlog
        self.parse()
        #for scan in self.scans:
        #    print scan.date, scan.filename_template
    # __init__()

    def parse(self):
        re_scanstart = re.compile("Diffraction scan \(([0-9A-Za-z/:\[\] ]+)\)")
        re_nums = re.compile("([0-9]+) +([\-0-9\.]+|dummy) +([\-0-9\.]+|dummy) +([\-0-9\.]+|dummy)")
        re_point_step = re.compile("point: *([0-9]+) *step: *([0-9\.]+)")
        re_osc = re.compile("Oscillation start: ([-0-9\.]+) \[deg\], step: ([-0-9\.]+) \[deg\]")
        re_att = re.compile("Attenuator: +([^ ]+) +([0-9]+)um")
        re_att2 = re.compile("Attenuator transmission: +([^ ]+) +\(([^ ]+) attenuator: ([0-9]+)\[um\]\)")
        re_exp = re.compile("Exp\. time: ([\.0-9]+) ")
        re_beam = re.compile("hor\. beam size: +([\.0-9]+)\[um\], ver\. beamsize: +([\.0-9]+)\[um\]") # old bss, wrong (need to swap h/v)
        re_beam2 = re.compile("horizontal size: +([\.0-9]+)\[um\], vertical size: +([\.0-9]+)\[um\]") # new bss (2015-Apr), correct
        re_fixed_spindle = re.compile("Fixed spindle angle: ([-0-9\.]+) \[deg\]")
        re_frame_rate = re.compile("Frame rate: ([-0-9\.]+) \[frame/s\]")

        self.scans = []
        self.filename_gonio_gc = OrderedDict()

        for l in open(self.scanlog):
            r_scanstart = re_scanstart.search(l)
            if r_scanstart:
                datestr = re.sub("\[[A-Za-z]+\] ", "", r_scanstart.group(1))
                if datestr.count(":") == 3: datestr = datestr[:datestr.rindex(":")] # from bss jan29-2015, millisec time is recorded; but we discard it here.
                date = datetime.datetime.strptime(datestr, "%Y/%m/%d %H:%M:%S")
                self.scans.append(ScanInfo())
                self.scans[-1].date = date
                continue

            r_osc = re_osc.search(l)
            if r_osc:
                self.scans[-1].osc_start = float(r_osc.group(1))
                self.scans[-1].osc_step = float(r_osc.group(2))
                continue

            r_exp = re_exp.search(l)
            if r_exp:
                self.scans[-1].exp_time = float(r_exp.group(1))
                continue

            r_beam = re_beam.search(l)
            if r_beam:
                self.scans[-1].beam_hsize = float(r_beam.group(2))
                self.scans[-1].beam_vsize = float(r_beam.group(1))
                continue

            r_beam2 = re_beam2.search(l)
            if r_beam2:
                self.scans[-1].beam_hsize = float(r_beam2.group(1))
                self.scans[-1].beam_vsize = float(r_beam2.group(2))
                continue

            r_fspindle = re_fixed_spindle.search(l)
            if r_fspindle:
                self.scans[-1].fixed_spindle = float(r_fspindle.group(1))
                continue

            r_frate = re_frame_rate.search(l)
            if r_frate:
                self.scans[-1].frame_rate = float(r_frate.group(1))
                continue

            if "No dummy image generated in each scan" in l:
                # this scan does not have dummy images.
                self.scans[-1].has_extra_images = False
                self.scans[-1].need_treatment_on_image_numbers = False
                continue

            if "Dummy images generated in each scan!" in l:
                # this scan actually has dummy images but diffscan.log has "dummy" lines
                self.scans[-1].has_extra_images = True
                self.scans[-1].need_treatment_on_image_numbers = False
                continue

            if "Scan direction:" in l:
                self.scans[-1].scan_direction = l[l.index(":")+1:].strip()
                continue                

            if "Scan path:" in l:
                self.scans[-1].scan_path = l[l.index(":")+1:].strip()
                continue

            if "Wavelength: " in l:
                self.scans[-1].wavelength = float(l.strip().split()[1])
                continue

            if "Attenuator" in l:
                r_att = re_att.search(l)
                r_att2 = re_att2.search(l)
                if r_att:
                    self.scans[-1].attenuator = (r_att.group(1), int(r_att.group(2)))
                elif r_att2:
                    self.scans[-1].attenuator = (r_att2.group(2), int(r_att2.group(3))) # (1) is transmisttance
                continue

            if "Cameralength: " in l:
                self.scans[-1].distance = float(l.strip().split()[1])
                continue

            if "FILE_NAME = " in l:
                self.scans[-1].filename_template =  os.path.basename(l[l.index("FILE_NAME = ")+12:].strip())

                # VERY DIRTY FIX!
                # Explanation: if diffscan.log contains ".h5", we assume it's Eiger hdf5 and hit-finding is done by streaming mode.
                if self.scans[-1].filename_template.endswith(".h5"):
                    self.scans[-1].filename_template = self.scans[-1].filename_template[:-3] + ".img"
                continue

            if "Vertical   scan" in l:
                vpoint, vstep = re_point_step.search(l.strip()).groups()
                self.scans[-1].vpoints, self.scans[-1].vstep = int(vpoint), float(vstep)
                continue

            if "Horizontal scan" in l:
                hpoint, hstep = re_point_step.search(l.strip()).groups()
                self.scans[-1].hpoints, self.scans[-1].hstep = int(hpoint), float(hstep)
                continue

            r = re_nums.search(l)
            if r:
                if "dummy" in r.groups():
                    assert not self.scans[-1].need_treatment_on_image_numbers
                    continue
                else:
                    gonio = tuple([float(x) for x in r.groups()[1:]])

                num = int(r.group(1))
                if self.scans[-1].need_treatment_on_image_numbers and self.scans[-1].is_shutterless() and self.scans[-1].hpoints > 1:
                    # This if section is not needed *if* BSS no longer creates such non-sense files.
                    # this should be an option.
                    num += int(math.ceil(float(num)/self.scans[-1].hpoints)) - 1

                filename = template_to_filenames(self.scans[-1].filename_template, num, num)[0]

                grid_coord = self.get_grid_coord_internal(self.scans[-1].vpoints, self.scans[-1].vstep,
                                                          self.scans[-1].hpoints, self.scans[-1].hstep,
                                                          num, self.scans[-1].is_shutterless() and self.scans[-1].has_extra_images,
                                                          self.scans[-1].scan_direction, self.scans[-1].scan_path)

                # XXX Need to check if the filename is already registered (file is overwritten!!)
                self.scans[-1].filename_coords.append((filename, (gonio, grid_coord)))
                self.scans[-1].filename_idxes.append((filename, num))
                if len(l.split())==5: # for SACLA
                    self.scans[-1].filename_tags.append((filename, int(l.split()[-1])))
                self.filename_gonio_gc[os.path.basename(filename)] = (gonio, grid_coord)
    # parse()

    def __getitem__(self, filename):
        for scan in reversed(self.scans):
            print(scan.filename_coords)
            if filename in scan.filename_coords:
                return scan
        raise KeyError
    # __getitem__()

    def get_grid_coord(self, filename):
        h = self.filename_gonio_gc.get(filename, None)
        if h is not None:
            return h[1]
    # get_grid_coord()

    def get_gonio_xyz(self, filename):
        h = self.filename_gonio_gc.get(filename, None)
        if h is not None:
            return h[0]
    # get_gonio_xyz()

    @staticmethod
    def get_grid_coord_internal(vpoint, vstep, hpoint, hstep, num, include_extra_files, scan_direction, scan_path):
        voffset = 0 if vpoint % 2 == 1 else -0.5
        hoffset = 0 if hpoint % 2 == 1 else -0.5
        
        if vpoint == 1:
            # right -> left
            return (((hpoint//2 + 1) - num)*hstep, 0)
        if hpoint == 1:
            # top -> bottom
            return (0, ((vpoint//2 + 1)-num)*vstep)
        # XXX In case of hpoint==1 and 2D scan, dummy images are created for each row.
        # XXX However, we cannot know whether it is 2D scan or 1D scan from diffscan.log.
        # XXX Of course we can know it after scan (by counting dummy images).. But, in the first place, this should not be done!


        # right->left, top->bottom
        
        if include_extra_files and hpoint > 1:
            # This if section is not needed *if* BSS no longer creates such non-sense files.
            # this should be an option.
            if scan_direction in ("horizontal", None):
                v = (num-1)//(hpoint+1)
                h = (num-1)%(hpoint+1)
                if h==hpoint: h = float("nan")
            else: # vertical
                v = (num-1)%(vpoint+1)
                h = (num-1)//(vpoint+1)
                if v==vpoint: v = float("nan")
        else:
            if scan_direction in ("horizontal", None):
                v = (num-1)//hpoint
                h = (num-1)%hpoint
            else: # vertical
                v = (num-1)%vpoint
                h = (num-1)//vpoint

        x, y = (hpoint//2 - h + hoffset)*hstep, (vpoint//2 - v + voffset)*vstep

        # scan_path can be zig-zag or normal (in case of streaming) or None
        if scan_path == "zig-zag":
            if scan_direction in ("horizontal", None):
                if int(v%2) == 0: return x, y
                else: return -x, y
            else:
                if int(h%2) == 0: return x, y
                else: return x, -y
        else:
            return x, y

    # get_grid_coord()

    def calc_grid_coord(self, prefix=None, num=None, filename=None):
        if filename is not None:
            assert (prefix, num).count(None)==2
            bf = os.path.basename(filename)
            prefix = bf[:bf.rindex("_")+1]
            num = int(bf[bf.rindex("_")+1:bf.rindex(".")])
        else:
            assert None not in (prefix, num)
            
        matched = [x for x in self.scans if x.get_prefix()==prefix]
        if not matched: return None
        scan = matched[-1]
        return self.get_grid_coord_internal(scan.vpoints, scan.vstep,
                                            scan.hpoints, scan.hstep,
                                            num, scan.is_shutterless() and scan.has_extra_images,
                                            scan.scan_direction, scan.scan_path)
    # calc_grid_coord()

    def remove_overwritten_scans(self):
        table = {}
        for i, scan in enumerate(self.scans):
            if scan.filename_coords:
                table.setdefault(scan.filename_template, []).append(i)

        table = [x for x in list(table.items()) if len(x[1])>1]
        rem_idxes = [x[1][:-1] for x in table] # last indexes will be alive
        if not rem_idxes: return
        rem_idxes = reduce(lambda x,y:x+y, rem_idxes) 
        for i in sorted(rem_idxes, reverse=True):
            del self.scans[i]
    # remove_overwritten_scan()
# class BssDiffscanLog


def parse_bss_diffscan_log(scanlog):
    ret = []
    bdl = BssDiffscanLog(scanlog)

    for scan in bdl.scans:
        for f, coords in scan.filename_coords:
            ret.append((f, coords))

    return ret
# parse_bss_diffscan_log()

class ImageStatus(object):
    def __init__(self, line=None):
        items_str = ("imagenum", "wavelength", "exp_time", "osc_start", "osc_step", "gonio_xyz", "filename", "ring_current", "I0", "time", "date")
        for s in items_str:
            setattr(self, s, None)
        self.overwritten = False

        if line is not None:
            self.parse_line(line)
    # __init__()

    def parse_line(self, line):
        sp = [s.strip() for s in line.split(",")]
        if len(sp) == 14:
            junk, self.imagenum, self.wavelength, self.exp_time, self.osc_start, self.osc_step, gx,gy,gz, self.filename, self.ringcurrent, self.I0, self.time, self.date = sp
            self.gonio_xyz = (gx, gy, gz)
        elif len(sp) == 11:
            junk, self.imagenum, self.wavelength, self.exp_time, self.osc_start, self.osc_step, self.filename, self.ringcurrent, self.I0, self.time, self.date = sp

        self.datetime = datetime.datetime.strptime("%s %s"%(self.date, self.time), "%Y/%m/%d %H:%M:%S")
    # parse_line()
# class ImageStatus


re_sample = re.compile("Tray:(.*), Well:([0-9]+)")
re_scan = re.compile("scan_from = ([-0-9\.]+) *\[deg\], *scan_to = ([-0-9\.]+) *\[deg\], *scan_step = ([-0-9\.]+) *\[deg\]")
re_advcen_mode = re.compile("mode = ([^ ]*), type = ([^ ]*)")
re_advcen_pts = re.compile("adv_npoint = ([0-9]*), adv_step = ([0-9\.]*)\[mm\], adv_interval = ([0-9]*)")
re_beam_size = re.compile("Beam shape = ([^ ]*), Beam size = ([0-9\.]*)\[um\] x ([0-9\.]*)\[um\] \(h x w\)")
re_sampling = re.compile("sampling_interval = ([0-9]*), number_of_images = ([0-9]*) x ([0-9]*)")
re_delay = re.compile("delay_time = ([0-9\.]*)\[msec\], cameralength = ([0-9\.]*)\[mm\], attenuator = (.*)")
re_wavelen = re.compile("wavelength = ([0-9\.]*), expose_time = ([0-9\.]*)\[sec\]")
re_gonio_center = re.compile("center #([0-9]*): +([-0-9\.]*) +([-0-9\.]*) +([-0-9\.]*)")

def interpret_attenuator_label(label):
    if label == "None":
        return ("None", 0)
    elif "T=" in label:
        r = re.search("([-A-Za-z0-9\. ]+)\[([a-z]+)\] \(T= *([-0-9\.]+)\)", label)
        if r:
            unit = r.group(2)
            fac = 1.
            if unit == "mm": fac = 1.e3
            mat, thick = "??", float("nan")
            sp = r.group(1).split()
            if len(sp) == 2:
                mat, thick = sp[0], float(sp[1])*fac
            else:
                thick = float(sp[0])*fac
                
            return (mat, thick)
        else:
            return ("??", float("nan"))
    else:
        # old style
        sp = label.split()
        if len(sp) == 2 and sp[1].endswith(("mm","um")):
            fac = 1e3 if sp[1].endswith("mm") else 1
            thickness = float(sp[1][:-2]) * fac
            return (sp[0], thickness)
        else:
            return ("??", float("nan"))
# interpret_attenuator_label()

class JobInfo(object):
    """
 JOB_ID#   = 568
 Beamline  = SACLABL3
 Sample  = Tray:1518, Well:4
 FILE_NAME = /data/psii/20131122/27S/27S052_d0000/27S052_d0000_??????.img
 JOB_MODE  = Crystal Check
 DETECTOR  = CCD (MX225HE) :
  2theta = 0.0, binning = 2, fastadc = 1, dezinger = 0, transform = 1
 DATA COLLECTION_PARAMETERS :
  scan_from = 0.00[deg], scan_to = 12.00[deg], scan_step = 0.20[deg]
  sampling_interval = 1, number_of_images = 60 x 1
  delay_time = 100.0[msec], cameralength = 150.0[mm], attenuator = None
  ver. offset = 0.0[mm], hor. offset = 0.0[mm]
  wavelength = 1.23860, expose_time = 1.0[sec]

 Advanced Centering Parameters
  mode = vector_centering, type = auto_step
  adv_npoint = 60, adv_step = 0.0304[mm], adv_interval = 1
  center #1:  0.0538 -3.1154 -0.5815
  center #2:  0.0642 -4.8906 -0.8235

               DataSet, ImageNum, Wavelength, ExposeTime, StartAngle, OscAngle, Gonio X, Gonio Y, Gonio Z, FileName, RingCurrent, I0, Time ,Date
 IMAGE_STATUS000001 = 1,    1, 1.23860,  1.0,   0.0,   0.2, 0.054, -3.115, -0.582, 27S052_d0000_000001.img, 100.00, -10.00, 07:55:18, 2013/11/24
"""
    def __init__(self, logfilename):
        self.logfilename = logfilename

        items_str = ("job_id", "beamline", "sample", "filename", "prefix", "job_mode",
                     "beam_shape", "beam_size",
                     "detector", "twotheta", "binning", "fastadc", "dezinger", "transform",
                     "osc_start", "osc_end", "osc_step",
                     "osc_sample_interval", "n_images", "n_wavelengths", "delay_time", "distance", "attenuator",
                     "ver_offset", "hor_offset", "wavelength", "exp_time",
                     )
        for s in items_str:
            setattr(self, s, None)

        self.advanced_centering = {}
        self.images, self.scans = [], []
        self.beam_size = float("nan"), float("nan")
        self.status = None # Programmer can define
    # __init__()

    def parse_line(self, line):
        if line.startswith(" JOB_ID#   ="):
            self.job_id = line[line.index("=")+1:].strip()
        if line.startswith(" Sample  ="):
            r = re_sample.search(line)
            if r: self.sample = (r.group(1), int(r.group(2)))
        elif line.startswith(" Beamline  ="):
            self.beamline = line[line.index("=")+1:].strip()
        elif line.startswith(" FILE_NAME ="):
            self.filename = line[line.index("=")+1:].strip()
            tmp = os.path.basename(self.filename)
            if "?" in tmp: # no ? included when XAFS mode
                tmp = tmp[:tmp.index("?")]
                self.prefix = tmp[:-1] if tmp.endswith("_") else tmp
        elif "JOB_MODE  = " in line:
            # can be "Crystal Check", "XAFS", "Single wavelength", "Multi wavelength", "Unknown"
            self.job_mode = line[line.index("=")+1:].strip()
        elif "DETECTOR  =" in line:
            tmp = line[line.index("=")+1:].strip()
            if ":" in tmp: tmp = tmp[:tmp.index(":")].strip()
            self.detector = tmp
        elif "Beam shape =" in line:
            r = re_beam_size.search(line)
            if r:
                self.beam_shape = r.group(1)
                self.beam_size = (float(r.group(2)), float(r.group(3)))
        elif "scan_from" in line:
            r = re_scan.search(line)
            if r:
                self.osc_start, self.osc_end, self.osc_step = list(map(float, r.groups()))
        elif "sampling_interval = " in line:
            r = re_sampling.search(line)
            if r:
                self.osc_sample_interval = int(r.group(1))
                self.n_images = int(r.group(2))
                self.n_wavelengths = int(r.group(3))
        elif " delay_time =" in line:
            r = re_delay.search(line.strip())
            if r:
                self.delay_time = float(r.group(1)) # msec
                self.distance = float(r.group(2)) # mm
                self.attenuator = interpret_attenuator_label(r.group(3))
        elif "wavelength =" in line: # TODO support wavelength1 = ,...
            r = re_wavelen.search(line)
            if r:
                self.wavelength = float(r.group(1))
                self.exp_time = float(r.group(2))
        elif " mode =" in line:
            r = re_advcen_mode.search(line.strip())
            if r:
                self.advanced_centering["mode"] = r.group(1)
                self.advanced_centering["type"] = r.group(2)
            else:
                self.advanced_centering["mode"] = line[line.index("=")+1:].strip()
                
        elif " adv_npoint = " in line:
            r = re_advcen_pts.search(line.strip())
            if r:
                self.advanced_centering["adv_npoint"] = int(r.group(1))
                self.advanced_centering["adv_step"] = float(r.group(2))*1e3 # to um unit
                self.advanced_centering["adv_interval"] = int(r.group(3))
        elif "center #" in line: # NEED to also take care of "skipped point:"
            r = re_gonio_center.search(line)
            if r:
                self.advanced_centering.setdefault("centers", []).append(list(map(float, r.groups()[1:])))
        elif line.startswith(" IMAGE_STATUS"):
            self.images.append(ImageStatus(line))
        elif line.startswith(" SCAN_STATUS"):
            self.scans.append(ImageStatus(line))

    # parse_line()

    def get_frame_num_range(self):
        if len(self.images) > 0 :
            r_st, r_en = re_pref_num_ext.search(self.images[0].filename), re_pref_num_ext.search(self.images[-1].filename)
            if None not in (r_st, r_en):
                return (int(r_st.group(2)), int(r_en.group(2)))

        return (None, None)

    def get_frame_num_ranges_for_h5(self):
        """
        Returns frame number ranges assuming h5 file (multiple frames in single file)
        Not sure non-Eiger h5 file will follow this rule (Currently only Eiger writes h5 files)
        """
        if self.advanced_centering.get("mode", "") != "multiple_crystals":
            return [(1, self.n_images)]
        
        ret = []
        for i in range(len(self.advanced_centering.get("centers", []))):
            ret.append((i*self.n_images+1, (i+1)*self.n_images))
        return ret
    # get_frame_num_ranges_for_h5()

    def get_master_h5_if_exists(self):
        if self.prefix is None: return None
        masterh5 = os.path.join(os.path.dirname(self.filename), self.prefix+"_master.h5")
        if os.path.isfile(masterh5): return masterh5
        return None
    # get_master_h5_if_exists()

    def all_image_files_exist(self, nr=None, debug=True):
        from yamtbx.dataproc import eiger
        import h5py
        master_h5 = self.get_master_h5_if_exists()
        if master_h5: # XXX check only for nr?
            if nr:
                frames_requested = set(range(nr[0], nr[1]+1))
                frames_available = set(eiger.get_available_frame_numbers(master_h5))
                return frames_requested.issubset(frames_available)
            else:
                if debug: print("Checking if related files exist for %s" % master_h5)
                try:
                    files = eiger.get_masterh5_related_filenames(master_h5)
                except:
                    if debug:
                        print("Error when reading master h5 (%s)" % master_h5)
                        print(traceback.format_exc())
                    return False

                flag_ng = False
                for f in files:
                    if os.path.isfile(f):
                        try: h5py.File(f, "r")
                        except:
                            if debug: print(" file incomplete or broken: %s" % f)
                            flag_ng = True
                    else:
                        if debug: print(" not exists: %s" % f)
                        flag_ng = True

                return not flag_ng
        else:
            if nr:
                filenames = template_to_filenames(self.filename, nr[0], nr[1])
                return all([os.path.exists(x) for x in filenames])
            else:
                return True # XXX 

# class JobInfo

class BssJobLog(object):
    def __init__(self, joblog=None, remove_overwritten=False):
        self.jobs = []
        if joblog is not None:
            self.parse(joblog)


        self.annotate_overwritten_images(remove=remove_overwritten)
    # __init__()

    def parse(self, joblog):

        for l in open(joblog):
            if l.startswith(" JOB_ID#   ="):
                self.jobs.append(JobInfo(joblog))
            if len(self.jobs) > 0:
                self.jobs[-1].parse_line(l)
    # parse()

    def annotate_overwritten_images(self, remove=False):
        del_indices = {} # {i: [j,...]}
        # Assume all files in the same directory and the common prefix is in single file. Is it true??
        # Assume files in the same job are never overwritten. Only check inter-job files.

        filenames = [] # (i, j, filename)
        for i, job in enumerate(self.jobs):
            this_job_filenames = []
            for j, img in enumerate(job.images):
                fltr = [x for x in filenames if x[2]==img.filename]
                #assert len(fltr) in (0,1)
                if len(fltr) > 0:
                    #del_indices.append((fltr[-1][0],fltr[-1][1]))
                    del_indices.setdefault(fltr[-1][0], []).append(fltr[-1][1])
                    self.jobs[fltr[-1][0]].images[fltr[-1][1]].overwritten = True

                this_job_filenames.append((i,j,img.filename))
            filenames.extend(this_job_filenames)

        ow_flag = False
        if remove:
            for i in del_indices:
                for j in sorted(del_indices[i], reverse=True):
                    print("overwritten:", self.jobs[i].logfilename, self.jobs[i].job_id, self.jobs[i].images[j].filename)
                    del self.jobs[i].images[j]
                    ow_flag = True

        if ow_flag:
            print()
    # remove_overwritten_images()
# class BssJobLog
