"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
import os, glob, re
from yamtbx.dataproc import XIO
import traceback

IMG_EXTENSIONS = ".img", ".osc", ".cbf", ".mccd", ".mar1600", "_master.h5"
COMPRESS_EXTENSIONS = ".bz2", ".gz"

re_pref_num_ext = re.compile("(.*[^0-9])([0-9]+)\.(.*)")
re_pref_num = re.compile("(.*\.)([0-9]+)(\.)?((?(3).+))$") # .0001 style

def find_img_files(parentdir, recursive=True, skip_symlinks=False):
    """
    find files ending with IMG_EXTENSIONS and including digits
    @return image file names
    """
    matches = []
    re_including_digits = re.compile("[0-9]")
    re_noext = re.compile(".*\.[0-9]+(?:\.gz|\.bz2)?$")
    possible_extensions = tuple([ i+c for i in IMG_EXTENSIONS for c in COMPRESS_EXTENSIONS+("",) ])

    for root, dirnames, filenames in os.walk(parentdir, followlinks=not skip_symlinks):
        if not recursive and root != parentdir:
            continue

        for filename in filenames:
            if (filename.endswith(possible_extensions) and re_including_digits.search(filename)) or re_noext.search(filename):
                if skip_symlinks and os.path.islink(os.path.join(root, filename)):
                    continue
                matches.append(os.path.join(root, filename))

    return matches
# find_img_files()

def group_img_files_template(img_files, skip_0=False):
    """
    generate template string by image files
    This grouping method is only by "names" and doesn't read the header.

    @return [(template_string, min_frame, max_frame), ... ]
             e.g. [(A11_1A_??????.img', 1, 360), ...]
    """

    first_selection = {}

    # At first, grouping by pref+numbers+ext style
    for f in img_files:
        r = re_pref_num.search(f)
        if r: 
            first_selection.setdefault((r.group(1), r.group(4)), []).append(r.group(2))
            continue
        r = re_pref_num_ext.search(f)
        if r:
            first_selection.setdefault((r.group(1), r.group(3)), []).append(r.group(2))

    group = []

    for pref, ext in first_selection:
        # check all digits are the same
        numbers = first_selection[(pref,ext)]
        numbers_digit = [ len(d) for d in numbers ]
        numbers_int = [ int(d) for d in numbers ]
        if skip_0: numbers_int = filter(lambda x: x!=0, numbers_int)
        if len(set(numbers_digit)) != 1:
            pass

        template_str = "%(pref)s%(digit)s%(ext)s" % dict(pref= pref,
                                                         digit= "?"*numbers_digit[0],
                                                         ext= "."+ext if ext else ""
                                                         )

        min_frame, max_frame = min(numbers_int), max(numbers_int)

        group.append((template_str, min_frame, max_frame))

    return group
# group_img_files_template ()

def template_to_filenames(img_template, min_frame, max_frame):
    """
    e.g. template_to_filenames("aaa_????.img", 1, 150)
         you get: [ "aaa_0001.img", ...., "aaa_0150.img" ]
    """

    if "_??????.h5" in img_template: return [img_template.replace("_??????.h5", "_master.h5")]

    re_var = re.compile("\?+")

    # like: "????"
    var = re_var.search(img_template).group()
    repl = "%%.%dd" % len(var)

    # replace e.g. ???? => 0001
    filenames = [ img_template.replace(var, repl % d) for d in range(min_frame, max_frame+1) ]

    return filenames

# template_to_filenames()

def find_existing_files_in_template(img_template, min_frame, max_frame, datadir=None, check_compressed=True):
    # Change dirname
    if datadir is not None:
        img_template = os.path.join(datadir, os.path.basename(img_template))

    filenames = template_to_filenames(img_template, min_frame, max_frame)
    found_files = []

    for f in filenames:
        if os.path.isfile(f):
            found_files.append(f)
        elif check_compressed:
            for ext in (".bz2", ".gz", ".xz"):
                if os.path.isfile(f+ext):
                    found_files.append(f+ext)
                    break
    return found_files
# find_existing_files_in_template()

def is_dataset(img_template, min_frame, max_frame, _epsilon=1e-5, quiet=False):
    """
    Check whether image files belong to single dataset
     - Successive oscillation (missing frames OK)
     - Identical wavelength, distance
    """

    img_indexes = []
    wavelengths = []
    distances   = []
    start_angles = []
    end_angles  = []
    ang_widths = []
    filenames = []

    img_files = template_to_filenames(img_template, min_frame, max_frame)

    # Read header
    for i, f in enumerate(img_files):
        if os.path.isfile(f):
            try:
                im = XIO.Image(f)
            except Exception, ex:
                if not quiet:
                    print traceback.format_exc()
                return False

            start_angles.append( im.header["PhiStart"] )
            end_angles.append( im.header["PhiEnd"] )
            ang_widths.append( im.header["PhiWidth"] )
            wavelengths.append( im.header["Wavelength"] )
            distances.append( im.header["Distance"] )
            img_indexes.append(i+1)
            filenames.append(f)

    ret_val = True

    # Check wavl
    wavelength_set = set(wavelengths)
    if len(wavelength_set) > 1:
        if not quiet:
            print "Dataset has %d wavelengths:"%(len(wavelength_set)),
            print wavelength_set
        ret_val = False

    # Check distances
    distance_set = set(distances)
    if len(distance_set) > 1:
        if not quiet:
            print "Dataset has %d distances:"%(len(distance_set)),
            print distance_set
        ret_val = False

    # Check osc. angle widths
    width_set = set(ang_widths)
    if len(width_set) > 1:
        if not quiet:
            print "Dataset has %d osc widths:"%(len(width_set)),
            print width_set
        ret_val = False

    # Check whether continuous images
    first_start_angle = start_angles[0]
    first_image_index = img_indexes[0]
    osc_width = ang_widths[0]

    for i, (f, img_i, s_ang, e_ang) in enumerate(zip(filenames,img_indexes, start_angles, end_angles)):
        expected_s_ang = (img_i - first_image_index) * osc_width + first_start_angle
        if abs(s_ang - expected_s_ang) >= _epsilon:
            if not quiet:
                print os.path.basename(f),"is Non-continuous frame?",
                print s_ang, "!=", expected_s_ang, "diff= expected+", s_ang - expected_s_ang
            ret_val = False

    return ret_val

# is_dataset()

def takeout_datasets(img_template, min_frame, max_frame, _epsilon=1e-5,
                     check_wavelength=True, check_distance=True, check_oscwidth=True):
    img_indexes = []
    wavelengths = []
    distances   = []
    start_angles = []
    end_angles  = []
    ang_widths = []
    filenames = []

    img_files = template_to_filenames(img_template, min_frame, max_frame)

    # Read header
    for i, f in enumerate(img_files):
        if os.path.isfile(f):
            try:
                im = XIO.Image(f)
            except Exception, ex:
                print "Error on reading", f
                print traceback.format_exc()
                return []

            start_angles.append( im.header["PhiStart"] )
            end_angles.append( im.header["PhiEnd"] )
            ang_widths.append( im.header["PhiWidth"] )
            wavelengths.append( im.header["Wavelength"] )
            distances.append( im.header["Distance"] )
            img_indexes.append(min_frame+i)
            filenames.append(f)

    borders = [] # e.g. if 1: [:1]+[1:]

    for i in xrange(1, len(filenames)):
        # Check wavlength
        if check_wavelength and abs(wavelengths[i] - wavelengths[i-1]) > 0.001:
            borders.append(i)

        # Check distances
        if check_distance and abs(distances[i] - distances[i-1]) > 0.01:
            borders.append(i)

        # Check osc. angle widths
        if check_oscwidth and abs(ang_widths[i] - ang_widths[i-1]) > 0.0001:
            borders.append(i)

    # Check whether continuous images
    borders.append(len(filenames))
    ex_borders = []
    i = 0
    for b in borders:
        diffs = [0]
        for j in xrange(i+1, b):
            expected_s_ang = (img_indexes[j] - img_indexes[i]) * ang_widths[i] + start_angles[i]
            diffs.append(expected_s_ang - start_angles[j])

        for j in xrange(1, len(diffs)):
            if abs(diffs[j] - diffs[j-1]) > 0.002:
                ex_borders.append(i+j)
            
        i = b

    borders.extend(ex_borders)
    borders.sort()

    ranges = []
    i = 0
    for b in borders:
        ranges.append((img_indexes[i], img_indexes[b-1]))
        i = b

    return ranges
# takeout_datasets()

def find_data_sets(wdir, skip_symlinks=True, skip_0=False, split_hdf_miniset=True):
    """
    Find data sets in wdir
    """

    img_files = find_img_files(wdir, skip_symlinks=skip_symlinks)
    img_files.sort()

    h5files = filter(lambda x:x.endswith(".h5"), img_files)

    ret = []

    # group images
    group = group_img_files_template(img_files, skip_0=skip_0)
    print group

    for img_template, min_frame, max_frame in group:
        if min_frame == max_frame:
            continue

        for minf, maxf in takeout_datasets(img_template, min_frame, max_frame):
            ret.append([img_template, minf, maxf])

    for f in h5files:
        im = XIO.Image(f)

        if not split_hdf_miniset:
            ret.append([f.replace("_master.h5","_??????.h5"), 1, im.header["Nimages"]])
            continue

        for i in xrange(im.header["Nimages"]//im.header["Nimages_each"]+1):
            nr0, nr1 = im.header["Nimages_each"]*i+1, im.header["Nimages_each"]*(i+1)
            if nr1 > im.header["Nimages"]: nr1 = im.header["Nimages"]
            ret.append([f.replace("_master.h5","_??????.h5"), nr0, nr1])
            if nr1 == im.header["Nimages"]: break

    return ret
# find_data_sets()

if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        wdir = sys.argv[1]
    else:
        wdir = os.getcwd()

    for img_template, min_frame, max_frame in find_data_sets(wdir):
        print "NAME_TEMPLATE_OF_DATA_FRAMES= %s" % img_template
        print "DATA_RNAGE= %d %d" % (min_frame, max_frame)
        print
