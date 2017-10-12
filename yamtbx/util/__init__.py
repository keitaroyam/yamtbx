"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
import os
import sys
import re
import shutil
import subprocess
import commands
import glob
import tempfile
from libtbx.utils import null_out
import libtbx.load_env

def call(cmd, arg="",
         stdin=None, stdout=subprocess.PIPE,
         wdir=None,
         expects_in=[], expects_out=[]):

    ##
    # call the external program using subprocess.
    #
    # @param expects_in expected files before running
    # @param expects_out expected files after running
    #
    # expected_in/out must be written as relative path from wdir or absolute path.
    #

    def check_exist(files):
        is_exist = [os.path.isfile(f) for f in files]

        if sum(is_exist) != len(is_exist):
            not_founds = [ f for f, e in zip(files, is_exist) if not e ]
            raise Exception("Expected file(s) not found: " + " ".join(not_founds))

    # check_exist()

    if wdir is None:
        wdir = os.getcwd()

    # Go to working directory
    cwd = os.getcwd()
    os.chdir(wdir)

    # check before run
    check_exist(expects_in)

    # call the program
    p = subprocess.Popen("%s %s" % (cmd, arg),
                         shell=True,
                         stdin=subprocess.PIPE,
                         stdout=stdout,
                         stderr=stdout
                         )

    if stdin is not None:
        p.stdin.write(stdin)    

    if stdout == subprocess.PIPE:
        out, err = p.communicate()
    else:
        out, err = None, None
        p.stdin.close()
        p.wait()

    if p.returncode < 0:
        print >>sys.stderr, cmd, ": returncode is", p.returncode

    # check after run
    check_exist(expects_out)

    # go back to the previous working directory
    os.chdir(cwd)

    return p.returncode, out, err
# call()

def rotate_file(filename, copy=False):
    """
    Rotate file like logrotate.
    If given filename already exists, rename it to "filename".n, n=1...
    Filename with larger n is older one.
    """

    # If not exist,
    if not os.path.isfile(filename):
        return

    # make list [ [filename, number], ... ]
    old_list = []
    dot_files = glob.glob(filename + ".*")
    for f in dot_files:
        suffix = f.replace(filename+".", "")
        try:
            i = int(suffix)
            if str(i) == suffix: # ignore if suffix was such as 003...
                old_list.append([f, i])
        except ValueError, e:
            continue

    old_list.sort(lambda x,y: x[1]-y[1])

    # rotate files
    for f, i in reversed(old_list):
        os.rename(f, "%s.%d" % (f[:f.rfind(".")], i+1))

    if copy:
        shutil.copyfile(filename, filename + ".1")
    else:
        os.rename(filename, filename + ".1")

    return filename + ".1"
# rotate_file()

def safe_copy(src, dst, move=False):
    """
    Don't reveal file before copy completed.
    """

    src_name = os.path.basename(src)
    if os.path.isdir(dst): dst = os.path.join(dst, src_name)

    tmpfd, tmp = tempfile.mkstemp(prefix="."+src_name, dir=os.path.dirname(dst))
    os.close(tmpfd)

    shutil.copy2(src, tmp)
    os.rename(tmp, dst)

    if move and os.path.isfile(dst) and not os.path.islink(dst) and os.path.getsize(src)==os.path.getsize(dst):
        os.remove(src)

# safe_copy()

def commonalize(Is):
    new_Is = []
    Is0 = Is[0]
    for I in Is[1:]:
        Is0, I = Is0.common_sets(I, assert_is_similar_symmetry=False)
        new_Is.append(I)

    Is = []

    for I in new_Is:
        I = I.common_set(Is0, assert_is_similar_symmetry=False)
        assert len(Is0.data()) == len(I.data())
        Is.append(I)

    return [Is0,] + Is
# commonalize()

def get_number_of_processors(default=4):
    nproc = default

    if os.path.isfile("/proc/cpuinfo"):
        nproc = len(filter(lambda x:x.startswith("processor"), open("/proc/cpuinfo")))
    else:
        try:
            nproc = int(commands.getoutput("sysctl -n hw.ncpu"))
        except:
            pass

    return nproc
# get_number_of_processors()

def safe_float(v):
    try:
        return float(v)
    except ValueError:
        return float("nan")
# safe_float()

def num_th_str(v):
    s = str(v)
    if s[-1] == "1": return s+"st"
    if s[-1] == "2": return s+"nd"
    if s[-1] == "3": return s+"rd"
    return s+"th"
# num_th_str()

def directory_included(path, topdir, include_dir=[], exclude_dir=[]):
    l1 = filter(lambda x: x, path.split(os.sep))
    l2 = filter(lambda x: x, topdir.split(os.sep))

    lc = os.path.commonprefix([l1,l2])

    if len(lc) != len(l2): return False

    if include_dir == exclude_dir == []:
        return True

    if include_dir != []:
        for d in include_dir:
            if directory_included(path, d): return True
        return False

    if exclude_dir != []:
        for d in exclude_dir:
            if directory_included(path, d): return False
        return True
# directory_included()

def read_path_list(lstin, comment_strs=["#"], only_exists=False, as_abspath=False, err_out=null_out()):
    ret = []

    for l in open(lstin):
        for c in comment_strs:
            if c in l: l = l[:l.index(c)]
        
        l = l.strip()
        if not l: continue
        if only_exists and not os.path.exists(l):
            err_out.write("Error: file not found: %s\n"%l)
            continue
        ret.append(os.path.abspath(l) if as_abspath else l)

    return ret
# read_path_list()

def return_first_found_file(files, wd=None):
    for f in files:
        if wd is not None: f = os.path.join(wd, f)
        if os.path.isfile(f): return f
# return_first_found_file()

def expand_wildcard_in_list(fdlst, err_out=null_out()):
    ret = []
    for d in fdlst:
        gd = glob.glob(d)
        
        if len(gd) == 0:
            print >>err_out, "Error: No match!!: %s" % d
            continue
        ret.extend(gd)
    return ret
# expand_wildcard_in_list()

def check_disk_free_bytes(d):
    try:
        x = os.statvfs(d)
        return x.f_frsize * x.f_bavail
    except:
        return -1
# check_disk_free_bytes()

def get_temp_local_dir(prefix, min_bytes=None, min_kb=None, min_mb=None, min_gb=None, additional_tmpd=None):
    assert (min_bytes, min_kb, min_mb, min_gb).count(None) >= 2

    min_free_bytes = 0

    if min_bytes is not None: min_free_bytes = min_bytes
    if min_kb is not None: min_free_bytes = min_kb * 1024
    if min_mb is not None: min_free_bytes = min_mb * 1024**2
    if min_gb is not None: min_free_bytes = min_gb * 1024**3

    ramdisk = "/dev/shm"

    if os.path.isdir(ramdisk): tmpdirs = [ramdisk, tempfile.gettempdir()]
    else: tmpdirs = [tempfile.gettempdir()]

    if type(additional_tmpd) is str:
        tmpdirs.append(additional_tmpd)
    elif type(additional_tmpd) in (list, tuple):
        tmpdirs.extend(additional_tmpd)

    for tmpdir in tmpdirs:
        if check_disk_free_bytes(tmpdir) >= min_free_bytes:
            return tempfile.mkdtemp(prefix=prefix, dir=tmpdir)

    return None
# get_temp_local_dir()

def get_temp_filename(prefix="tmp", suffix="", wdir=None):
    tmpfd, tmp = tempfile.mkstemp(prefix=prefix, suffix=suffix, dir=wdir)
    os.close(tmpfd)
    return tmp
# get_temp_filename()

def replace_forbidden_chars(filename, repl="-"):
    return re.sub(r"[/><\*\\\?%:]", repl, filename)
# replace_forbidden_chars()

def human_readable_bytes(bytes):
    if bytes < 1024:
        return bytes, "B"
    elif bytes < 1024**2:
        return bytes/1024., "KB"
    elif bytes < 1024**3:
        return bytes/1024.**2, "MB"
    elif bytes < 1024**4:
        return bytes/1024.**3, "GB"
    elif bytes < 1024**5:
        return bytes/1024.**4, "TB"
    else:# if bytes < 1024**6:
        return bytes/1024.**5, "PB"
    
# human_readable_bytes()

def yamtbx_module_root():
    """
    Possible locations: modules/yamtbx or modules/yamtbx/yamtbx
    """
    tmp = libtbx.env.find_in_repositories("yamtbx/yamtbx")
    if tmp: return tmp
    tmp = libtbx.env.find_in_repositories("yamtbx")
    if tmp: return tmp
# yamtbx_module_root()
