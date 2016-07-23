#!/usr/bin/env phenix.python

# The original code is iotbx/command_line/emma.py
from __future__ import division

from iotbx import crystal_symmetry_from_any
from iotbx.option_parser import option_parser
from cctbx import euclidean_model_matching as emma
from iotbx.command_line.emma import get_emma_model
import sys, os, re

def get_emma_models_from_lst(file_name, crystal_symmetry):
  read_flag = False
  positions = []
  re_lst_header = re.compile("Try  *([0-9]+), CPU  *([0-9]+), CC All/Weak  *([-0-9\.]+)  */  *([-0-9\.]+)")
  
  for l in open(file_name):
    if l.startswith("    x       y       z"):
      read_flag = True
      positions = []
    elif read_flag and l.strip() == "":
      read_flag = False
    elif read_flag:
      site = map(float, (l[:8], l[8:16], l[16:24]))
      positions.append(emma.position(str(len(positions)+1), site))
    elif l.startswith(" Try "):
      r = re_lst_header.search(l)
      itry, ccall, ccweak = r.group(1), r.group(3), r.group(4)

      ret = emma.model(crystal_symmetry.special_position_settings(), positions)
      ret.label = "Try%s_CCall_%s_CCweak_%s" % (itry, ccall, ccweak)
      yield (ret, itry, ccall, ccweak)
# get_emma_models_from_lst


def run(args, command_name="emma_shelxd_lst.py"):
  command_line = (option_parser(
    usage=command_name + " [options]"
         +" reference_coordinates shelxd-lst-file",
    description="Example: %s model1.pdb sad_fa.lst" % command_name)
    .enable_symmetry_comprehensive()
    .option(None, "--tolerance",
      action="store",
      type="float",
      default=.5,
      help="match tolerance",
      metavar="FLOAT")
    .option(None, "--diffraction_index_equivalent",
      action="store_true",
      help="Use only if models are diffraction-index equivalent.")
  ).process(args=args, nargs=2)
  crystal_symmetry = command_line.symmetry
  if (   crystal_symmetry.unit_cell() is None
      or crystal_symmetry.space_group_info() is None):
    for file_name in command_line.args:
      crystal_symmetry = crystal_symmetry.join_symmetry(
        other_symmetry=crystal_symmetry_from_any.extract_from(
          file_name=file_name),
        force=False)

  tolerance = command_line.options.tolerance
  print "Tolerance:", tolerance
  if (tolerance <= 0.):
    raise ValueError, "Tolerance must be greater than zero."
  print
  diffraction_index_equivalent = \
    command_line.options.diffraction_index_equivalent
  if (diffraction_index_equivalent):
    print "Models are diffraction index equivalent."
    print
  emma_ref = get_emma_model(file_name=command_line.args[0],
                            crystal_symmetry=crystal_symmetry)

  emma_ref.show("Reference model")

  emma_others = get_emma_models_from_lst(command_line.args[1], crystal_symmetry)


  print "try CCall CCweak nmatch rms order.min order.max"
  
  for emma_other, itry, ccall, ccweak in emma_others:
    model_matches = emma.model_matches(model1=emma_ref,
                                       model2=emma_other,
                                       tolerance=tolerance,
                                       models_are_diffraction_index_equivalent=diffraction_index_equivalent)
    print itry, ccall, ccweak,

    if (model_matches.n_matches() == 0):
      print "0 nan nan nan"
    else:
      max_n_pairs = None
      first=True
      for match in model_matches.refined_matches:
        if (max_n_pairs is None or len(match.pairs) > max_n_pairs*0.2):
          orders = map(lambda x: int(x[1]), match.pairs)
          print "%3d %.5f %3d %3d" % (len(match.pairs), match.rms, min(orders), max(orders))
          #match.show()
          #first=False
          break
        if (max_n_pairs is None):
          max_n_pairs = len(match.pairs)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
