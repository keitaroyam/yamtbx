"""
(c) RIKEN 2017. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from dxtbx.model.experiment_list import ExperimentListFactory
from dxtbx.model.experiment_list import ExperimentListDumper
from dxtbx.model import Crystal
from yamtbx.dataproc.xds import integrate_hkl_as_flex
from yamtbx.dataproc.xds.idxreflp import SpotXds
from yamtbx.dataproc.xds.xparm import prep_xparm_objects_from_integrate_lp
from dials.array_family import flex
from dials.algorithms.centroid import centroid_px_to_mm_panel
import os
import copy

def import_integrated(integrate_hkl, min_ios=3):
    reader = integrate_hkl_as_flex.reader(integrate_hkl, "IOBS,SIGMA,XCAL,YCAL,ZCAL,XOBS,YOBS,ZOBS,ISEG".split(","))

    # reference: dials/command_line/import_xds.py
    table = flex.reflection_table()
    table["id"] = flex.int(len(reader.hkl), 0)
    table["panel"] = flex.size_t(len(reader.hkl), 0) # only assuming single panel
    table["miller_index"] = reader.hkl
    table["xyzcal.px.value"] = flex.vec3_double(reader.data["XCAL"], reader.data["YCAL"], reader.data["ZCAL"])
    table["xyzobs.px.value"] = flex.vec3_double(reader.data["XOBS"], reader.data["YOBS"], reader.data["ZOBS"])
    table["intensity.cor.value"] = reader.data["IOBS"]
    table["intensity.cor.sigma"] = reader.data["SIGMA"] # not valid name, just for making selection
    table["flags"] = flex.size_t(len(table), table.flags.indexed | table.flags.strong)

    table = table.select(table["intensity.cor.sigma"] > 0)
    table = table.select(table["intensity.cor.value"]/table["intensity.cor.sigma"] >= min_ios)
    table = table.select(table["xyzobs.px.value"].norms() > 0) # remove invalid xyzobs

    # dummy
    table["xyzobs.px.variance"] = flex.vec3_double(len(table), (1,1,1)) # TODO appropriate variance value
    table["s1"] = flex.vec3_double(len(table), (0,0,0)) # will be updated by set_obs_s1()

    del table["intensity.cor.value"]
    del table["intensity.cor.sigma"]

    return table
# import_integrated()

def import_spot_xds(spot_xds):
    sx = SpotXds(spot_xds)
    spots = filter(lambda x: x[-1][0] is not None and not (x[-1][0]==x[-1][1]==x[-1][2]==0), sx.items)

    # reference: dials/command_line/import_xds.py
    table = flex.reflection_table()
    table["id"] = flex.int(len(spots), 0)
    table["panel"] = flex.size_t(len(spots), 0) # only assuming single panel
    table["miller_index"] = flex.miller_index(map(lambda x: x[-1], spots))
    table["xyzobs.px.value"] = flex.vec3_double(map(lambda x: x[0], spots))
    table["flags"] = flex.size_t(len(table), table.flags.indexed | table.flags.strong)
    

    # dummy
    table["xyzobs.px.variance"] = flex.vec3_double(len(table), (1,1,1)) # TODO appropriate variance value
    table["s1"] = flex.vec3_double(len(table), (0,0,0)) # will be updated by set_obs_s1()

    return table
# import_spot_xds()

def px_to_mm(experiment, table):
    # reference: dials/algorithms/indexing/indexer.py map_spots_pixel_to_mm_rad()
    centroid_position, centroid_variance, _ = centroid_px_to_mm_panel(experiment.detector[0], experiment.scan,
                                                                      table['xyzobs.px.value'],
                                                                      table['xyzobs.px.variance'],
                                                                      flex.vec3_double(len(table), (1,1,1)))
    table['xyzobs.mm.value'] = centroid_position
    table['xyzobs.mm.variance'] = centroid_variance
# px_to_mm()

def run(xds_inp, xparm=None, integrate_lp=None, integrate_hkl=None, spot_xds=None,
        space_group=None, reindex_op=None, out_prefix=None, out_dir="."):
    out_prefix = out_prefix+"_" if out_prefix else ""

    if integrate_lp is not None:
        xparm_objs = prep_xparm_objects_from_integrate_lp(integrate_lp, xparm_ref=xparm)
        rr, xp = xparm_objs[0]
        xparm = os.path.join(os.path.dirname(xds_inp), "XPARM.XDS_%.6d-%.6d"%rr)
        open(xparm, "w").write(xp.xparm_str())

    # FIXME template of experiment.imageset could be wrong when relative path
    #       and ######.h5 may need to be replaced with master.h5
    experiments = ExperimentListFactory.from_xds(xds_inp, xparm) # XDS.INP needed for sweep info

    assert len(experiments) == 1
    experiment = experiments[0]

    # I don't know what invalid X/Y/ZOBS values should be when multi-panel detector
    assert len(experiment.detector) == 1 

    if None not in (space_group, reindex_op):
        cryst_orig = copy.deepcopy(experiment.crystal)
        cryst_reindexed = cryst_orig.change_basis(reindex_op)
        a, b, c = cryst_reindexed.get_real_space_vectors()
        cryst_reindexed = Crystal(a, b, c, space_group=space_group)
        experiment.crystal.update(cryst_reindexed)

    # Very dirty fix.. but no way to change template after object creation??
    json_str = ExperimentListDumper(experiments).as_json().replace("_######.h5", "_master.h5")
    open(os.path.join(out_dir, out_prefix+"experiments.json"), "w").write(json_str)

    if integrate_hkl is not None:
        table = import_integrated(integrate_hkl)
        px_to_mm(experiment, table)
        if None not in (space_group, reindex_op): table["miller_index"] = reindex_op.apply(table["miller_index"])
        table.as_pickle(os.path.join(out_dir, out_prefix+"integrate_hkl.pickle"))

    if spot_xds is not None:
        table = import_spot_xds(spot_xds)
        px_to_mm(experiment, table)
        if None not in (space_group, reindex_op): table["miller_index"] = reindex_op.apply(table["miller_index"])
        table.as_pickle(os.path.join(out_dir, out_prefix+"spot_xds.pickle"))

# run()

if __name__ == "__main__":
    import sys
    xds_inp, xparm, integrate_hkl = sys.argv[1:]

    """
    from cctbx import sgtbx
    space_group=sgtbx.space_group_info("c2").group()
    reindex_op=sgtbx.change_of_basis_op("h+2*l,h,k"))
    """

    run(xds_inp, integrate_lp=xparm, integrate_hkl=integrate_hkl)
