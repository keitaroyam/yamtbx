"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
import iotbx.phil
from cctbx import uctbx
from cctbx import sgtbx
from cctbx import crystal
from cctbx.crystal import reindex
from cctbx.sgtbx import pointgroup_tools
from yamtbx.dataproc.xds.xparm import XPARM
from yamtbx.dataproc.pointless import Pointless
from yamtbx.dataproc.xds import correctlp

import os
import sys
import networkx as nx
import numpy

master_params_str = """
topdir = None
 .type = path
xdsdir = None
 .type = path
 .multiple = true
 .help = Either topdir= or (multiple) xdsdir= should be specified.
tol_length = 0.1
 .type = float
 .help = relative_length_tolerance
tol_angle = 5
 .type = float
 .help = absolute_angle_tolerance in degree
do_pointless = False
 .type = bool
 .help = Run pointless for largest group data to determine symmetry
"""

class CheckMulti:
    def __init__(self, topdir=None, xdsdirs=None, out=sys.stdout):
        assert (topdir, xdsdirs).count(None) == 1

        self.topdir, self.xdsdirs = topdir, xdsdirs
        self.dirs, self.p1cells, self.symms = [], [], []
        self.reference_symmetries = []
        self.G = None
        self.cosets = {}
        self.groups = []
        self.grouped_dirs = []
        self.out = out
    # __init__()

    def get_symms_from_xds_results(self):
        if self.xdsdirs is not None:
            xdsdirs = filter(lambda x: os.path.isfile(os.path.join(x, "GXPARM.XDS")), self.xdsdirs)
        else:
            xdsdirs = map(lambda x: x[0], filter(lambda x: "GXPARM.XDS" in x[2], os.walk(self.topdir)))

        symms = []
        print >>self.out, "Idx Dir Cell P1Cell"
        idx = 0
        for root in xdsdirs:
            print >>self.out, "%.3d"%idx,
            print >>self.out, os.path.relpath(root, self.topdir) if self.topdir is not None else root,
            gxparm_xds = os.path.join(root, "GXPARM.XDS")
            p1cell = correctlp.get_P1_cell(os.path.join(root, "CORRECT.LP"))
            try:
                xparm = XPARM(gxparm_xds)
            except ValueError:
                print >>self.out, "Invalid xparm format:", gxparm_xds
                continue
            xs = xparm.crystal_symmetry()

            self.dirs.append(root)
            self.p1cells.append(p1cell)
            self.symms.append(xs)
            print >>self.out, xs.space_group_info(), xs.unit_cell(), p1cell
            idx += 1

        assert len(self.dirs) == len(self.symms) == len(self.p1cells)
    # get_cells_from_xds_results()

    def construct_graph(self, tol_length, tol_angle):
        self.G = nx.Graph()

        for i in xrange(len(self.p1cells)):
            self.G.add_node(i)

        for i in xrange(len(self.p1cells)-1):
            for j in xrange(i+1, len(self.p1cells)):
                cell_i = self.p1cells[i]
                cell_j = self.p1cells[j]
                if cell_i.is_similar_to(cell_j, tol_length, tol_angle):
                    self.G.add_edge(i, j)
                else:
                    cosets = reindex.reindexing_operators(crystal.symmetry(cell_i, 1),
                                                          crystal.symmetry(cell_j, 1),
                                                          tol_length, tol_angle)
                    self.cosets[(i,j)] = cosets
                    if cosets.double_cosets is not None:
                        self.G.add_edge(i, j)
        #nx.write_dot(self.G, "compatible_cell_graph.dot")
    # construct_graph()

    def _average_p1_cell(self, idxes):
        cells = [self.p1cells[idxes[0]].parameters()]
        for j in idxes[1:]:
            key = tuple(sorted((j, idxes[0])))
            if key in self.cosets and self.cosets[key].double_cosets is not None:
                #print "debug:: using cosets", self.p1cells[j].parameters()
                cbop = self.cosets[key].combined_cb_ops()[0]
                cells.append(self.p1cells[j].change_basis(cbop).parameters())
            else:
                cells.append(self.p1cells[j].parameters())

        cells = numpy.array(cells)
        return map(lambda i: cells[:,i].mean(), xrange(6))
    # _average_p1_cell()

    def group_xds_results(self, show_details=True):
        self.groups = map(lambda g: list(g), nx.connected_components(self.G))
        self.groups.sort(key=lambda x:-len(x))
        self.grouped_dirs = []
        self.reference_symmetries = []

        for i, sg in enumerate(self.groups):
            self.reference_symmetries.append([])
            avg_cell = uctbx.unit_cell(self._average_p1_cell(sg))
            print >>self.out, "[%2d]"%(i+1), len(sg), "members:"
            print >>self.out, " Averaged P1 Cell=", " ".join(map(lambda x:"%.2f"%x, avg_cell.parameters()))
            print >>self.out, " Members=", sg
            if show_details:
                # by explore_metric_symmetry
                print >>self.out, " Possible symmetries:"
                sg_explorer = pointgroup_tools.space_group_graph_from_cell_and_sg(avg_cell,  sgtbx.space_group_info("P1").group(), max_delta=10)
                tmp = []
                for obj in sg_explorer.pg_graph.graph.node_objects.values():
                    pg = obj.allowed_xtal_syms[0][0].space_group().build_derived_reflection_intensity_group(True).info()
                    cbop = obj.allowed_xtal_syms[0][1]
                    trans_cell = avg_cell.change_basis(cbop)

                    # If beta<90 in monoclinic system, force it to have beta>90
                    if pg.group().crystal_system() == "Monoclinic" and trans_cell.parameters()[4] < 90:
                        op = sgtbx.change_of_basis_op("-h,-k,l")
                        cbop = op * cbop
                        trans_cell = trans_cell.change_basis(op)

                    self.reference_symmetries[i].append((pg, trans_cell))
                    tmp.append(["   %-10s %s %s" % (pg, " ".join(map(lambda x:"%6.2f"%x, trans_cell.parameters())), cbop), pg.type().number()])
                for t in sorted(tmp, key=lambda x:x[1]):
                    print >> self.out, t[0]
                
                self.reference_symmetries[i].sort(key=lambda x:x[0].type().number())

                # by pointless results
                sgs = map(lambda x: self.symms[x].space_group(), sg)
                sgs_set = set(sgs)
                print >>self.out, " Space group and frequency in processed results:"
                for s in sorted(sgs_set, key=lambda x:x.type().number()):
                    print >>self.out, "   %-10s: %4d files" % (s.info(), sgs.count(s))
                print >>self.out, ""

            dirs = map(lambda x: self.dirs[x], sg)
            self.grouped_dirs.append(dirs)
    # group_xds_results()

    def get_reference_symm(self, group_idx, rs_idx):
        # XXX should be able to specify space group with screws

        if group_idx >= len(self.reference_symmetries):
            return None
        if rs_idx >= len(self.reference_symmetries[group_idx]):
            return None

        pg, cell = self.reference_symmetries[group_idx][rs_idx]

        return crystal.symmetry(cell,
                                space_group_info=pg,
                                assert_is_compatible_unit_cell=False)
    # get_reference_symm()

    def get_selectable_symms(self, group_idx):
        if group_idx >= len(self.reference_symmetries):
            return []

        return self.reference_symmetries[group_idx]
    # get_selectable_symms()


    def get_most_frequent_symmetry(self, group_idx):
        # Should call after self.group_xds_results()
        sym_freq = {}
        grp = self.groups[group_idx]
        sgs = map(lambda x: str(self.symms[x].space_group().build_derived_reflection_intensity_group(True).info()), grp)
        for s in sgs:
            sym_freq[s] = sym_freq.get(s,0) + 1

        if len(sym_freq) > 0:
            sym_freq = sym_freq.items()
            sym_freq.sort(key=lambda x:-x[1])
            if len(sym_freq) > 1 and sym_freq[0][0] == "P 1":
                return sym_freq[1][0]
            else:
                return sym_freq[0][0]
        else:
            return ""
    # get_most_frequent_symmetry()

def run(params, out=sys.stdout):
    cm = CheckMulti(topdir=params.topdir, xdsdirs=params.xdsdir, out=out)
    cm.get_symms_from_xds_results()
    cm.construct_graph(params.tol_length, params.tol_angle)
    cm.group_xds_results()
    print

    ret = cm.grouped_dirs

    if len(ret) == 0:
        return ret
    
    print >>out
    print >>out, "About the largest group:"
    for idx, wd in enumerate(ret[0]):
        gxparm_xds = os.path.join(wd, "GXPARM.XDS")
        correct_lp = os.path.join(wd, "CORRECT.LP")
        print >>out, "%.3d %s" % (idx, os.path.relpath(wd, params.topdir) if params.topdir is not None else wd),
        if not os.path.isfile(gxparm_xds):
            print >>out, "Unsuccessful"
            continue
        
        sg = sgtbx.space_group_info(XPARM(gxparm_xds).spacegroup)
        clp = correctlp.CorrectLp(correct_lp)
        if "all" in clp.table:
            cmpl = clp.table["all"]["cmpl"][-1]
        else:
            cmpl = float("nan")
        ISa = clp.a_b_ISa[-1]
        print >>out, "%10s ISa=%5.2f Cmpl=%5.1f " % (sg, ISa, cmpl)

    if params.do_pointless:
        worker = Pointless()
        files = map(lambda x: os.path.join(x, "INTEGRATE.HKL"), ret[0])
        #print files
        files = filter(lambda x: os.path.isfile(x), files)
        
        print >>out, "\nRunning pointless for the largest member."
        result = worker.run_for_symm(xdsin=files, 
                                  logout="pointless.log",
                                  tolerance=10, d_min=5)
        if "symm" in result:
            print >>out, " pointless suggested", result["symm"].space_group_info()


    if 0:
        import pylab
        pos = nx.spring_layout(G)
        #pos = nx.spectral_layout(G)
        #pos = nx.circular_layout(G)

        #nx.draw_networkx_nodes(G, pos, node_size = 100, nodelist=others, node_color = 'w')
        nx.draw_networkx_nodes(G, pos, node_size = 100, node_color = 'w')
        nx.draw_networkx_edges(G, pos, width = 1)
        nx.draw_networkx_labels(G, pos, font_size = 12, font_family = 'sans-serif', font_color = 'r')

        pylab.xticks([])
        pylab.yticks([])
        pylab.savefig("network.png")
        pylab.show()

    return cm
# run()

def run_from_args(argv):
    cmdline = iotbx.phil.process_command_line(args=argv,
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args
    
    for arg in args:
        if os.path.isdir(arg) and params.topdir is None:
            params.topdir = arg

    if params.topdir is None:
        params.topdir = os.getcwd()

    run(params)
# run_from_args()

if __name__ == "__main__":
    import sys
    run_from_args(sys.argv[1:])
