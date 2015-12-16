#!/usr/bin/env yamtbx.python
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""

from yamtbx.dataproc.xds.correctlp import CorrectLp
import iotbx.phil
from collections import OrderedDict
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

master_params_str="""\
plot = *ios *rmeas *ccano *cmpl sigano cchalf red
 .type = choice(multi=True)
 .help = What to plot
output = "plot.pdf"
 .type = path
rdataout = "for_R.dat"
 .type = path
"""

def run(params, args):

    ofs = open(params.rdataout, "w")
    ofs.write("name s2max variable value\n")

    for_plot = OrderedDict() 
    for p in params.plot:
        print "Preparing", p
        for_plot[p] = OrderedDict()

    trans_table = dict(ios="i_over_sigma",
                       rmeas="r_meas",
                       ccano="cc_ano",
                       cmpl="cmpl",
                       sigano="sig_ano",
                       cchalf="cc_half",
                       red="redundancy")


    for lpfile, label in ((args[2*i],args[2*i+1]) for i in xrange((len(args))//2)):
        lp = CorrectLp(lpfile)
        print label, lpfile, lp.space_group.info(), "anomalous=%s"%lp.anomalous_flag
        ofs.write("# %s %s %s anomalous=%s\n" % (label, lpfile, lp.space_group.info(), lp.anomalous_flag))
        plot_x = map(lambda x:1/x**2, lp.table["all"]["dmin"][:-1])
        for p in params.plot:
            plot_y = lp.table["all"][trans_table[p]][:-1]
            for_plot[p][label] = plot_x, plot_y

            for px, py in zip(plot_x, plot_y):
                ofs.write("%s %.5f %s %f\n" % (label, px, trans_table[p], py))

    fig, ax = plt.subplots()
    #plt.title("Comapring xds results")

    s2_formatter = lambda x,pos: "inf" if x == 0 else "%.2f" % (1./math.sqrt(x))

    for i, p in enumerate(params.plot):
        ax = plt.subplot(len(params.plot),1,i+1)
        for lab in for_plot[p]:
            plot_x, plot_y = for_plot[p][lab]
            ax.plot(plot_x, plot_y, label=lab, marker="o")
        ax.set_ylabel(trans_table[p])
        ax.xaxis.set_major_formatter(FuncFormatter(s2_formatter))

    plt.xlabel("Resolution [A]")
    #plt.legend()
    leg = plt.legend(loc='center left', bbox_to_anchor=(1,0.5), numpoints=1)
    fig.subplots_adjust(top=0.8)
    fig.savefig(params.output, bbox_extra_artists=(leg,), bbox_inches='tight')
    plt.show()
    
    print """

Instruction for R:

R
library(ggplot2)
d <- read.table("%s", h=T)
number_ticks <- function(n) {function(limits) pretty(limits, n)}
ggplot(d, aes(x=s2max, y=value, colour=factor(name))) + geom_point() + geom_line() + facet_grid(variable~., scale="free") + scale_x_continuous(label=function(x)sprintf("%%.2f", 1/sqrt(x)), breaks=number_ticks(10)) 


""" % params.rdataout
# run()

if __name__ == "__main__":
    import sys

    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args

    run(params, args)
