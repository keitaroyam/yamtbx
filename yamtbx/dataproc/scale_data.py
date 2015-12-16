"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
import scipy.optimize
import numpy
from cctbx.array_family import flex
import sys

class kBdecider3:
    def __init__(self, obs, calc):
        self.obs, self.calc = obs.common_sets(calc, assert_is_similar_symmetry=False)
    # __init__()

    def run(self):
        from mmtbx.scaling import absolute_scaling
        iso_scale_and_b_obs = absolute_scaling.ml_iso_absolute_scaling(miller_array=self.obs,
                                                                       n_residues=100)
        iso_scale_and_b_calc = absolute_scaling.ml_iso_absolute_scaling(miller_array=self.calc,
                                                                       n_residues=100)
        b_obs = iso_scale_and_b_obs.b_wilson
        b_calc = iso_scale_and_b_calc.b_wilson

        B = -2*(b_calc - b_obs)
        k = kBdecider2.get_linear_scale(self.obs, self.calc, B)
        print "k,B=",k,B
        return k, B

# Does not work. why?
class kBdecider2:
    def __init__(self, obs, calc):
        self.obs, self.calc = obs.common_sets(calc, assert_is_similar_symmetry=False)

        self.x = numpy.array([0.]) # B
    # __init__()

    @staticmethod
    def get_linear_scale(obs, calc, B):
        #return flex.sum(obs*calc) / flex.sum(flex.pow2(calc))
        o = obs.data()
        c = flex.exp(-B*calc.d_star_sq().data())*calc.data()
        return flex.sum(o*c) / flex.sum(flex.pow2(c))

    def callback(self, xk):
        pass
    # callback()

    def run(self):
        f = self.f(self.x)
        print "# initial  f= %.6e" % f
        print "# initial  x=", self.x

        status = scipy.optimize.minimize(fun=self.f,
                                         x0=self.x,
                                         method="CG",
                                         callback=self.callback)
        self.x = status.x

        f = self.f(self.x)
        print "#   final  f= %.6e" % f
        print "#   final  x=", self.x
        B = float(self.x[0])
        return self.get_linear_scale(self.obs, self.calc, B), B
    # run()

    def f(self, x):
        print x
        B = float(x[0])
        #d_star_sq = self.calc.d_star_sq().data()
        #obs = self.obs.data()
        #calc = flex.exp(-B*d_star_sq)*self.calc.data()
        k = self.get_linear_scale(self.obs, self.calc, B)
        return flex.sum(flex.pow2(self.obs.data() - k*self.calc.data()))
    # f()

# class kBdecider2


class kBdecider:
    def __init__(self, obs, calc):
        self.obs, self.calc = obs.common_sets(calc, assert_is_similar_symmetry=False)

        ini_scale = flex.sum(self.obs.data()*self.calc.data()) / flex.sum(flex.pow2(self.calc.data()))

        self.x = numpy.array([ini_scale, 0.]) # k, B
    # __init__()

    def callback(self, xk):
        pass
    # callback()

    def run(self):
        f, df = self.f(self.x), self.df(self.x)
        print "# initial  f= %.6e" % f
        print "# initial df=", tuple(df)

        status = scipy.optimize.minimize(fun=self.f,
                                         x0=self.x,
                                         method="L-BFGS-B",
                                         jac=self.df,
                                         callback=self.callback)
        self.x = status.x

        f, df = self.f(self.x), self.df(self.x)
        print "#   final  f= %.6e" % f
        print "#   final df=", tuple(df)
        print "#   final  x=", self.x

        return float(self.x[0]), float(self.x[1])
    # run()

    def f(self, x):
        k, B = float(x[0]), float(x[1])
        d_star_sq = self.calc.d_star_sq().data()
        return flex.sum(flex.pow2(self.obs.data() - k*flex.exp(-B*d_star_sq)*self.calc.data()))
    # f()

    def df(self, x):
        k, B = float(x[0]), float(x[1])
        d_star_sq = self.calc.d_star_sq().data()
        tmp = self.obs.data() - k * flex.exp(-B*d_star_sq) * self.calc.data()
        dfdk = flex.sum(-2. * tmp * flex.exp(-B*d_star_sq) * self.calc.data())
        dfdB = flex.sum(2. * tmp * k * d_star_sq * flex.exp(-B*d_star_sq) * self.calc.data())
        return numpy.array([dfdk, dfdB])
    # df()
# class kBdecider

def calc_R(obs, calc, do_scale=True):
    #obs, calc = obs.common_sets(calc, assert_is_similar_symmetry=False)
    if do_scale:
        scale = flex.sum(obs.data()*calc.data()) / flex.sum(flex.pow2(calc.data()))
    else:
        scale = 1.

    R = flex.sum(flex.abs(obs.data() - scale*calc.data())) / flex.sum(0.5 * obs.data() + 0.5 * scale*calc.data())

    return R, scale
# calc_CC()

def dump_R_in_bins(obs, calc, scale_B=True, log_out=sys.stdout, n_bins=20):
    #obs, calc = obs.common_sets(calc, assert_is_similar_symmetry=False)

    if scale_B:
        scale, B = kBdecider(obs, calc).run()
        d_star_sq = calc.d_star_sq().data()
        calc = calc.customized_copy(data = scale * flex.exp(-B*d_star_sq) * calc.data())

    binner = obs.setup_binner(n_bins=n_bins)
    count=0
    log_out.write("dmax - dmin: R (nref) <I1> <I2> scale\n")

    for i_bin in binner.range_used():
        tmp_obs = obs.select(binner.bin_indices() == i_bin)
        tmp_calc = calc.select(binner.bin_indices() == i_bin)

        low = binner.bin_d_range(i_bin)[0]
        high = binner.bin_d_range(i_bin)[1]

        if scale_B:
            scale = 1.
        else:
            scale = flex.sum(tmp_obs.data()*tmp_calc.data()) / flex.sum(flex.pow2(tmp_calc.data()))

        R = flex.sum(flex.abs(tmp_obs.data() - scale*tmp_calc.data())) / flex.sum(0.5 * tmp_obs.data() + 0.5 * scale*tmp_calc.data())

        log_out.write("%5.2f - %5.2f: %.5f (%d) %.1f %.1f %.3e\n" % (low, high, R, len(tmp_obs.data()),
                                                                 flex.mean(tmp_obs.data()), flex.mean(tmp_calc.data()),
                                                                 scale))

    log_out.write("Overall R = %.5f (scale=%.3e, %%comp=%.3f)\n\n" % (calc_R(obs, calc, do_scale=not scale_B) + (obs.completeness()*100.,)) )

# dump_CC_in_bins()
