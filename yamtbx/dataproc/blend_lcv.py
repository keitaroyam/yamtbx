from __future__ import division
from __future__ import unicode_literals
import numpy

"""
Calculate LCV (linear cell variation in percentage) and aLCV (absolute LCV in Angstrom),
which are used in BLEND program (Foadi et al. 2013, Acta Cryst. D)
"""


def diagonals(cells):
    a, b, c, al, be, ga = [cells[:,i] for i in range(6)]

    Dab = numpy.sqrt(a**2 + b**2 -2.*a*b*numpy.cos(numpy.pi-numpy.deg2rad(ga)))
    Dbc = numpy.sqrt(b**2 + c**2 -2.*b*c*numpy.cos(numpy.pi-numpy.deg2rad(al)))
    Dca = numpy.sqrt(c**2 + a**2 -2.*c*a*numpy.cos(numpy.pi-numpy.deg2rad(be)))
    return Dab, Dbc, Dca
# diagonals()

def aldists(d):
    tmp = numpy.zeros(dtype=numpy.float, shape=(d.size,d.size))
    tmp[:,] = d
    adist = numpy.abs(tmp - tmp.transpose()) # redundant!! don't want to do this.. but for-for loops could be slow..
    ldist = adist / numpy.minimum(tmp, tmp.transpose())

    return adist, ldist
# aldists()

def calc_lcv(cells):
    Dab, Dbc, Dca = diagonals(numpy.array(cells))

    Mab, Lab = aldists(Dab)
    Mbc, Lbc = aldists(Dbc)
    Mca, Lca = aldists(Dca)


    # take max of M series
    iab = numpy.where(Mab==numpy.amax(Mab))
    ibc = numpy.where(Mbc==numpy.amax(Mbc))
    ica = numpy.where(Mca==numpy.amax(Mca))
    Lab = Lab[iab][0]
    Lbc = Lbc[ibc][0]
    Lca = Lca[ica][0]
    Mab = Mab[iab][0]
    Mbc = Mbc[ibc][0]
    Mca = Mca[ica][0]
    
    midx = numpy.argmax([Mab,Mbc,Mca])
    lcv = [Lab,Lbc,Lca][midx]
    alcv = [Mab,Mbc,Mca][midx]

    return lcv*100., alcv
# calc_lcv()
