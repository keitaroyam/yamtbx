"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
import os, sys, glob, subprocess

def make_xy_hkl_list(hklfile):

	xy_hkl = {}
	for l in open(hklfile):
		if l.startswith("!"): continue

		H,K,L,IOBS,SIGMA,XCAL,YCAL,ZCAL,RLP,PEAK,CORR,MAXC,XOBS,YOBS,ZOBS,ALF0,BET0,ALF1,BET1,PSI,ISEG = l.split()
		
		#reflections.setdefault((H,K,L), []).append((num, XOBS, YOBS, ZOBS))
		#xy_hkl[(float(XOBS), float(YOBS))] = (int(H), int(K), int(L) ) 
		xy_hkl[(float(XCAL), float(YCAL))] = (int(H), int(K), int(L) ) 
		
	return xy_hkl
# make_xy_hkl_list()

def make_adx_lines(xy_hkl):
	s = ""
	for xy, hkl in xy_hkl.items():
		s += "%d %d %d %d %d\n"%(xy[0], xy[1], hkl[0], hkl[1], hkl[2])

	return s
# make_adx_lines()

def make_adxfile(hklfile):
	out = os.path.splitext(hklfile)[0] + ".adx"
	ofs = open(out, "w")
	xy_hkl = make_xy_hkl_list(hklfile)
	ofs.write(make_adx_lines(xy_hkl))
	return out
# make_adxfile()

if __name__ == "__main__":


	hklfile = sys.argv[1]
	adxfile = hklfile.replace("INTEGRATE_", "FRAME_").replace(".HKL", ".adx")

	xy_hkl = make_xy_hkl_list(hklfile)

	ofs = open(adxfile, "w")
	ofs.write(make_adx_lines(xy_hkl))

