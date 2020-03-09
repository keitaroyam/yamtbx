from yamtbx.dataproc.xds.xparm import XPARM

def run(xds_inp):
    xp = XPARM()
    xp.set_info_from_xdsinp_or_inpstr(xdsinp=xds_inp)
    xp.update_cell_based_on_axes()
    print xp.xparm_str().rstrip()

if __name__ == "__main__":
    import sys
    xds_inp = sys.argv[1]
    run(xds_inp)
