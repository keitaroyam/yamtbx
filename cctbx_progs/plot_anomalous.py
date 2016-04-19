from cctbx.eltbx import sasaki, henke
import iotbx.phil

master_params_str = """
  element = None
    .type = str
  energy_range = [6000,20000]
    .type = floats
    .help = range of energy in eV
  energy_step = 1.
    .type = float
    .type = step in eV
"""

def run(params):
    assert len(params.energy_range) == 2
    assert params.energy_range[0] < params.energy_range[1]
    er = params.energy_range
    es = params.energy_step

    print "lambda eV fp fdp table"

    for table in ["sasaki", "henke"]:
        for ev in (er[0]+es*i for i in xrange(int((er[1]-er[0])/es+1))):
            t = eval(table).table(params.element)
            f = t.at_ev(ev)

            fp = f.fp() if f.is_valid_fp() else float("nan")
            fdp = f.fdp() if f.is_valid_fdp() else float("nan")

            print "%.4f %.1f %.3f %.3f %s" % (12398.4/ev, ev, fp, fdp, table)

# run()
if __name__ == "__main__":
    import sys
    cmdline = iotbx.phil.process_command_line(args=sys.argv[1:],
                                              master_string=master_params_str)
    params = cmdline.work.extract()
    args = cmdline.remaining_args

    run(params)
