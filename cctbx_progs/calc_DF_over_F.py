import iotbx.file_reader
from cctbx.array_family import flex

def run(hklin):
    arrays = iotbx.file_reader.any_file(hklin).file_server.miller_arrays
    for arr in arrays:
        if not arr.anomalous_flag():
            continue

        print arr.info()
        if arr.is_complex_array():
            arr = arr.as_amplitude_array() # must be F

        ano = arr.anomalous_differences()
        ave = arr.average_bijvoet_mates()
        
        ano, ave = ano.common_sets(ave)
        print "   <d''/mean>=", flex.mean(flex.abs(ano.data()) / ave.data())
        print " <d''>/<mean>=", flex.mean(flex.abs(ano.data())) / flex.mean(ave.data())
        print

if __name__ == "__main__":
    import sys
    run(sys.argv[1])
