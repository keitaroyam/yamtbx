import cPickle as pickle
import os

def run(streamin, pklin, key, streamout=None):
    if streamout is None:
        streamout = os.path.splitext(os.path.basename(streamin))[0] + "_sort_%s.stream" % key

    assert key[-1] in ("+","-") # + for increasing order, - for decreasing order
    rev_order = key[-1] == "-"
    key = key[:-1]

    stats = pickle.load(open(pklin))
    sorted_indices = sorted(range(len(stats["chunk_ranges"])),
                            key=lambda x: stats[key][x],
                            reverse=rev_order)

    ifs = open(streamin)
    ofs = open(streamout, "w")
    ofs.write(ifs.readline())

    for idx in sorted_indices:
        print "writing", stats[key][idx]
        s, e = stats["chunk_ranges"][idx]
        ifs.seek(s-1)
        ofs.write(ifs.read(e-s+1))

# run()

if __name__ == "__main__":
    import sys
    run(sys.argv[1], sys.argv[2], sys.argv[3])
    
