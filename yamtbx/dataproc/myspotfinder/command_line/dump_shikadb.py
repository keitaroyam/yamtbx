import  sqlite3
import cPickle as pickle
import datetime

def read_db(dbfile):
    con = sqlite3.connect(dbfile, timeout=10)

    try:
        c = con.execute("select dirname, time from updates")
        results = c.fetchall()
        if results:
            print "TABLE updates"
            for d, t in results:
                t = datetime.datetime.fromtimestamp(t)
                print "  %s %s" % (t.strftime("%Y-%m-%d %H:%M:%S.%f"), d)
            return
    except sqlite3.OperationalError:
         print "TABLE updates does not exist\n"

    print "TABLE spots"
    c = con.execute("select filename,spots from spots")
    results = dict(map(lambda x: (str(x[0]), pickle.loads(str(x[1]))), c.fetchall()))
    for filename in sorted(results):
        msg = results[filename]
        spots = msg["spots"]
        startt = datetime.datetime.fromtimestamp(msg["starttime"])
        print "  %s %s %d spots" % (startt.strftime("%Y-%m-%d %H:%M:%S.%f"), 
                                    filename, len(spots))
        #print "    %s" % msg["header"]

    print
    print "TABLE stats"
    c = con.execute("select imgf,nspot,total,mean from stats")
    results = dict(map(lambda x: (str(x[0]), x[1:]), c.fetchall()))
    for filename in sorted(results):
        stats = results[filename]
        print "  %s %s" % (filename, stats)

    print
    print "TABLE status"
    c = con.execute("select filename from status")
    results = c.fetchall()
    for filename in results:
        print "  %s" % (filename)


# read_db()

def run_from_args(argv):
    return read_db(argv[0])

if __name__ == "__main__":
    import sys
    #scanlog, shika_db = sys.argv[1:]

    run_from_args(sys.argv[1:])
