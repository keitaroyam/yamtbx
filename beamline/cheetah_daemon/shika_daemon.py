"""
Usage:

Copy this script and start_shika.sh to secure location (only root can edit)
Setup /etc/init.d/shika_cheetah_daemon

From client:
wget -O - http://ramchan:8080/start/$(id -u)/$(id -g) 2>/dev/null
wget -O - http://ramchan:8080/stop 2>/dev/null
"""

import daemon.runner
import lockfile
import web
import psutil
import signal
import os
import sys
import subprocess
import pwd
import time

EXEC_SH = "/usr/local/cheetah_daemon/start_shika.sh"
PID_FILE = "/var/run/shikaeiger_sh.pid"

def get_all_children(pid, ret):
    children = psutil.Process(pid).children()
    ret.extend(map(lambda x: x.pid, children))
    for cp in children:
        get_all_children(cp.pid, ret)
# get_all_children()

def change_user(uid, gid):
    os.setgid(gid)
    os.setuid(uid)

    pw = pwd.getpwuid(uid)
    os.environ["HOME"] = pw.pw_dir
    os.environ["USER"] = pw.pw_name
    os.environ["LOGNAME"] = pw.pw_name
# set_uidgid()

class start_program:
    def GET(self, uid, gid):
        web.header("Content-type","text/plain")
        uid, gid = map(int, (uid, gid))

        # Security
        if uid < 500 or gid < 500:
            yield "Invalid UID (%d) or GID (%d)\n" % (uid, gid)
            return

        try: pwd.getpwuid(uid)
        except KeyError:
            yield "UID (%d) does not exist\n" % uid
            return

        for k in stop_program().GET(): yield k
        print "Starting program with %d/%d (now= %s)" % (uid, gid, time.ctime())
        yield "Starting program with %d/%d\n" % (uid, gid)
        #p = subprocess.Popen(EXEC_SH, shell=True,
        #                     preexec_fn=lambda: change_user(uid,gid)) # This fails when running as daemon

        rpipe, wpipe = os.pipe() # Reference: http://ameblo.jp/oyasai10/entry-10615215673.html
        pid = os.fork()
        if pid == 0: # Child
            os.close(rpipe)
            wpipe = os.fdopen(wpipe, "w")
            change_user(uid,gid)
            p = subprocess.Popen(EXEC_SH, shell=True)
            wpipe.write("%d\n"%p.pid)
            sys.exit()
        else: # Parent
            os.close(wpipe)
            rpipe = os.fdopen(rpipe, "r")
            pid = int(rpipe.readline().strip())
            open(PID_FILE, "w").write("%d"%pid)
            os.wait() # Wait child
    # GET()
# class start_program

class stop_program:
    def GET(self, from_web=True):
        if from_web: web.header("Content-type","text/plain")
        yield "Stop detected\n"

        if not os.path.isfile(PID_FILE):
            return
        
        pid = int(open(PID_FILE).read())
        print "Killing children of %d (now= %s)" % (pid, time.ctime())
        yield "Killing children of %d\n" % pid

        # Check process exists
        try: os.kill(pid, 0)
        except OSError:
            print "PID %d does not exist" % pid
            yield "PID %d does not exist\n" % pid
            return

        # Check process name
        proc = psutil.Process(pid)
        print "  proc name of %d = '%s' (user: %s)" % (pid, proc.name(), proc.username())
        if proc.name() != "start_shika.sh":
            print "  Not start_shika.sh! Quit."
            return
        # Kill all children
        children = [pid]
        get_all_children(pid, children)
        for p in reversed(children):
            try:
                proc = psutil.Process(p)
                print "  killing %d (name= %s, user= %s, ctime= %s)" % (p, proc.name(), proc.username(), time.ctime(proc.create_time()))
            except psutil.NoSuchProcess: continue
            os.kill(p, signal.SIGKILL)
        os.kill(pid, signal.SIGKILL)

        print "All child processes killed (now= %s)" % time.ctime()
        yield "All child processes killed\n"
    # GET()
# class stop_program

class App:
    # Reference: http://www.gavinj.net/2012/06/building-python-daemon-process.html
    def __init__(self):
        self.stdin_path = '/dev/null'
        self.stdout_path = '/dev/null'
        self.stderr_path = '/dev/null'
        self.pidfile_path =  "/var/run/shikaeiger.pid"
        self.pidfile_timeout = 5
    # __init__()

    def run(self):
        app = web.application(("/start/(.*)/(.*)", "start_program",
                               "/stop", "stop_program"),
                              dict(start_program=start_program, stop_program=stop_program))
        #app.run() # This sees sys.argv
        web.httpserver.runsimple(app.wsgifunc(), ("0.0.0.0", 8080))
    # run()
# class App

def terminate_daemon(signum, frame):
    print "Signal:", signum, frame
    try:
        for k in stop_program().GET(from_web=False): print k
    finally:
        sys.exit()
# terminate_daemon()

def run():
    app = App()
    daemon_runner = daemon.runner.DaemonRunner(app)
    daemon_runner.daemon_context.signal_map = {signal.SIGTERM: terminate_daemon}
    out = open("/var/log/shikaeiger.log", "a")
    daemon_runner.daemon_context.stdout = out
    daemon_runner.daemon_context.stderr = out
    #daemon_runner.daemon_context.files_preserve = [important_file, interesting_file]

    daemon_runner.do_action()
# run()

if __name__ == "__main__":
    run()
