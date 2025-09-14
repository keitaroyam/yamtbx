# -*- coding: utf-8 -*-
"""
(c) RIKEN 2015. All rights reserved. 
Author: Keitaro Yamashita

This software is released under the new BSD License; see LICENSE.
"""
from __future__ import print_function
from __future__ import unicode_literals

import os, subprocess, re, threading, time, stat
import shlex
import shutil

# JobState
#  0: previous job is not finished.
#  1: already submitted and waiting for execution
#  2: just running
#  3: failed. next job won't start
#  4: successufully finished.
STATE_WAITING, STATE_SUBMITTED, STATE_RUNNING, STATE_FAILED, STATE_FINISHED = "waiting", "submitted", "running", "failed", "finished"


job_header = """\
#!/bin/sh
#$ -cwd
#$ -S /bin/bash

echo started at `date "+%Y-%m-%d %H:%M:%S"`
echo "host: `hostname -s` (`uname`) user: `whoami`"
echo
echo
"""

job_footer = """

echo
echo
echo finished at `date "+%Y-%m-%d %H:%M:%S"`
"""

class SgeError(Exception):
    pass
class SlurmError(Exception):
    pass

class JobManager(object): # interface
    def __init__(self): pass
    def submit(self, j): pass
    def update_stat(self, j): pass # update j's state to RUNNING/FINISHED
    def stop_all(self):pass
    def wait_all(self, jobs, interval=5, timeout=-1):
        acc = 0
        while True:
            for job in jobs: self.update_state(job)
            if all([job.state==STATE_FINISHED for job in jobs]):
                return True
            
            time.sleep(interval)
            acc += interval
            if timeout > 0 and acc > timeout: return False
    # wait_all()
# class JobManager
class LocalThread(threading.Thread):
    def __init__(self, num_jobs):
        self._stopevent = threading.Event()
        self._sleepperiod = 1.0

        self.num_jobs = num_jobs
        self.waiting_jobs = [] # [Job, ...]
        self.p_list = [] # running process list [(Job, subprocess.Popen), ..]

        threading.Thread.__init__(self)
        self.setDaemon(True)
    # __init__()

    def start_job(self, j):
        p = subprocess.Popen(os.path.join(".", j.script_name), shell=True, cwd=j.wdir,
                             stdout=open(os.path.join(j.wdir, j.script_name + ".out"), "w"),
                             stderr=open(os.path.join(j.wdir, j.script_name + ".err"), "w"),
                             universal_newlines=True)
        return p
    # start_job()
    
    def run(self):
        while not self._stopevent.isSet():
            # Change state if finished
            for j, p in self.p_list:
                if p.poll() is not None:
                    j.state = STATE_FINISHED

            # Keep unfinished jobs
            self.p_list = [p for p in self.p_list if p[1].poll() is None]

            # Register new jobs
            for i in range(self.num_jobs - len(self.p_list)):
                # Start conversion
                if len(self.waiting_jobs) > 0:
                    j = self.waiting_jobs.pop(0)
                    self.p_list.append( (j, self.start_job(j)) )
                    j.state = STATE_RUNNING
            
            time.sleep(0.1)

        if self._stopevent.isSet():
            for j, p in self.p_list:
                p.kill()
                j.state = STATE_FAILED
    # run()

    def join(self, timeout=None):
        self._stopevent.set()
        threading.Thread.join(self, timeout)
    # join()
    
# class LocalThread()
 
class ExecLocal(JobManager):
       
    def __init__(self, max_parallel):
        JobManager.__init__(self)
        self.num_jobs = max_parallel # referred by control tower when pickling
        self._thread = LocalThread(num_jobs=self.num_jobs)
        self._thread.start()
        
    # __init__()

    def submit(self, j):
        self._thread.waiting_jobs.append(j)
        j.state = STATE_SUBMITTED
    # submit()
    
    def update_state(self, j):
        # if running locally, state is changed during execution loop
        pass
                                            
    def stop_all(self):
        self._thread.join()

# class ExecLocal
        

class SGE(JobManager):
    def __init__(self, pe_name="par"):
        JobManager.__init__(self)
        self.pe_name = pe_name

        qsub_found, qstat_found = False, False
        
        for d in os.environ["PATH"].split(":"):
            if os.path.isfile(os.path.join(d, "qsub")):
                qsub_found = True
            if os.path.isfile(os.path.join(d, "qstat")):
                qstat_found = True

        if not( qsub_found and qstat_found ):
            if shutil.which("sbatch") and shutil.which("scancel") and shutil.which("squeue"):
                print("cannot find sge engine but detected slurm engine. use slurm engine instead of sge.")
                self.Slurm = Slurm(pe_name=self.pe_name)
                self.qstat = self.Slurm.qstat
                self.submit = self.Slurm.submit
                self.update_state = self.Slurm.update_state
                self.qstat = self.Slurm.qstat
                self.stop_all = self.Slurm.stop_all
                self.qdel = self.Slurm.qdel
            else:
                raise SgeError("cannot find qsub or qstat command under $PATH")
                
        self.job_id = {} # [Job: jobid]
    # __init__()

    def submit(self, j):
        ##
        # submit script
        # @return jobID 

        script_name = j.script_name
        wdir = j.wdir

        if j.nproc > 1:
            cmd = "qsub -j y -pe %s %d %s" % (self.pe_name, j.nproc, script_name)
        else:
            cmd = "qsub -j y %s" % script_name

        p = subprocess.Popen(cmd, shell=True, cwd=wdir, 
                             stdout=subprocess.PIPE, universal_newlines=True)
        p.wait()
        stdout = p.stdout.readlines()
        
        if p.returncode != 0:
            raise SgeError("qsub failed. returncode is %d.\nstdout:\n%s\n"%(p.returncode,
                                                                            stdout))
        
        r = re.search(r"^Your job ([0-9]+) ", stdout[0])
        job_id = r.group(1)
        if job_id == "":
            raise SgeError("cannot read job-id from qsub result. please contact author. stdout is:\n" % stdout)
        
        self.job_id[j] = job_id
        print("Job %s on %s is started. id=%s"%(j.script_name, j.wdir, job_id))

    # submit()

    def update_state(self, j):
        # if job_id is unknown (waiting or finished), state won't be changed
        if j in self.job_id:
            status = self.qstat(self.job_id[j])

            if status is None: # if qsub failed, flagged as FINISHED
                j.state = STATE_FINISHED
                self.job_id.pop(j)
                #j.check_after_run() # j.state may be changed to FAILED

            else: # if qsub succeeded, RUNNING or WAITING.
                j.state = STATE_RUNNING
            
    # update_state()

    def qstat(self, job_id):
        cmd = "qstat -j %s" % job_id
        p = subprocess.Popen(cmd, shell=True,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        p.wait()
        stdout = p.stdout.readlines()

        if p.returncode != 0:
            print("job %s finished (qstat returned %s)." % (job_id, p.returncode))
            return None
        #raise SgeError("qstat failed. returncode is %d.\nstdout:\n%s\n"%(p.returncode,
        #                                                                     stdout))

        status = {}

        for l in [s for s in stdout if ":" in s]:
            splitted = l.split(":")
            key = splitted[0].strip()
            val = "".join(splitted[1:]).strip()
            status[key] = val
            
        return status
    # qstat()

    def stop_all(self):
        for i in list(self.job_id.values()):
            self.qdel(i)
    # stop_all()
    
    def qdel(self, job_id):
        cmd = "qdel %s" % job_id
        p = subprocess.Popen(cmd, shell=True,
                             stdout=subprocess.PIPE, universal_newlines=True)
        p.wait()
        stdout = p.stdout.readlines()

        if p.returncode != 0:
            print("qdel %s failed."%job_id) 
            return None
        #raise SgeError("qstat failed. returncode is %d.\nstdout:\n%s\n"%(p.returncode,
        #                                                                     stdout))

    # qdel()

# class SGE

class Slurm(JobManager):
    def __init__(self, pe_name="default",mem_per_cpu="default"):
        JobManager.__init__(self)
        self.pe_name = pe_name
        self.mem_per_cpu = mem_per_cpu
        sbatch_found, squeue_found = False, False

        for d in os.environ["PATH"].split(":"):
            if os.path.isfile(os.path.join(d, "sbatch")):
                sbatch_found = True
            if os.path.isfile(os.path.join(d, "squeue")):
                squeue_found = True

        if not( sbatch_found and squeue_found ):
            raise SlurmError("cannot find sbatch or squeue command under $PATH")
                
        self.job_id = {} # [Job: jobid]
    # __init__()

    def submit(self, j):
        ##
        # submit script
        # @return jobID 

        script_name = j.script_name
        wdir = j.wdir
        partition = ""
        memory = ""
        if self.pe_name != "default" and self.pe_name != "par":
            partition = f"-p {self.pe_name}"
        if self.mem_per_cpu != "default":
            memory = f"--mem-per-cpu={self.mem_per_cpu}"
        cmd = "sbatch %s %s -c %d %s" % (partition, memory, j.nproc, script_name)

        p = subprocess.Popen(cmd, shell=True, cwd=wdir, 
                             stdout=subprocess.PIPE, universal_newlines=True)
        p.wait()
        stdout = p.stdout.readlines()
        
        if p.returncode != 0:
            raise SlurmError("sbatch failed. returncode is %d.\nstdout:\n%s\n"%(p.returncode,
                                                                            stdout))
        
        r = re.search(r"^Submitted batch job ([0-9]+)", stdout[0])
        job_id = r.group(1)
        if job_id == "":
            raise SlurmError("cannot read job-id from sbatch result. please contact author. stdout is:\n" % stdout)
        
        self.job_id[j] = job_id
        print("Job %s on %s is started. id=%s"%(j.script_name, j.wdir, job_id))
        self.submit_time = time.time()
        self.registered_job = False


    # submit()

    def update_state(self, j):
        # if job_id is unknown (waiting or finished), state won't be changed
        if j in self.job_id:
            status = self.qstat(self.job_id[j])

            if status is None: # if sbatch failed, flagged as FINISHED?
                j.state = STATE_FINISHED
                self.job_id.pop(j)
                #j.check_after_run() # j.state may be changed to FAILED

            else: # if sbatch succeeded, RUNNING or WAITING.
                j.state = STATE_RUNNING
            
    # update_state()

    def qstat(self, job_id):
        # for the problem that jobs are not displayed if squeue before the database is updated.
        if shutil.which("sacct"):
            cmd = "sacct --job %s -n --format=STATE" % job_id
            p = subprocess.Popen(cmd, shell=True,
                                 stdout=subprocess.PIPE, universal_newlines=True)
            p.wait()
            result_sacct = p.stdout.read()
            sacct_stop_state = ["COMPLETED","FAILED","CANCELED","STOPPED"]
            for state in sacct_stop_state:
                if state in result_sacct:
                    print("job %s finished (sacct showed %s)." % (job_id, state))
                    return None
            return "running"  # better to parse ST
        else:
            # Waiting period of 15 seconds on systems that do not support sacct
            if self.submit_time and ((time.time() - self.submit_time) > 15):
                cmd = "squeue --job %s" % job_id
                p = subprocess.Popen(cmd, shell=True,
                                     stdout=subprocess.PIPE, universal_newlines=True)
                p.wait()
                stdout = p.stdout.read()

                if p.returncode == 1:
                    # maybe job is finished. squeue does not retain job information for long periods.
                    return None
                elif p.returncode != 0:
                    raise SlurmError("squeue failed. returncode is %d.\nstdout:\n%s\n" % (p.returncode,
                                                                                          stdout))
                """
                example:
                     JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
                    944900       cpu  junk.sh kyamashi PD       0:00      1 (Priority)
                """

                status = {}

                if job_id not in stdout and self.registered_job:
                    print("job %s finished (squeue showed nothing)." % job_id)
                    return None
            return "running"  # better to parse ST


    # qstat()

    def stop_all(self):
        for i in list(self.job_id.values()):
            self.qdel(i)
    # stop_all()
    
    def qdel(self, job_id):
        cmd = "scancel %s" % job_id
        p = subprocess.Popen(cmd, shell=True,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        p.wait()
        stdout = p.stdout.readlines()

        if p.returncode != 0:
            print("scancel %s failed."%job_id) 
            return None

    # qdel()

# class Slurm

class Job(object):
    ##
    # This class will be overridden
    #
    def __init__(self, wdir, script_name, nproc=1, copy_environ=True):
        self.wdir = wdir
        self.state = STATE_WAITING
        self.script_name = script_name
        self.nproc = nproc
        self.expects_out = []
        self.copy_environ = copy_environ
    # __init__()

    def write_script(self, script_text):
        env = ""
        re_allowed_env = re.compile(r"^[a-zA-Z_][a-zA-Z0-9_]*$")

        if self.copy_environ:
            for k in os.environ:
                # h5lib seems to fail when invalid directory included in HDF5_PLUGIN_PATH...
                if k == "HDF5_PLUGIN_PATH":
                    sh = shlex.shlex(os.environ[k])
                    sh.whitespace=":"
                    sh.whitespace_split = True
                    env += 'export %s='%k
                    env += ":".join(['"%s"'%x for x in [x for x in sh if os.path.isdir(x)]])
                    env += "\n"
                else:
                    if re_allowed_env.match(k):
                        env += 'export %s="%s"\n' % (k, os.environ[k].replace('"', r'\"'))
                
            env += "\n"

        script = job_header + env + script_text + job_footer 
        
        open(os.path.join(self.wdir, self.script_name), "w").write(script)
        os.chmod(os.path.join(self.wdir, self.script_name) , stat.S_IXUSR + stat.S_IWUSR + stat.S_IRUSR + stat.S_IRGRP + stat.S_IROTH)
        print("job_file=", os.path.join(self.wdir, self.script_name))
    # write_script()

# class Job
class AutoJobManager(JobManager):
    def __init__(self, pe_name="par",mem_per_cpu="default"):
        engine = detect_engine()
        if engine == "sge":
            self.engine = SGE(pe_name=pe_name)
        elif engine == "slurm":
            self.engine = Slurm(pe_name=pe_name, mem_per_cpu=mem_per_cpu)
        else:
            self.engine =ExecLocal()
        self.qstat = self.engine.qstat
        self.submit = self.engine.submit
        self.update_state = self.engine.update_state
        self.qstat = self.engine.qstat
        self.stop_all = self.engine.stop_all
        self.qdel = self.engine.qdel
        
        


def detect_engine():
    if shutil.which("squeue") and shutil.which("sbatch"):
        print("slurm detected. batch.engine=slurm")
        return "slurm"
    elif shutil.which("qstat"):
        proc = subprocess.check_output(["qstat", "-help"], stderr=subprocess.PIPE)
        if " -pe " in proc:
            print("sge detected. batch.engine=sge ")
            return "sge"
        else:
            print("pbs detected. batch.engine=pbs ")
            return "pbs"
    print("job scheduler was not found. batch.engine=sh")
    return "sh"
    
