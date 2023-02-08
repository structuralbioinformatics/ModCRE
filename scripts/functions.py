import os, sys, re
import gzip
import pwd
import subprocess
import time

#-------------#
# Parsers     #
#-------------#

def parse_file(file_name, gz=False):
    """
    This function parses any file and yields lines one by one.
    
    @input:
    file_name {string}
    @return:
    line {string}

    """
    if os.path.exists(file_name):
        # Initialize #
        f = None
        # Open file handle #
        if gz: f = gzip.open(file_name, "rt")
        else: f = open(file_name, "rt")
        # For each line... #
        for line in f:
            yield line.strip("\n")
        f.close()
    else:
        raise ValueError("Could not open file %s" % file_name)

def parse_fasta_file(file_name, gz=False, clean=True):
    """
    This function parses any FASTA file and yields sequences as a tuple
    of the form (identifier, sequence).

    @input:
    file_name {string}
    @return:
    line {tuple} header, sequence

    """
    # Initialize #
    identifier = ""
    sequence = ""
    # For each line... #
    for line in parse_file(file_name, gz):
        if line == "": continue
        if line[0] == ">":
            if sequence != "":
                if clean:
                    sequence = re.sub("\W|\d", "X", sequence)
                yield (identifier, sequence)
            m = re.search("^>(.+)", line)
            identifier = m.group(1)
            sequence = ""
        else:
            sequence += line.upper()
    if clean:
        sequence = re.sub("\W|\d", "X", sequence)

    yield (identifier, sequence)


def fileExist(file):
    '''
    Check existing files
    '''
    if file is not None:
        return os.path.exists(file) and os.path.isfile(file)
    else:
        return False


#-------------#
# Write       #
#-------------#

def write(file_name=None, content=None):
    """
    This function writes any {content} to a file or to stdout if no
    file is provided. If the file already exists, it pushed the {content}
    at the bottom of the file.

    @input:
    file_name {string}
    content {string}

    """
    if file_name is not None:
        try:
            f = open(file_name, "a")
            f.write("%s\n" % content)
        except:
            raise ValueError("Could not create file %s" % file_name)
    else:
        sys.stdout.write("%s\n" % content)

#-------------#
# Cluster     #
#-------------#

def submit_command_to_queue(command, queue=None, max_jobs_in_queue=None, queue_file=None, dummy_dir="/tmp", submit="qsub", qstat="qstat"):
    """
    This function submits any {command} to a cluster {queue}.

    @input:
    command {string}
    queue {string} by default it submits to any queue
    max_jobs_in_queue {int} limits the number of jobs in queue
    queue_file is a file with information specific of the cluster for running a queue

    """
    import hashlib


    if max_jobs_in_queue is not None:
      try:
        while number_of_jobs_in_queue(qstat) >= max_jobs_in_queue: time.sleep(5)
      except Exception as e:
        print "Queue error to get the number of jobs (%s). Continue submitting jobs"%e

    cwd = os.path.join(dummy_dir,"sh")
    if not os.path.exists(cwd): os.makedirs(cwd)
    script= os.path.join(cwd,"submit_"+hashlib.sha224(command).hexdigest()+".sh")
    if queue_file is not None:
      fd=open(script,"w")
      with open(queue_file,"r") as queue_standard:
        data=queue_standard.read()
        fd.write(data)
        fd.write("%s\n\n"%(command))
      fd.close()
      queue_standard.close()
      if queue is not None:
        try:
          process = subprocess.check_output([submit," -q ",queue,script])
        except Exception as e:
          print("Failed execution %s"%e)
          try: 
            os.system("%s -q %s %s" % (submit, queue,script))
          except:
            print("Failed submission")
      else:
        try:
          process = subprocess.check_output([submit,script])
        except Exception as e:
          print("Failed execution %s"%e)
          try: 
            os.system("%s %s"% (submit,script))
          except:
            print("Failed submission")

    else:
      if queue is not None:
        try:
          process = subprocess.check_output(["echo",command,"|",submit," -q ",queue])
        except Exception as e:
          print("Failed execution %s"%e) 
          try:
            os.system("echo \"%s\" | %s -q %s" % (command, submit, queue))
          except:
            print("Failed submission")
      else:
        try:
          process = subprocess.check_output(["echo",command,"|",submit])
        except Exception as e:
          print("Failed execution %s"%e) 
          try:
            os.system("echo \"%s\" | %s" % (submit,command))
          except:
            print("Failed submission")

def submit_command_to_server(command, config, dummy_dir="/tmp"):
    """
    This function submits any {command} to a cluster {queue} on a host server.

    @input:
    command {string}
    configuration data config
    dummy directory
    """
    import hashlib
    import paramiko

    #create client object
    ssh_client=paramiko.SSHClient()
    ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    #Collect data from configuration
    queue             = config.get("Cluster","cluster_queue")
    max_jobs_in_queue = config.get("Cluster","max_jobs_in_queue")
    submit            = config.get("Cluster","cluster_submit")
    qstat             = config.get("Cluster","cluster_qstat")
    host              = config.get("Cluster","server_host")
    user              = config.get("Cluster","server_user")
    passwd            = config.get("Cluster","server_passwd")
    remote_folder     = config.get("Cluster","server_directory")
    remote_outdir     = config.get("Cluster","server_output")
    scripts_path      = config.get("Paths","scripts_path")
    queue_file        = os.path.join(scripts_path,config.get("Cluster","command_queue"))

    #open connection 
    ssh_client.connect(hostname = host, username = user, password = passwd)
    
    #open a dummy folder in server
    remote_dummy_dir      = os.path.join(remote_outdir,dummy_dir)
    remote_cwd            = os.path.join(remote_dummy_dir,"sh")
    print("REMOTE SUBMISSION: remote dummy directory ..."+remote_dummy_dir)
    print("REMOTE SUBMISSION: remote command directory  ..."+remote_cwd)
    #because the address is the same
    os.system("mkdir %s"%remote_dummy_dir)
    os.system("chmod -R 777 %s"%remote_dummy_dir)
    os.system("mkdir %s"%remote_cwd)
    os.system("chmod -R 777 %s"%remote_cwd)
    #try also from the server
    stdin, stdout, stderr = ssh_client.exec_command("mkdir "+remote_dummy_dir)
    stdin, stdout, stderr = ssh_client.exec_command("chmod -R 777 "+remote_dummy_dir)
    stdin, stdout, stderr = ssh_client.exec_command("mkdir "+remote_cwd)
    stdin, stdout, stderr = ssh_client.exec_command("chmod -R 777 "+remote_cwd)
    stdin, stdout, stderr = ssh_client.exec_command("cd "+remote_cwd)

    if not os.path.exists(dummy_dir):
        os.makedirs(dummy_dir)
    os.system("chmod -R 777 "+dummy_dir)
    cwd = os.path.join(dummy_dir,"sh")
    if not os.path.exists(cwd): 
        os.makedirs(cwd)
    os.system("chmod -R 777 "+cwd)
    script= os.path.join(cwd,"submit_"+hashlib.sha224(command).hexdigest()+".sh")
    remote_script= os.path.join(remote_cwd,"submit_"+hashlib.sha224(command).hexdigest()+".sh")
    print("SERVER QUEUE FILE: "+queue_file)
    print("SERVER SUBMISSION: "+script)
    print("REMOTE SUBMISSION: "+remote_script)
    try:
      if queue_file is not None:
       #Home files
       print("REMOTE QUEUE: "+script)
       fd=open(script,"w")
       with open(queue_file,"r") as queue_standard:
        data=queue_standard.read()
        fd.write(data)
        fd.write("%s\n\n"%(command))
       fd.close()
       queue_standard.close()
       #Remote upload files
       print("REMOTE CONNECTING SSH CLIENT")
       ftp_client = ssh_client.open_sftp()
       print("REMOTE CONNECTING COPY QUEUE-SCRIPT ON CLIENT")
       ftp_client.put(script,remote_script)
       print("REMOTE CLOSE CONNECTION")
       ftp_client.close()
       #Create remote command
       if queue != "None":
          exec_command="%s -p %s %s" % (submit, queue,remote_script)
       else:
          exec_command="%s %s"% (submit,remote_script)
       print("COMMAND = "+exec_command)
      else:
       if queue != "None":
          exec_command="echo \"%s\" | %s -p %s" % (command, submit, queue)
       else:
          exec_command="echo \"%s\" | %s" % (submit,command)
       print("COMMAND = "+exec_command)

      #Run remote command
      print("SEND JOB: %s "%exec_command)
      stdin, stdout, stderr = ssh_client.exec_command(exec_command)
      print("CHECKING ERRORS OF SERVER EXECUTION")
      for line in stderr:
         error=line
         print(line)
      print("CHECKING OUTPUT OF SERVER EXECUTION")
      for line in stdout:
        print line
        if line.startswith("Submitted"): job=line.split()[-1]

      try:
        return job
      except:
        raise error

    except Exception as e:

      print("Error on server connection: "+str(e))
      return "Error on server connection: "+str(e)


def execute_in_remote(command,parameters,config,output_dir,logfile,waiting=True):    
    """
    NOTE the output_dir must exist and be the same address in local and remote host
    """
        
    remote_scripts_path = os.path.join(config.get("Cluster","server_directory"),"scripts")
    remote_python       = os.path.join(config.get("Cluster","server_python"),"python")
    shell_file=os.path.join(output_dir,os.path.basename(output_dir)+"_"+command.rstrip(".py")+".sh")
    submit_file=open(os.path.join(output_dir,os.path.basename(output_dir)+"_"+command.rstrip(".py")+".sh"),"w")
    submit_file.write("echo '#START JOB %s  ' >> %s\n"%(os.path.basename(output_dir), logfile))
    submit_file.write("chmod 777 %s \n"%logfile)
    submit_file.write("%s %s %s \n"%(remote_python,os.path.join(remote_scripts_path,command),parameters))
    submit_file.write("echo '%s  %s   DONE' >> %s\n"%(command.rstrip(".py"), os.path.basename(output_dir), logfile))
    submit_file.write("echo '#FINISHED JOB  %s' >> %s\n"%( os.path.basename(output_dir), logfile))
    submit_file.write("chmod 777 %s \n"%logfile)
    submit_file.write("chmod -R 777 %s \n"%output_dir)
    submit_file.close()
    dummy_dir = "dummy_"+command.rstrip(".py")

    keep=[]
    if os.path.exists(logfile) and os.path.isfile(logfile):
        logread=open(logfile,"r")
        for line in logread:
            data=line.split()
            if command.rstrip(".py") == data[0]  and os.path.basename(output_dir) == data[1]: continue
            if command.rstrip(".py") in line and os.path.basename(output_dir) in line: continue
            keep.append(line.strip())
        logread.close()

    if len(keep)>0:
       log=open(logfile,"w")
       for line in keep:
           log.write("%s\n"%line)
       log.close()
       os.system("chmod 777 %s"%logfile)
    else:
       if os.path.exists(logfile) and os.path.isfile(logfile): os.remove(logfile)

    iterate = True
    #timing to write in output
    n_steps=0

    print("\n\n... execute in remote .... (iterate %s , waiting %s , JOB_ID %s)"%(iterate,waiting,os.path.basename(output_dir)))
    
    #job = submit_command_to_server(" sh "+shell_file,config,dummy_dir)
    job = submit_command_to_server(" bash "+shell_file,config,dummy_dir)

    if str(job).startswith("Error"):
       return job

    print("...starting job .... (iterate %s , waiting %s , JOB_ID %s , JOB_SUBMIT %s )"%(iterate,waiting,os.path.basename(output_dir),str(job)))

    if waiting:
     print("...starting iteration.... (iterate %s , waiting %s , JOB_ID %s , JOB_SUBMIT %s )"%(iterate,waiting,os.path.basename(output_dir),str(job)))
     while(iterate):
      n_steps = n_steps + 1
      time.sleep(5)
      if n_steps>1.e+2:
        n_steps=0
        print("...waiting to finish...(JOB_ID %s , LOG %s , JOB %s)"%(os.path.basename(output_dir),logfile,str(job)))
      if os.path.exists(logfile) and os.path.isfile(logfile):
          os.system("chmod 777 %s"%logfile)
          time.sleep(2)
          if n_steps==0: print("...open log file...(JOB_ID %s , LOG %s )"%(os.path.basename(output_dir),logfile))
          try:
            logread=open(logfile,"r")
            for line in logread:
                data=line.split()
                if n_steps==0: print("...LOGFILE DATA: search command %s and JOB_ID %s"%(command.rstrip(".py"),os.path.basename(output_dir)))
                if n_steps==0: print("...LOGFILE DATA: "+line)
                if command.rstrip(".py") == data[0]  and os.path.basename(output_dir) == data[1]: 
                    print("...LOGFILE DATA FOUND: "+line)
                    iterate=False
                if command.rstrip(".py") in line and os.path.basename(output_dir) in line: 
                    print("...LOGFILE DATA FOUND: "+line)
                    iterate=False
            logread.close()
          except Exception as e:
            print("Error while reading %s: %s"%(logfile,e))
            iterate=False
            continue
    return job

def submit_command_to_server3(command, config, dummy_dir="/tmp"):
    """
    This function submits any {command} to a cluster {queue} on a host server.

    @input:
    command {string}
    configuration data config
    dummy directory
    """
    import hashlib
    import paramiko

    #create client object
    ssh_client=paramiko.SSHClient()
    ssh_client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    #Collect data from configuration
    queue             = config.get("Cluster","cluster_queue")
    max_jobs_in_queue = config.get("Cluster","max_jobs_in_queue")
    submit            = config.get("Cluster","cluster_submit")
    qstat             = config.get("Cluster","cluster_qstat")
    host              = config.get("Cluster","server_host")
    user              = config.get("Cluster","server_user")
    passwd            = config.get("Cluster","server_passwd")
    remote_folder     = config.get("Cluster","server_directory")
    remote_outdir     = config.get("Cluster","server_output")
    scripts_path      = config.get("Paths","scripts_path")
    queue_file        = os.path.join(scripts_path,config.get("Cluster","command_queue3"))

    #open connection 
    ssh_client.connect(hostname = host, username = user, password = passwd)
    print("REMOTE SUBMISSION: USING PYTHON3 ") 
    #open a dummy folder in server
    remote_dummy_dir      = os.path.join(remote_outdir,dummy_dir)
    remote_cwd            = os.path.join(remote_dummy_dir,"sh")
    print("REMOTE SUBMISSION: remote dummy directory ..."+remote_dummy_dir)
    print("REMOTE SUBMISSION: remote command directory  ..."+remote_cwd)
    #because the address is the same
    os.system("mkdir %s"%remote_dummy_dir)
    os.system("chmod -R 777 %s"%remote_dummy_dir)
    os.system("mkdir %s"%remote_cwd)
    os.system("chmod -R 777 %s"%remote_cwd)
    #try also from the server
    stdin, stdout, stderr = ssh_client.exec_command("mkdir "+remote_dummy_dir)
    stdin, stdout, stderr = ssh_client.exec_command("chmod -R 777 "+remote_dummy_dir)
    stdin, stdout, stderr = ssh_client.exec_command("mkdir "+remote_cwd)
    stdin, stdout, stderr = ssh_client.exec_command("chmod -R 777 "+remote_cwd)
    stdin, stdout, stderr = ssh_client.exec_command("cd "+remote_cwd)

    if not os.path.exists(dummy_dir):
        os.makedirs(dummy_dir)
    os.system("chmod -R 777 "+dummy_dir)
    cwd = os.path.join(dummy_dir,"sh")
    if not os.path.exists(cwd): 
        os.makedirs(cwd)
    os.system("chmod -R 777 "+cwd)
    script= os.path.join(cwd,"submit_"+hashlib.sha224(command).hexdigest()+".sh")
    remote_script= os.path.join(remote_cwd,"submit_"+hashlib.sha224(command).hexdigest()+".sh")
    print("SERVER SUBMISSION: "+script)
    print("REMOTE SUBMISSION: "+remote_script)
    try:
      if queue_file is not None:
       #Home files
       print("REMOTE QUEUE: "+script)
       fd=open(script,"w")
       with open(queue_file,"r") as queue_standard:
        data=queue_standard.read()
        fd.write(data)
        fd.write("%s\n\n"%(command))
       fd.close()
       queue_standard.close()
       #Remote upload files
       print("REMOTE CONNECTING SSH CLIENT")
       ftp_client = ssh_client.open_sftp()
       print("REMOTE CONNECTING COPY QUEUE-SCRIPT ON CLIENT")
       ftp_client.put(script,remote_script)
       print("REMOTE CLOSE CONNECTION")
       ftp_client.close()
       #Create remote command
       if queue != "None":
          exec_command="%s -p %s %s" % (submit, queue,remote_script)
       else:
          exec_command="%s %s"% (submit,remote_script)
       print("COMMAND = "+exec_command)
      else:
       if queue != "None":
          exec_command="echo \"%s\" | %s -p %s" % (command, submit, queue)
       else:
          exec_command="echo \"%s\" | %s" % (submit,command)
       print("COMMAND = "+exec_command)

      #Run remote command
      print("SEND JOB: %s "%exec_command)
      stdin, stdout, stderr = ssh_client.exec_command(exec_command)
      print("CHECKING ERRORS OF SERVER EXECUTION")
      for line in stderr:
         error=line
      #   print(line)
      print("CHECKING OUTPUT OF SERVER EXECUTION")
      for line in stdout:
        print line
        if line.startswith("Submitted"): job=line.split()[-1]

      try:
         return job
      except:
         print("ERRORS OF SERVER EXECUTION")
         raise error

    except Exception as e:

      print("Error on server connection: "+str(e))

      return "Error on server connection: "+str(e)


def execute_in_remote3(command,parameters,config,output_dir,logfile,waiting=True):    
    """
    NOTE the output_dir must exist and be the same address in local and remote host
    """
        
    remote_scripts_path = os.path.join(config.get("Cluster","server_directory"),"scripts")
    remote_python       = os.path.join(config.get("Cluster","server_python3"),"python")
    shell_file=os.path.join(output_dir,os.path.basename(output_dir)+"_"+command.rstrip(".py")+".sh")
    submit_file=open(os.path.join(output_dir,os.path.basename(output_dir)+"_"+command.rstrip(".py")+".sh"),"w")
    submit_file.write("echo '#START JOB %s  ' >> %s\n"%(os.path.basename(output_dir), logfile))
    submit_file.write("chmod 777 %s \n"%logfile)
    submit_file.write("%s %s %s \n"%(remote_python,os.path.join(remote_scripts_path,command),parameters))
    submit_file.write("echo '%s  %s   DONE' >> %s\n"%(command.rstrip(".py"), os.path.basename(output_dir), logfile))
    submit_file.write("echo '#FINISHED JOB  %s' >> %s\n"%( os.path.basename(output_dir), logfile))
    submit_file.write("chmod 777 %s \n"%logfile)
    submit_file.write("chmod -R 777 %s \n"%output_dir)
    submit_file.close()
    dummy_dir = "dummy_"+command.rstrip(".py")

    keep=[]
    if os.path.exists(logfile) and os.path.isfile(logfile):
        logread=open(logfile,"r")
        for line in logread:
            data=line.split()
            if command.rstrip(".py") == data[0]  and os.path.basename(output_dir) == data[1]: continue
            if command.rstrip(".py") in line and os.path.basename(output_dir) in line: continue
            keep.append(line.strip())
        logread.close()

    if len(keep)>0:
       log=open(logfile,"w")
       for line in keep:
           log.write("%s\n"%line)
       log.close()
       os.system("chmod 777 %s"%logfile)
    else:
       if os.path.exists(logfile) and os.path.isfile(logfile): os.remove(logfile)

    iterate = True
    #timing to write in output
    n_steps=0

    print("\n\n... execute in remote .... (iterate %s , waiting %s , JOB_ID %s)"%(iterate,waiting,os.path.basename(output_dir)))
    print("\n\n... execute IN PYTHON3 .... (iterate %s , waiting %s , JOB_ID %s)"%(iterate,waiting,os.path.basename(output_dir))) 
    #job = submit_command_to_server3(" sh "+shell_file,config,dummy_dir)
    job = submit_command_to_server3(" bash "+shell_file,config,dummy_dir)

    if str(job).startswith("Error"):
       return job

    print("...starting job .... (iterate %s , waiting %s , JOB_ID %s , JOB_SUBMIT %s )"%(iterate,waiting,os.path.basename(output_dir),str(job)))

    if waiting:
     print("...starting iteration.... (iterate %s , waiting %s , JOB_ID %s , JOB_SUBMIT %s )"%(iterate,waiting,os.path.basename(output_dir),str(job)))
     while(iterate):
      n_steps = n_steps + 1
      time.sleep(5)
      if n_steps>1.e+2:
        n_steps=0
        print("...waiting to finish...(JOB_ID %s , LOG %s , JOB %s)"%(os.path.basename(output_dir),logfile,str(job)))
      if os.path.exists(logfile) and os.path.isfile(logfile):
          os.system("chmod 777 %s"%logfile)
          time.sleep(2)
          if n_steps==0: print("...open log file...(JOB_ID %s , LOG %s )"%(os.path.basename(output_dir),logfile))
          try:
            logread=open(logfile,"r")
            for line in logread:
                data=line.split()
                if n_steps==0: print("...LOGFILE DATA: search command %s and JOB_ID %s"%(command.rstrip(".py"),os.path.basename(output_dir)))
                if n_steps==0: print("...LOGFILE DATA: "+line)
                if command.rstrip(".py") == data[0]  and os.path.basename(output_dir) == data[1]: 
                    print("...LOGFILE DATA FOUND: "+line)
                    iterate=False
                if command.rstrip(".py") in line and os.path.basename(output_dir) in line: 
                    print("...LOGFILE DATA FOUND: "+line)
                    iterate=False
            logread.close()
          except Exception as e:
            print("Error while reading %s: %s"%(logfile,e))
            iterate=False
            continue
    return job



def number_of_jobs_in_queue(qstat="qstat"):
    """
    This functions returns the number of jobs in queue for a given
    user.

    """

    # Initialize #
    user_name = get_username()

    process = subprocess.check_output([qstat, "-u", user_name])

    return len([line for line in process.split("\n") if user_name in line])

def get_username():
    """
    This functions returns the user name.

    """

    return pwd.getpwuid(os.getuid())[0]



#--------------------------------#
# Parse with configuration       #
#--------------------------------#


def parse_best_orthologs(input_file,pdb_dir,config,rank):

    orthologs_list = []

    # Get parameters from configuration
    max_orthologs = int(config.get("Parameters","max_orthologs"))
    #store families per pdb_chain
    families={}
    for line in parse_file(os.path.join(pdb_dir, "families.txt")):
            if line.startswith("#"): continue
            pdb_chain, family = line.split(";")
            families[pdb_chain] = family


    inp = open(input_file,"r")
    for line in inp:
        if line.startswith(">"): 
            #initialize
            ortholog_dict={}
            fragment_dict={}
            monomers     =[]
            dimers       =[]
            #store data
            data=line.strip().split("|")
            interval = data[1]
            start    = int(interval.split("-")[0])
            fragment_dict.setdefault("start",start)
            end      = int(interval.split("-")[1])
            fragment_dict.setdefault("end",end)
            pval     = float(data[2])
            fragment_dict.setdefault("pval",pval)
            hit_name = data[3]
            fragment_dict.setdefault("hit_name",hit_name)
            fragment_dict.setdefault("proteins",set())
            fragment_dict.setdefault("monomer",set())
            fragment_dict.setdefault("dimer",set())
            fragment_dict.setdefault("families",set())
            fragment_dict.setdefault("pdb_chains",set())
            continue
        if line.startswith("//"):
            #print "NEW ORTHOLOG LIST"
            if len(monomers)>0:
               #print monomers
               monomer_sorted = [monomer[1] for monomer in sorted(monomers,key=lambda x: x[0])]
               if not rank: max_orthologs=len(monomer_sorted)
               for thread in monomer_sorted[:min(max_orthologs,len(monomer_sorted))]:
                  uid,gene,thread_file,score,d_score = thread.split(";")
                  fragment_dict["proteins"].add((uid,gene))
                  fragment_dict["monomer"].add((thread_file,score,d_score))
                  pdb_Hchain = thread_file.split(".")[1]
                  pdb_chain  = pdb_Hchain[0:4] + "_" + pdb_Hchain[6:7]
                  fragment_dict["pdb_chains"].add(pdb_chain)
                  if families.has_key(pdb_chain): fragment_dict["families"].add(families[pdb_chain])
            if len(dimers)>0:
               #print dimers
               dimer_sorted   = [dimer[1] for dimer in sorted(dimers,key=lambda x: x[0])]
               if not rank: max_orthologs=len(dimer_sorted)
               for thread in dimer_sorted[:min(max_orthologs,len(dimer_sorted))]:
                  #read thread files 
                  monomer_A,monomer_B = thread
                  uid_A,gene_A,thread_A,score_A,d_score_A = monomer_A.split(";")
                  uid_B,gene_B,thread_B,score_B,d_score_B = monomer_B.split(";")
                  #define fragment_dict for uidA
                  fragment_dict["proteins"].add((uid_A,gene_A))
                  fragment_dict["monomer"].add((thread_A,score_A,d_score_A))
                  pdb_Hchain = thread_A.split(".")[1]
                  pdb_chain  = pdb_Hchain[0:4] + "_" + pdb_Hchain[6:7]
                  fragment_dict["pdb_chains"].add(pdb_chain)
                  if families.has_key(pdb_chain): fragment_dict["families"].add(families[pdb_chain])
                  #define fragment_dict for uidB
                  fragment_dict["proteins"].add((uid_B,gene_B))
                  fragment_dict["monomer"].add((thread_B,score_B,d_score_B))
                  pdb_Hchain = thread_B.split(".")[1]
                  pdb_chain  = pdb_Hchain[0:4] + "_" + pdb_Hchain[6:7]
                  fragment_dict["pdb_chains"].add(pdb_chain)
                  if families.has_key(pdb_chain): fragment_dict["families"].add(families[pdb_chain])
                  #define fragment_dict for dimer
                  dimer=((thread_A,score_A,d_score_A),(thread_B,score_B,d_score_B))
                  idimer=((thread_B,score_B,d_score_B),(thread_A,score_A,d_score_A))
                  if dimer in fragment_dict["dimer"]: continue
                  if idimer in fragment_dict["dimer"]: continue
                  fragment_dict["dimer"].add(dimer)
            if len(dimers)>0 or len(monomers)>0:
                  fragment_dict["dimer"]      = [x for x in fragment_dict["dimer"] ]
                  fragment_dict["monomer"]    = [x for x in fragment_dict["monomer"] ]  
                  fragment_dict["proteins"]   = [x for x in fragment_dict["proteins"] ]
                  fragment_dict["pdb_chains"] = [x for x in fragment_dict["pdb_chains"] ]
                  fragment_dict["families"]   = [x for x in fragment_dict["families"] ]
                  ortholog_dict = fragment_dict
            if ortholog_dict != {}: orthologs_list.append(ortholog_dict)
            continue
        if not line.startswith("//") and not line.startswith(">"):
            monomer_A,monomer_B = line.strip().split()
            if monomer_A == "None" or len(monomer_A.split(";"))<5: monomer_A=None
            if monomer_B == "None" or len(monomer_B.split(";"))<5: monomer_B=None
            if monomer_A is not None:
                uid_A,gene_A,thread_A,score_A,d_score_A = monomer_A.split(";")
                if score_A != 0 and d_score_A != 0 : monomers.append((float(score_A),monomer_A))
            if monomer_B is not None:
                uid_B,gene_B,thread_B,score_B,d_score_B = monomer_B.split(";")
                if score_B != 0 and d_score_B != 0 : monomers.append((float(score_B),monomer_B))
            if monomer_B is not None and monomer_A is not None:
                if score_B != 0 and d_score_B != 0 and score_A != 0 and d_score_A != 0 :
                   dimers.append((min(float(score_A),float(score_B)),(monomer_A,monomer_B)))

    inp.close()

    return orthologs_list





#-------------#
# Iteration   #
#-------------#


def done_jobs(info_file):
    done=set()
    if not fileExist(info_file): return done
    info=open(info_file,"r")
    for line in info:
        if line.startswith("#"):continue
        done.add(line.strip().split()[0])
    info.close()
    return done

def check_done(done,set_of_jobs):
    if len(set_of_jobs)<=0:return False
    done_set=set()
    for data in done:
        done_set.add(os.path.basename(data))
    testing_set=set()
    for data in set_of_jobs:
        testing_set.add(os.path.basename(data))
    if len(done)>0: return ( testing_set!= testing_set.intersection(done_set))
    else: return True


def check_submitted(submitted,pdb_files):
    submitted_set=set()
    for data in submitted:
        submitted_set.add(os.path.basename(data))
    pdb_set=set()
    for data in pdb_files:
        pdb_set.add(os.path.basename(data))
    if len(submitted)>0: return (pdb_set!= pdb_set.intersection(submitted_set))
    else: return True



