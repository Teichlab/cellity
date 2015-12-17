#!/usr/local/bin/ipython

__author__ = 'Tomislav Ilicic'
__copyright__ = 'Copyright (C) 2015 Ambrose Carr'
__license__ = 'MIT'
__version__ = 0.0


import os
import sys
import glob
import re
import ConfigParser
import subprocess

import os, fnmatch
import numpy
import pandas as pd
from operator import itemgetter
from itertools import groupby

ROOT_DIR = ""
TEMP_DIR = ""
REF_DIR = ""
TOOLS_MAPPING = ""
TOOLS_QUANTIFICATION = ""

#LIST ALL DIRECTORIES IN SELECTED DIRECTORY. NOT RECURSIVE
def list_dirs(dir):
    walk = os.walk(dir)
    
    if (len(list(os.walk(dir))) == 0):
        return []
    else : 
        return next(os.walk(dir))[1]



def print_error (message):
    print "[ERROR]: " + message;
    sys.exit(1);


#CHECK IF OBJECT HAS BEEN SPECIFIED
def check_object_specified (type, object, objects_available):

    if (len(objects_available) == 0):
        print_error("No " + type + " are available. Please create one.");
    else:
        #IF OBJECT IS NOT AVAILABLE
        if (object not in objects_available):
            counter = 0;
            objects_available_num = [];
            for key in objects_available:
                counter =  counter + 1
                objects_available_num.append(str(counter) + ". " + str(key))
            #AND SUGGEST TO USER
            object = str(object)
            print_error(type + " \"" + object + "\" not available or not specified. Create or choose from:\n" + "\n".join(objects_available_num));

def run(args):
    import subprocess
    if (not len(sys.argv) > 1):
        return False
    
    #AVAILABLE SOFTWARE AND REFERENCE GENOMES
    ROOT_DIR = ConfigSectionMap("DIRECTORIES")['root_dir']
    REF_DIR =  ConfigSectionMap("DIRECTORIES")['ref_dir']
    TEMP_DIR =  ConfigSectionMap("DIRECTORIES")['temp_dir']
    TOOLS_QUANTIFICATION = ConfigSectionMap("DIRECTORIES")['tools_quantification']
    TOOLS_MAPPING = ConfigSectionMap("DIRECTORIES")['tools_mapping']

    available_genomes = list_dirs(REF_DIR);
    available_quants = list_dirs(TOOLS_QUANTIFICATION);
    available_mapper = list_dirs(TOOLS_MAPPING);  
    commands = [] 
    commands.append("perl " + os.path.dirname(sys.argv[0]) + "/pipeline.pl")
    commands.append("-c " + args.config)
    files_process = set()


    #USER WANTS TO CREATE REFERENCE GENOME INDEX FILE
    if (args.fasta != None):
    
        #CHECK IF GTF AND GENOME FASTA FILE HAVE BEEN SET
        fasta_file =  args.fasta;
        
        if (not os.path.isfile(fasta_file)):
            print_error("FASTA sequencing file\" " + fasta_file + "\" does not exist");
        
        gtf = args.gtf
        if (gtf == None):
            print_error("GTF file (-gtf) not specified.");
        
        if (not os.path.isfile(gtf)):
            print_error("GTF file \'" + gtf + "\' does not exist.");
            
        if (args.mapping == None):
            print_error("Mapping tool (-m) not specified.");
        check_object_specified("Mapping tool (-m)", args.mapping, available_mapper);
        commands.append("-build_db")
        commands.append(fasta_file)
        commands.append("-gtf")
        commands.append(gtf)
        commands.append("-m")
        commands.append(args.mapping)
	files_process.add("1")
	log_file = REF_DIR + "/" + os.path.basename(fasta_file) + "_" + args.mapping + ".log"
     # return(subprocess.call(" ".join(commands), shell=True))
    #ALL OTHER THINGS
    else:
        #INPUT NEEDS TO BE SET AT THIS POINT
        if (args.input != None):
            #THROW ERROR IF INPUT DIR NOT EXISTEND
            input_dir = ROOT_DIR + "/" + args.input + "/raw"
            if (not os.path.exists(input_dir)):
                print_error("Input folder \'" + input_dir + "\' does not exist. Please create and fill it with raw data.")
            supported_extensions = " | ".join(ConfigSectionMap("EXTENSIONS"))
            files = glob.glob(str(input_dir) + '/*[' + supported_extensions + ' ]')
    	  #  print (str(input_dir) + '/*[' + supported_extensions + ' ]')
            if (len(files) == 0):
                print_error("No input files found in: " + input_dir)
            type = []
            for file in files:
                file_name = "."+str(file.split(os.extsep, 1)[0])
                file_name =  os.path.basename(file_name)
                split=file_name.split("#");
                no_file_name=split[len(split)-1]
                descr=no_file_name.split("_")

                if len(descr) == 2:
                    type.append(descr[1])

            types=numpy.unique(type)
            for file in files:
                cell_num = -1
                lane_num = -1 
                pair = -1
                file_name = "."+str(file.split(os.extsep, 1)[0])
                file_name =  os.path.basename(file_name)
                ext = "."+str(file.split(os.extsep, 1)[1])
                #REMOVE FILE NAME
                split=file_name.split("#");
                no_file_name=split[len(split)-1]
                descr=no_file_name.split("_")
                print no_file_name
                if len(descr) == 1:
                    cell_num = descr[0]
                elif (len(descr) == 3 ):
                    cell_num = descr[0]
                    lane_num =  descr[1]
                    pair = descr[2]
                elif len(descr) == 2:
                    cell_num = descr[0]
                    if ("1" in types and "2" in types and len(types) == 2):
                        pair = descr[1]
                    else:
                        lane_num = descr[1]

                print cell_num
                regexp = re.compile(r'\d+')
                if regexp.search(cell_num) is not None:
                    files_process.add(int(cell_num))

            commands.append("-i")
            commands.append(args.input)
            commands.append("-o")
            commands.append(args.output)
            if (args.overwrite != None):
                commands.append("-r")
            log_file = ROOT_DIR + "/" + args.input + "/" + args.input + "_" + args.output
        else:
            print_error("Input folder (-i) not specified")
        user_out=TEMP_DIR+args.input+"/"+args.output
        #PIPELINE MAIN OUTPUT DIRECTORIES
        preprocessing       = user_out+"/preprocessed";
        mapping_root        = user_out+"/mapped";
        counts_dir          = user_out+"/counts";
        #MAPPINGOUTPUT DIRECTORIES
        mapping_dir         = mapping_root+"/sam";
        sorted_dir          = mapping_root+"/sorted_bam";
        mapping_stats       = mapping_root+"/stats";
        deduplication_dir   = mapping_root+"/deduplicated_bam";

        #QUANTIFICATION OUTPUT DIRECTORIES
        standard_counts = counts_dir+"/standard";
        de_counts_dir = counts_dir+"/deduplicated";
 
        
        #MAPPING CHOOSED
        if (args.mapping != None):
            mapper = args.mapping
            #THE PROBLEM WITH TOPHAT IS THAT IT RELIES ON BOWTIE
            #HENCE IT CAN'T BE TREATED INDIVIDUALLY.
            #THEREFORE IF TOPHAT WAS SELECTED, IT LOOKS UP WHICH
            #BOWTIE VERSION IS THE NEWEST AND GETS ALL AVAILABLE
            #REFERENCE GENOME FOR BOWTIE
            #if "tophat" not in mapper:
            #mapper = get_latest_mapper(available_mapper, "bowtie2", "tophat", 0);
        
            check_object_specified("Mapping tool (-m)", mapper, available_mapper)
            
            #ONLY LIST GENOMES THAT ARE AVAIALBLE FOR A SPECIFIC MAPPING TOOL
            mapper_available_genomes = []
            for key in available_genomes:
                g = list_dirs(REF_DIR + "/" + key);
                if (mapper in g):
                    mapper_available_genomes.append(key)

            #THROW ERROR IF GENOME NOT AVAILABLE
            #if (args.genome == None):
             #   print_error("Reference genome (-g) option not specified")
            used_genome =  args.genome;
            check_object_specified("Reference genome (-g)", used_genome,mapper_available_genomes);
            commands.append("-m")
            commands.append(mapper)
            if (args.mapping_args != None):
                arg = args.mapping_args[0].split(" ")
                del args.mapping_args[0]
                arg = arg + args.mapping_args
                commands.append("-margs "+ ','.join(arg))

            commands.append("-g")
            commands.append(used_genome)
        #QUANTIFICATION CHOOSED
        if (args.quantification != None):
            check_object_specified("Quantification tool (-q)", args.quantification, available_quants)

            commands.append("-q")
            commands.append(args.quantification)
            if (args.quant_args != None):
                commands.append("-qargs "+ ','.join(args.quant_args))
            #THROW ERROR IF GENOME NOT AVAILABLE
            if (args.genome == None):
                print_error("Reference genome (-g) option not specified")
            used_genome =  args.genome;
            check_object_specified("Reference genome (-g)", used_genome,available_genomes);
            commands.append("-g")
            commands.append(used_genome)

        if (args.range  != None):
            range = args.range.split(",")
            range = map(int, range)
            files_process = set(range).intersection(files_process)
    code = execute_pipeline(files_process, commands, log_file, args.cluster)

    if (args.quantification != None and code == 0):
        output_file = standard_counts + "/" + args.input + "." + args.output
        files = find('*.counts', standard_counts)
        merge(files, output_file)
    return code
    #return True


def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result

#RETURNS A NUMPY MATRIX (GENES (ROWS) x CELLS (COLUMNS))
#SAVES RAW MERGED FILE
def merge(files, output_file):
    print("Merge counts")
    cell_name = ""
    counts = []

    #INIT EMPTY ARRAYS
    genes = numpy.empty(1)

    #SPECIFICY HEADER
    #header = numpy.I#array(["Genes"])
    files.sort()
    header = []
    #header.append("Genes")
    for i in range(0, len(files)):
        file = files[i]
        path=file.split("/")
        cell_number = path[len(path)-1]
        cell_name = path[len(path)-4] + "." + cell_number
        matrix = numpy.genfromtxt(file, dtype = 'str')
        print cell_name
        #APPEND CELL TO MATRIX AND FILENAME TO HEADER
        if matrix.size:
            counts.append(matrix[:,1])
            #header = numpy.vstack((header, cell_name))
            header.append(cell_name)
        if (i == 0):
            genes = matrix[:,0]

    num_genes = 0
    num_cells = 0

    #IF MATRIX IS NOT EMPTY
    if (len(counts) > 0):
      
        print genes.shape
        counts = numpy.asarray(counts)
        counts = counts.T
        print counts.shape
        df = pd.DataFrame(counts, index=genes, columns=header)
        df.to_csv(output_file + ".txt", index=True, header=True, sep='\t')
        return (df)
    else:
        return (0)



def execute_pipeline(files_process, commands, log_file, cluster):
    if (len(commands) == 0 or len(files_process) == 0):
        return 0
    code = -1

    files_process =  sorted(files_process)
    files_process = list(files_process)
    
    print files_process
    ranges = []

    if (len(files_process) > 1):
        for k, g in groupby(enumerate(files_process), lambda (i,x):i-x):
            group = map(itemgetter(1), g)
            ranges.append(str(group[0]) + "-" + str(group[-1]))
    
        files_process = ranges
    files_process = [str(i) for i in files_process]

    if (cluster == "ebi"):
        code = run_on_EBI(files_process, commands, log_file, cluster)
    elif (cluster == "aws"):
        code = run_on_AWS(files_process, commands, log_file, cluster)
    return code

def run_on_EBI(files_process, commands, log_file, cluster):
        import time
        code = 0
        cluster_command = []
        #cluster_command.append("bsub -n 10 -M 32000 -R \"select[gpfs] rusage[mem=32000]\"")
        cluster_command.append("bsub -R \"select[gpfs]")
        if (args.ram != None):
                cluster_command.append(" rusage[mem=" + args.ram + "]\"")
                cluster_command.append("-M " + args.ram)
        else:
                cluster_command.append("\"")
        if (args.cpu != None):
                cluster_command.append("-n " + args.cpu)
        cluster_command.append("-J \"["+ ",".join(files_process) + "]\"")
        cluster_command.append("-oo " + log_file + "_%I.log")
        cluster_command.append("\"")
        cluster_command.append(" ".join(commands))
        cluster_command.append("\"")

        #code = subprocess.call(" ".join(cluster_command), shell=True)

        #cluster_command = [];
        #cluster_command.append("bsub")
        #cluster_command.append("-J \"[1,2]\"")
        #cluster_command.append("\"sleep 1m\"")
        print " ".join(cluster_command)
        proc = subprocess.Popen(" ".join(cluster_command), stdout=subprocess.PIPE, shell=True)
        output = proc.stdout.read()
        print output
        job_id_sub=output.split()[1]
        job_id_sub=job_id_sub[1:len(job_id_sub)-1]
        print "Submitted bsub ID: " + job_id_sub
         #CHECK STATUS
        quit();
        not_done = True
        return 0;
        while (not_done):
            cluster_command = []
            cluster_command.append("bjobs "+ job_id_sub)
            proc = subprocess.Popen(" ".join(cluster_command), stdout=subprocess.PIPE, shell=True)
            output = proc.stdout.read()

            status = output.split("\n")
            if (len(status) <=1):
                not_done = False
            else:
                status = status[1:len(status)-1]
                running = []
                completed = []
                exit = []
                uknown = []
                pending = [] 
                for i in range(0,len(status)):
                    stat =  status[i]
                    des = stat.split()
                    task_id = des[len(des)-4]
                    task_id = task_id[1:len(task_id)-1]
                    current_status = des[2]
                    if (current_status == "RUN"):
                        running.append(task_id)
                    elif (current_status == "EXIT"):
                        exit.append(task_id)
                    elif (current_status == "PEND"):
                        pending.append(task_id)
                    elif (current_status == "DONE"):
                        completed.append(task_id)
                    else: 
                        uknown.append(task_id)
                files_process = [2,2]
                if (len(completed) == len(files_process) or (len(running) == 0 and len(pending) == 0 and (len(exit) + len(uknown) ==  len(files_process)))): 
                    not_done = False
                    code = 0
                else :
                    not_done = True
                    if (len(running) > 0 ):
                        print "Running: \n" + ",".join(running)
                    if (len(pending) > 0 ): 
                        print "Pending: \n" + ",".join(pending)
                    if (len(exit) > 0 ):
                        print "Exit: \n" + ",".join(exit)
                        code = 1;
                    if (len(uknown) > 0):
                        code = 1;
                        print "Uknown\n" + ",".join(uknown)
                
                time.sleep(18)
    
        return 0

def run_on_AWS(files_process, commands, log_file, cluster):
    import tempfile
    import time
    cluster_command = []
    #GENERATE COMMAND
    with tempfile.NamedTemporaryFile(delete=False) as temp:
        temp.write('#!/bin/bash\n')
        temp.write('#$ -S /bin/bash\n')
        temp.write('#$ -cwd\n')
        temp.write('#$ -j y\n')
        temp.write('#$ -N ' + os.path.basename(log_file) + "\n")
        temp.write('#$ -t ' + ",".join(files_process)+"\n")
        resources = []
        if (args.ram != None):
                resources.append('s_vmem=' + args.ram)
        if (args.cpu != None):
               resources.append('s_core=' + args.cpu)
        if (len(resources) > 0):
                temp.write('#$ -l '+ ",".join(resources))                    
                temp.write("\n")
        #print log_file + '\${JOB_ID}.log'
        temp.write('echo \"job initiated at $(date)\"\n')
        temp.write(" ".join(commands) + " \n")
        temp.write('echo \"job ended at $(date)\"\n')
        temp.flush()

        #RUN COMMAND        
        cluster_command.append("qsub")
        cluster_command.append(temp.name)
        proc = subprocess.Popen(cluster_command, stdout=subprocess.PIPE)
        output = proc.stdout.read()
        print output
        job_id_sub=output.split(" ")[2]
        job_id_sub=job_id_sub.split(".")[0]
        print "Submitted qsub ID: " + job_id_sub
    
    #CHECK STATUS
    not_done = True
    while (not_done):
        cluster_command = []
        cluster_command.append("qstat")
        proc = subprocess.Popen(cluster_command, stdout=subprocess.PIPE)
        output = proc.stdout.read()

        status = output.split("\n")
        if (len(status) <=2):
            not_done = False
        else:
            status = status[2:len(status)-1]
            running = []
            for i in range(0,len(status)):
                stat =  status[i]
                des = stat.split()
                job_id =  des[0]
                task_id = des[len(des)-1]
                if (job_id == job_id_sub):
                    running.append(task_id)
            if (len(running) == 0):
                not_done = False
            else :
                not_done = True
                print "Running task-id: \n" + ",".join(running)
                time.sleep(180)
                
    return 0

def ConfigSectionMap(section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1
    
if __name__ == "__main__":
    # process command line args
    # TODO make a switch for various datatypes using **kwargs
    import argparse
    p = argparse.ArgumentParser()
    group1 = p.add_argument_group('Data processing', 'Modules to process your RNA sequencing data')
    group1.add_argument('-i','--input', help='Input directory')
    group1.add_argument('-m','--mapping', help='Mapping algorithm')
    group1.add_argument('-margs','--mapping_args', help='Mapping algorithm additional arguments', nargs="+")
    group1.add_argument('-q','--quantification', help='Quantification algorithm')
    group1.add_argument('-q_args','--quant_args', help='Quantification algorithm additional arguments', nargs="+")
    group1.add_argument('-o','--output', help='fastq file containing R1 pairs', default="output")
    group1.add_argument('-w','--overwrite', help='Overwrite files')
    group1.add_argument('-g','--genome', help='Reference genome')
    group1.add_argument('-r', '--range', help='Index range to run')
    group2 = p.add_argument_group('Ref index', 'Reference genome index files')
    group2.add_argument('-f','--fasta', help='Fasta files')
    group2.add_argument('-gtf', help='Gene annotation file')
    
    p.add_argument('-c', '--config', help='Config file')
    p.add_argument('-ram','--ram')
    p.add_argument('-cpu','--cpu')
    p.add_argument('-l', '--cluster', help='Cluster type') 
    args = p.parse_args()
    
    Config = ConfigParser.ConfigParser()
    if (args.config == None):
        print_error("Config file not specified");
    Config.read(args.config)
    code = run(args)

