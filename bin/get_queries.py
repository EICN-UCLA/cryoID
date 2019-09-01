#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
get_queries.py identifies multiple qualified segments from the unknown cryoEM density map and semi-automatically predicts the primary sequences for the segment.

Author: Xiaorun Li (Lee) @ Prof. Hong Zhou lab
University of California, Los Angeles (UCLA) &
University of Science and Technology of China (USTC)

Version 1.0, released in 9/1/2019
"""

import sys, subprocess, argparse, logging
from progress.spinner import Spinner
import time, datetime
import signal

#================================================================
# deal with kill signal by ctrl+c or kill pid
# code modified from https://stackoverflow.com/questions/18499497/how-to-process-sigterm-signal-gracefully
class Killer:
    kill_now = False
    def __init__(self):
        signal.signal(signal.SIGINT, self.kill)
        signal.signal(signal.SIGTERM, self.kill)
        
    def kill(self, signum, frame):
        self.kill_now = True

#from argparse import Namespace

#================================================================
# read input arguments
def Parser():

    parser = argparse.ArgumentParser(description='%(prog)s identifies multiple qualified segments from the density map and semi-automatically predicts the primary sequences for the segment.', 
                                     epilog='Example: %(prog)s -m UnkwnPro.mrc -r 3.2 -s T')
    parser.add_argument('-m', '--density_map', type=file, required=True, help='Input the cryoEM density map (mrc/ccp4 format).')
    parser.add_argument('-r', '--resolution', type=float, default=3.2,
                        help='Specify the high resolution limit (Å) of the density map. One may start with the average global resolution and then fine-tune it based on the local resolution of the selected regions.')
    parser.add_argument('-s', '--symmetry', default='ANY', 
                        help='Input reconstruction symmetry of the cryoEM density map. Default value is ANY (try everything and use te highest symmetry found).')
    parser.add_argument('-n', '--num_segments', type=int, default=10, help='Specify maximum number of queries to keep during query generation. Default value is 10.')
    parser.add_argument('-b', '--build_type', default='helices_strands', const='helices_strands', nargs='?', choices=['all', 'helices_strands', 'full_model', 'trace_and_build'], 
                        help='You can choose to build with trace_and_build, just helices/strands (quick) or a full model (takes more time)')
    parser.add_argument('-t', '--trim_from_ends', type=int, default=1, help='Trim the queries by N residues from each end')    	
    parser.add_argument('-o', '--output', help='Output file basename. Default value is the same as the density map basename')
    parser.add_argument('-p', '--nproc', type=int, default=1, help='Number of processors to use')
    
    params = parser.parse_args()

    return params


#================================================================
# check parameters if there are any errors
def Check(params):	
        
    if params.resolution >= 4.5:
        raise ValueError('Resolution parameter '+ params.resolution +' is beyond the workable range (< 4.5Å)!')
        
    # output file name 
    if not params.output:
        output = params.density_map.name.split('.mrc')[0]
        params.output = output
        
    return params

        
#================================================================
# build short sequence for good regions in the density map
def QueryBuild(params):
    
    query_pdb = params.output + '_queries.pdb'
    query_fasta = params.output + '_queries.fasta'
    query_dir = params.output + '_queries'
    logfile = params.output + '.log'
    
    # first print the command
    #print '\n********************************************************************'
    print 'Starting job get_query ' + ' '.join(sys.argv[1:]) + ' (' + str(datetime.datetime.now())[:-7] + ')\n'
    sys.stdout.flush()
    
    # save these messages to log file
    logging.basicConfig(filename=logfile, level=logging.INFO, format='%(message)s')
    logging.info('Starting job get_query ' + ' '.join(sys.argv[1:]) + ' (' + str(datetime.datetime.now())[:-7] + ')\n')
            
    # start timing
    job_start = time.time()
    
    if params.num_segments <= 10:
        max_segments = 10
    else:
        max_segments = params.num_segments
        
    if params.trim_from_ends == 0:
        trim_ends = None
    else:
        trim_ends = str(params.trim_from_ends)
        
    # call phenix to build query sequence model for the density map
    querybuild_arg = ['phenix.sequence_from_map', params.density_map.name, 'resolution='+str(params.resolution), 'symmetry='+params.symmetry, 'build_type='+params.build_type, 'pdb_out='+query_pdb, 
                      'temp_dir='+query_dir, 'max_segments='+str(params.num_segments), 'maximum_segments_in_sequencing='+str(max_segments), 'minimum_length=8',
                      'multiprocessing=multiprocessing', 'quick=True', 'scattering_table=electron', 'trim_from_ends='+trim_ends, 'chain_type=PROTEIN', 'nproc='+str(params.nproc)] # 
    #extract_n=15', 'extract_length=25', 
    
    #print querybuild_arg
    
    #print 'Call phenix.sequence_from_map to generate queries...'
    print '********************************************************************'
    print 'Getting queries from the density map...'
    # print output immediately
    sys.stdout.flush()
    logging.info('Getting queries from the density map...')
    
    
    # progress bar
    spinner = Spinner('')
    # deal with kill signal. Make sure to kill the subprocess before exiting
    killer = Killer()
    
    try:
        #subprocess.check_call(querybuild_arg, stderr=subprocess.STDOUT)
        proc = subprocess.Popen(querybuild_arg, stdout=subprocess.PIPE, stderr=subprocess.PIPE) #
        
        # add progress bar here
        while True:
            time.sleep(1)
            # if received the kill signal, kill the subprocess first
            if killer.kill_now:
                proc.kill()
            
            #if subprocess not finished, show progress bar
            if proc.poll() is None:
                #out = proc.stdout.readline()
                #print out.strip()
                spinner.next()
            else:
                break
        
        # get the output and error messages
        out, error = proc.communicate()
        
    except OSError:
        logging.exception('Could not find command phenix.sequence_from_map. Make sure you include it in your PATH!')
        raise OSError('Could not find command phenix.sequence_from_map. Make sure you include it in your PATH!')
        
    # exit the progress
    if killer.kill_now:
        print 'Job was terminated by the user'
        sys.stdout.flush()
        logging.info('Job was terminated by the user')
        sys.exit()        
    
    logging.info(out)    
        
    if error or proc.returncode != 0:
        logging.exception('Error when calling phenix.sequence_from_map to get queries\n' + error)
        raise OSError('Error when calling phenix.sequence_from_map to get queries\n' + error)
    
    # call phenix to print out the built query sequences
    print_queries_arg = ['phenix.print_sequence', query_pdb, '--letter_for_mse=X']
    
    print 'The generated query sequences are:'
    # print output immediately
    sys.stdout.flush()
    logging.info('The generated query sequences are:')
    
    try:
        #subprocess.call(print_queries_arg)
        proc = subprocess.Popen(print_queries_arg, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, error = proc.communicate()       
    except OSError:
        logging.exception('Could not find command phenix.print_sequence. Make sure you include it in your PATH!')
        raise OSError('Could not find command phenix.print_sequence. Make sure you include it in your PATH!')
                
    if error:
        logging.exception('Error when calling phenix.print_sequence to convert pdb file to fasta file\n' + error)
        raise OSError('Error when calling phenix.print_sequence to convert pdb file to fasta file\n' + error)
    
    print out
    sys.stdout.flush()
    logging.info(out)
    
    # print the sequences to fasta file
    fasta_out = open(query_fasta, 'w')    
    fasta_out.write(out)           
    fasta_out.close()
    
    # print the time (in seconds) used for the job
    job_end = time.time()
    print '********************************************************************'
    print 'Job complete in {0:.1f} minutes!'.format((job_end - job_start) / 60.0)
    sys.stdout.flush()    
    logging.info('Job complete in {0:.1f} minutes!'.format((job_end - job_start) / 60.0))
    
    #flog.close()
    
    # call coot to open the density map and query pdb for usr's inspection
    coot_arg = ['coot', '--pdb', query_pdb, '--map', params.density_map.name]
    
    print 'Open the density map and query pdb for usr inspection'
    # print output immediately
    sys.stdout.flush()
    logging.info('Open the density map and query pdb for usr inspection')
    
    try:
        #subprocess.call(coot_arg)
        #proc = subprocess.Popen(coot_arg, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #out, error = proc.communicate()
        proc = subprocess.Popen(coot_arg)
    except OSError:
        #raise OSError('Could not find command Coot. Make sure you include it in your PATH!')
        logging.warn('Warning: could not find Coot in your environment! You may want to open the query model and density map in coot manually.')
        print 'Warning: could not find Coot in your environment! You may want to open the query model and density map in coot manually.'
        sys.stdout.flush()
        
    #if error:
        #raise OSError('Error when calling coot for user inspection\n' + error)


#================================================================
if __name__ == "__main__":

    # read input arguments
    params=Parser()    

    # check if there are any conflicts
    params = Check(params)
    
    # build short sequence for good regions in the density map
    QueryBuild(params)
