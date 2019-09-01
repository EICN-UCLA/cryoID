#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
search_pool.py translates both the segment sequences and all the candidate protein sequences of the pool into the simplified “6-letter” code and performs a modified blastp search of the entire pool 
using the segment sequences as queries.

Author: Xiaorun Li (Lee) @ Prof. Hong Zhou lab
University of California, Los Angeles (UCLA) &
University of Science and Technology of China (USTC)

Version 1.0, released in 9/1/2019
"""

import os, operator, sys
import argparse, logging
from argparse import Namespace
import math, time, datetime
from matplotlib import pyplot as plt
from Bio import SeqIO
import search_support

# default blast format
blastformat7 = '# Fields: query acc.ver, subject acc.ver, identical, alignment length, query length, mismatches, gap opens, gaps, q. start, q. end, s. start, s. end, query seq, subject seq, evalue'

# save the results in class BlastSeq
class BlastSeq:
    
    def __init__(self, result):
        # read data from blast results
        self.data = result
        results = result.split('\t')
        self.query = results[0]
        self.candidate = results[1] # 'subject' in blast results corresponds to candidate in our case
        self.identity = int(results[2])
        self.alignlength = int(results[3])
        self.querylength = int(results[4])
        self.mismatch = int(results[5])
        self.gapopen = int(results[6])
        self.gaps = int(results[7])
        self.qstart = int(results[8])
        self.qend = int(results[9])
        self.sstart = int(results[10])
        self.send = int(results[11])
        self.qseq = results[12]
        self.sseq = results[13]
        self.evalue = float(results[14])
        #self.score = float(results[12])
        
        if len(results) == 17:
            #self.gene = results[15]
            self.length = int(results[15])
            self.composite_evalue = float(results[17])
    
    # add composite E-value     
    def set_composite_evalue(self, composite_evalue):
        self.composite_evalue = composite_evalue
        self.data +=('\t' + str(composite_evalue))
	
    # add gene names/length info	
    def set_prot_len(self, len):
        #self.gene = gene_len[0]
        self.length = int(len)
        self.data += ('\t' + str(self.length))
        
    # add overall identity % info
    def set_ave_identity(self, ave_identity):
        self.ave_identity = float(ave_identity)
        self.data += ('\t' + str(self.ave_identity))


#===========================================================================================
# read input arguments
def Parser():

    parser = argparse.ArgumentParser(description='%(prog)s translates both the segment sequences and all the candidate protein sequences of the pool into the simplified “6-letter” code \
                                     and performs a modified blast search of the entire pool using the segment sequences as queries.', 
                                     epilog='Example: %(prog)s -q queries.pdb -p candidate_pool.fasta')
    parser.add_argument('-q', '--queryfile', type=file, required=True, help='Queries (.fasta/.pdb file) to be searched')
    parser.add_argument('-p', '--candidate_pool', type=file, required=True, help='Candidate protein sequence pool file in fasta format)') 
                        # , or standard mass spetrum results (.xlsx format))')
    parser.add_argument('-o', '--outfile', default='', help='Output file basename. Default value is the same as query file basename')
    parser.add_argument('-r', '--reverse', action='store_true', default=False, help='If N/C terminal polarity of query sesquences are unknown, add -r flag')
    parser.add_argument('-l', '--length', type=int, nargs=2, default=[1, 10000], help='Set target sequence length filter for the final results. Note that this argument takes two parameters \
                        seperated by a space. Default value is 1 10000 for a sequence length range of 1~10000')
    parser.add_argument('-e', '--evalue', type=float, default=1.0, help='Expectation value (E) threshold for saving hits in blast search for each query. Default 1.0')
    parser.add_argument('-g', '--gapscore', type=int, nargs=2, default=[32767, 32767], help='Penalty score pairs for gapopen and gapextend. Note that this argument takes two parameters \
                        seperated by a space. Default value is 32767 32767 for ungapped searching. You may change it to 15 3 for gapped searching')
    parser.add_argument('-n', '--num_candidate', type=int, default=20, help='Maximum number of candidates to keep in the final result lists. Default 20')
    parser.add_argument('-v', '--verbose', action='store_true', default=False, help='Keep simplified query sequences files and all the significant blastp alignments')
    
    params = parser.parse_args()

    return params

#===========================================================================================
# check arguments if there are any errors
def Check(params):	

    # check query file format
    if not (params.queryfile.name.endswith('.fasta') or params.queryfile.name.endswith('.pdb')):
        raise TypeError('Query file '+ params.queryfile.name +' is not in .fasta/.pdb format. Exiting')
    
    # check candidate_pool file format       
    if not (params.candidate_pool.name.endswith('.fasta') or params.candidate_pool.name.endswith('.txt') or params.candidate_pool.name.endswith('.list')): # or params.candidate_pool.name.endswith('.xlsx')
        raise TypeError('Sequence pool file '+ params.candidate_pool.name +' format is not recognized. Exiting')
    
    # check the evalue            
    if params.evalue <= 0.0 or params.evalue >= 1000.0:
        raise ValueError('Your input E-value is beyond the reasonable range 0.0~1000.0. Check it again!')
    
    # check the length range     
    if not (params.length[0] > 0 and params.length[0] < params.length[1]):
        raise ValueError('Please input a reasonable length range. Check it again!')

    # check the gapscore
    if not (params.gapscore == [32767, 32767] or params.gapscore == [15, 3] or params.gapscore == [13, 3]):
        raise ValueError('Please input acceptable gapscores!(32767 32767 or 15 3 or 13 3)')
    
    if params.num_candidate <= 0:
        raise ValueError('Please input a positive interger for num_candidate!')
    
    # set output file basename    
    # if output is not specified, use the queryfile basename by default
    if params.outfile == '':
        out_base = os.path.basename(params.queryfile.name)
        out_base = os.path.splitext(out_base)[0]
        params.outfile = out_base
    #else:
        #output = params.outfile
    
    # check if output file already exists
    #if os.path.exists(output + '_alignment.txt') or os.path.exists(output + '.log'):
        #print('Output file already exists. Overwriting')
        # print output immediately
        #sys.stdout.flush()
        #raise Warning('Output file already exists. Exiting')
        
    return params


#===========================================================================================
# blast query sequences against local database
def Blast(params):
        
    #print 'Blasting the sequences...'
    #sys.stdout.flush()
    
    # save blast file name
    blast_files = []
    
    # loop over all query
    iblast = 0        
    while iblast < len(params.query_info):
        # to blast each query with SeqBlast       
        qsize = params.query_info[iblast][1]
        query_file = params.query_info[iblast][2]
        blast_file = query_file.split('_mapped.fasta')[0] + '_blast.txt'
        
        # blast each query using customized blastp search
        blastp_params = Namespace(dbsize=params.dbsize, qsize=qsize, query=query_file, outfile=blast_file, database=params.dbname, evalue=params.evalue, gapscore=params.gapscore, logfile=params.logfile)
        search_support.SeqBlast(blastp_params)
        # keep the blast result files
        blast_files.append(blast_file)
        # remove the mapped query file        
        os.remove(query_file)

        iblast += 1
    
    return blast_files


#===========================================================================================
# read blast result files and store the results
def Readfile(params):
    
    # logging setup
    logging.basicConfig(filename=params.logfile, level=logging.INFO, format='%(message)s')
    
    '''
    # first read in candidate info, including protein/gene names/length
    # dictionary to save protein, genes, length info for each candidate
    protein_info = {}
     
    # check if infofile file available
    if params.infofile and os.path.exists(params.infofile):    
        f1 = open(params.infofile, 'r')
        
        # loop over lines to get genes length info
        for i, line in enumerate(f1):
            # skip the first header line
            if line:
                row = line[:-1].split('\t')
        
            if i == 0 and row:
                if row[0] == 'Entry' and row[1] == 'Gene names' and row[2] == 'Length':
                    continue
                else:
                    print 'Warning: Infofile file %s has wrong format! Skipping it...' % params.infofile
                    logging.warn('Infofile file %s has wrong format! Skipping it...', params.infofile)
                    break
            else:
                # from left to right: protein, genes, length
                # there may be more than one sequences corresponding to current protein
                for protein in row[1].split(' '):                
                    protein_info[row[0]] = [protein, row[2]]
            
        f1.close()
    
    if not protein_info:
        print 'Warning: length info not available. Ignore the length factor of sequences and skip length filtering'
        logging.warn('Warning: length info not available. Ignore the length factor of sequences and skip length filtering')
    '''    
    
    # split files first
    blastfiles = params.blastfiles

    # list for all the quary results
    fields = []

    while blastfiles:
        
        # temp list/set to read and record the query results in blast result files
        tem = []
        f = blastfiles.pop(0)	
        f0 = open(f, 'r')					

        for i, line in enumerate(f0):
            
            # check if the file format is right
            if i == 3:
                if line[:-1] != blastformat7:
                    logging.exception('Error: ' + f + ' wrong format! Exiting')
                    raise OSError('Error: ' + f + ' wrong format! Exiting')

            # save the results in class and list
            elif line[0][0] != '#':
                field = BlastSeq(line[:-1])
                
                # add length to class BlastSeq
                if field.candidate in params.protlen:
                    field.set_prot_len(params.protlen[field.candidate])
                else:
                    field.set_prot_len('-1')
                
                tem.append(field)
                
            else:
                continue
        
        # sort the list according to candidate and evalue and then save the best result for each candidate
        tem.sort(key=operator.attrgetter('candidate', 'evalue'), reverse=False)
        
        # Filter(params) to place restrains on gap opens, align length etc.
        tem = Filter(tem, params)

        fields.extend(tem)
			
        f0.close() 
        
        # if verbosity not asked, delete the blast result file
        if not params.verbose:
            os.remove(f)
    
    
    return fields	


#=========================================================================================
# Filter(params) to place restrains on gap opens, align length etc.
def Filter(results, params):
	
    fields = []

    # save the first line for each candidate (best) and ignore the rest
    last_candidate = ''
    last_evalue = 1.0

    for field in results:
        
        # discard alignments if truncation at the ends (which are not treated as mismatches in blast) are more than 40%
        if field.querylength - (field.qend - field.qstart + 1) > 0.4 * field.querylength:
            continue
        
        # discard alignments if there are on average more than two gaps at each gap site
        #if field.gapopen > 0 and field.gaps > field.gapopen * 2:
            #continue
        
        # length filter applied here
        if field.length != -1:
            len_min, len_max = params.length
            
            if field.length < int(len_min) or field.length > int(len_max):
                continue   
        
        # if verbosity asked, keep all the aligned segments
        if params.verbose:
            fields.append(field)
        # If verbosity not asked, only keep the alignment(s) with best/lowest E-value for each query each candidate
        else:
            if field.candidate != last_candidate:
                    fields.append(field)
                    last_candidate = field.candidate
                    last_evalue = field.evalue                
            elif field.candidate == last_candidate and str(field.evalue) == str(last_evalue):
                    fields.append(field)

    return fields


#===========================================================================================
# calculate composite E-value for each candidate
def Calculate(results, params):

    # get the candidate set
    candidate = set()
    # candidate length dictionary to store the sequence length
    candidate_length = dict()    
                                     
    # get the candidates that contains segments aligning well with queries
    for field in results:
        candidate.add(field.candidate)
        # get the length for each candidates
        candidate_length[field.candidate] = field.length
        
    # candidate evalue dictionary to store the final e-value
    # initialized to the cutoff E-value
    candidate_evalue = dict.fromkeys(candidate, params.evalue)

    # candidate average identity dictionary to store the average alignment identity info
    # initialized to 0
    candidate_ave_identity = dict.fromkeys(candidate, 0.0)
    # normlized by the total length of query sequences.
    norm_factor = 0.0

    # loop over query to get the composite E-value and average identity % for each candidate
    # query_info contains query name, length, file name and number of X for all the query sequences
    for query in params.query_info:
        
        # temp evalue dictionary to get the lowest E-value from all the alignments between each candidate and current query
        # initialize the composite E-value to the cutoff E-value
        tem_evalue = dict.fromkeys(candidate, params.evalue)        
        # temp identity dictionary to get the highest identity from all the alignments between each candidate and current query
        # initialize to 0
        tem_identity = dict.fromkeys(candidate, 0.0)
        
        for field in results:
            #print field.query
            if field.query == query[0]:
                #tem_evalue[field.candidate] = min(field.evalue, tem_evalue[field.candidate])
                if field.evalue < tem_evalue[field.candidate]:
                    tem_evalue[field.candidate] = field.evalue
                    tem_identity[field.candidate] = field.identity

            elif (params.reverse) and (field.query == (query[0]+'_rev')):
                #tem_evalue[field.candidate] = min(field.evalue, tem_evalue[field.candidate])
                if field.evalue < tem_evalue[field.candidate]:
                    tem_evalue[field.candidate] = field.evalue
                    tem_identity[field.candidate] = field.identity
                
		# loop over all candidate to update the composite E-value and ave_identity         
        for icandidate in candidate_evalue:
            
            candidate_evalue[icandidate] *= tem_evalue[icandidate]

            # apply the length factor here
            if candidate_length[icandidate] > 0:
                # multiply by current minimum evalue for each candidate
                # warning: int to float!
                candidate_evalue[icandidate] *= 1.0 * candidate_length[icandidate] / 1000.0
            # if length factor not available, ignore it.    
            #else:
                #candidate_evalue[icandidate] *= 1.0
            
            # if there are segments aligned to query
            if tem_identity[icandidate] >= 1.0:
                candidate_ave_identity[icandidate] += tem_identity[icandidate]
            # if there is no aligned alignments, use half of the query length (exclude X) instead (usually blastp would fail if there is only 50% identity)
            else:
                # query contains query name, length, file name and number of X for the query
                candidate_ave_identity[icandidate] += 0.5 * (query[1] - query[3])
                
            
        norm_factor += query[1] - query[3]
    
     
    for icandidate in candidate_evalue:
        candidate_ave_identity[icandidate] /= norm_factor

    # save and sort the composite Evalue in results
    for field in results:       
        field.set_composite_evalue(candidate_evalue[field.candidate])
        field.set_ave_identity(candidate_ave_identity[field.candidate])

    results.sort(key=operator.attrgetter('composite_evalue', 'candidate'), reverse=False)
    
    return results


#=======================================================================================
# write output files
def Writefile(fields, params):

    # logging setup
    logging.basicConfig(filename=params.logfile, level=logging.INFO, format='%(message)s')
    
    # create the output file    
    out_align = params.outfile + '_alignment.txt'
    out_result = params.outfile + '_summary.txt'
    
    # write log/alignment results to file
    fout = open(out_align, 'w')
    fout1 = open(out_result, 'w')
    
    
    # save the brief summary in this list    
    brief_field = []
    last_candidate = None
    #last_query = None
    ranking = 0
    
    for i, field in enumerate(fields):
        
        # By default only keep the first 100 most promising candidates
        # may add an argument to specify a number by the user
        if ranking >= params.num_candidate:
            break
        
        # the first candidate is most promising
        if i == 0:
            #print '#############################################'
            
            # check the best candidate to see if it satisfies the criteria
            if field.ave_identity < 0.7:
                print 'Warning: No high identity (> 70%) protein candidates identified! You may want to optimize your query sequences first, otherwise the results may not be reliable.\n'
                sys.stdout.flush()
                logging.warn('Warning: No high identity (> 70%) protein candidates identified! You may want to optimize your query sequences first, otherwise the results may not be reliable.\n')
            
            out = 'Most likely candidate:\t' + field.candidate + '\nQuery set identity:\t{0:.0%}\nSequence length:\t'.format(field.ave_identity) +\
            str(field.length) + '\nComposite E-value:\t{0:1.1e}\n'.format(field.composite_evalue) + '\n'
            out += 'The candidates have been ranked by composite E-value, and the top candidates are shown in the graph in the pop-up window. Presence of a statistically significant "gap" in composite E-values between the \
top candidate and all other candidates indicates a high confidence identification. The aligned segments can be used for further verification/model building.'
                
            # print out to screen and log file
            print out
            sys.stdout.flush()
            logging.info(out)

            #print >> fout2, '\nBrief summary of the top candidates:\n'

            # print header to log/alignments files
            print >> fout1, 'Ranking, candidate, length, ave_identity %, composite_E'
            
            # check if gaps should be printed or not
            if params.gapscore == [32767, 32767]:
                print >> fout, 'Ranking, candidate, query, identity, query length, query_start, query_end, candidate_start, candidate_end, candidate_length, aligned_query, aligned_candidate, E-value, composite_E-value'
            else:
                print >> fout, 'Ranking, candidate, query, identity, query length, num_gaps, query_start, query_end, candidate_start, candidate_end, candidate_length, aligned_query, aligned_candidate, E-value, composite_E-value'
        
        
        # get the brief summary of all the promising candidates
        if field.candidate != last_candidate:
            ranking += 1
            brief_field.append(field)
            last_candidate = field.candidate
            
            print >> fout, ''
            
            # write brief summary to files.
            print >> fout1, str(ranking) + '\t' + field.candidate + '\t' + str(field.length)  + '\t' + '{0:.0%}\t{1:1.1e}'.format(field.ave_identity, field.composite_evalue)
        
        
        #if not (field.candidate == last_candidate and field.query == last.query):
        # revser mapping of aligned sequences from blast codes to degenerate groups. More readable for the users
        for k, v in search_support.aaRevMapping.iteritems():
                field.sseq = field.sseq.replace(k, v)
                field.qseq = field.qseq.replace(k, v)
        
        # write detailed alignments info to files
        if params.gapscore == [32767, 32767]:
            print >> fout, str(ranking)+'\t'+field.candidate+'\t'+field.query+'\t'+str(field.identity)+'\t'+str(field.querylength)+'\t'+str(field.qstart)+'\t'+str(field.qend)+\
            '\t'+str(field.sstart)+'\t'+str(field.send)+'\t'+str(field.length)+'\t'+field.qseq+'\t'+field.sseq+'\t'+'{0:1.1e}'.format(field.evalue)+'\t'+'{0:1.1e}'.format(field.composite_evalue)
        else:
            print >> fout, str(ranking)+'\t'+field.candidate+'\t'+field.query+'\t'+str(field.identity)+'\t'+str(field.querylength)+'\t'+str(field.gaps)+'\t'+str(field.qstart)+'\t'+str(field.qend)+\
            '\t'+str(field.sstart)+'\t'+str(field.send)+'\t'+str(field.length)+'\t'+field.qseq+'\t'+field.sseq+'\t'+'{0:1.1e}'.format(field.evalue)+'\t'+'{0:1.1e}'.format(field.composite_evalue)

    fout.close()
    fout1.close()
    #fout2.close()
    
    print '\nFor a brief summary of candidates\' information and their alignment details, please refer to %s and %s, respectively.' % (out_result, out_align)    
    print '********************************************************************'
    # print output immediately
    sys.stdout.flush()
  
    # print to log file
    logging.info('\nFor a brief summary of candidates\' information and their alignment details, please refer to %s and %s, respectively.', out_result, out_align)
    
    return brief_field


#===========================================================================================
# plot the figure of final results
def Plot(fields, params):
    
    # save the data in lists
    identity_list = []
    logE_list = []
    
    rank = 1
    
    for field in fields:
        # only show the results of the first 20 candiates at most
        if rank > 20: break
    
        # get the data for each candidates
        identity_list.append(field.ave_identity * 100.0)
        logE_list.append(math.log(field.composite_evalue, 10))
        
        rank += 1

    # plot the figures to show average identity/composite E-value of promising candidates
    # set x axis based on the data lists
    m = len(identity_list)
    xaxis = range(1, m+1)
            
    plt.figure(1)
    
    ax1 = plt.subplot(211)
    ax1.scatter(xaxis, logE_list, s=25)
    
    # show the composite E-value gap between the first candidate and the rest ones
    #line1 = array([logE_list[0] for i in xaxis])
    line1 = [logE_list[0]] * m
    ax1.plot(xaxis, line1, color='pink',  linestyle='--')
    #line2 = array([logE_list[1] for i in xaxis])
    line2 = [logE_list[1]] * m
    ax1.plot(xaxis, line2, color='pink',  linestyle='--')    
    
    ax1.set_xlabel('Ranking of candidates', fontsize=12)
    ax1.set_xlim(0.5, m + 0.5)
    ax1.set_xticks(range(1, m+1, 1))
    #ax1.set_ylim(-0, -50)
    #ax1.set_yticks(range(-0, -50, -10))
    #ax1.invert_yaxis()
    ax1.set_ylabel('Log of composite E-value', fontsize=12)       
    ax1.axis('tight')
    plt.title('Results for the first {0} candidates'.format(str(m)), fontsize=12)
    
    # plot two sub figures
    ax2 = plt.subplot(212)
    ax2.scatter(xaxis, identity_list, s=25)     
    
    # show the identity gap between the first candidate and the rest ones
    #line1 = array([identity_list[0] for i in xaxis])
    line1 = [identity_list[0]] * m
    ax2.plot(xaxis, line1, color='paleturquoise',  linestyle='--')
    
    # get the second highest identity statistics and draw a line
    identity = identity_list[1]
    for i in range(1, m):
        if identity < identity_list[i]:
            identity = identity_list[i]   
    #line2 = array([identity for i in xaxis])
    line2 = [identity ] * m
    ax2.plot(xaxis, line2, color='paleturquoise',  linestyle='--')
    
    # set axis ticks
    ax2.set_xlim(0.5, m + 0.5)
    ax2.set_xticks(range(1, m+1, 1))
    ax2.set_xlabel('Ranking of candidates', fontsize=12)
    #ax2.set_ylim(40, 100)
    #ax2.set_yticks(range(40, 110, 10))    
    ax2.set_ylabel('Query set identity (%)', fontsize=12)    
    ax2.axis('tight') 
    
    plt.subplots_adjust(hspace=0.3)

    # set figure name
    if params.outfile == '':
        fig_name = os.path.basename(params.queryfile)
        fig_name = os.path.splitext(fig_name)[0]
    else:
        fig_name = params.outfile
        
    fig_name += '_results.png'

    plt.savefig(fig_name, dpi=300)
    plt.show()
    
    '''
    try:
        subprocess.call(['xdg-open', fig_name])
    except OSError:
        print '\nCryoID failed to open the png image. You may open it yourself in your folder.'
    '''

#===========================================================================================
# process blast results and generate the final searching results
def Process(params):
    
    # read the blast results
    results = Readfile(params)
    
    # calculate compounted E-value
    results = Calculate(results, params)
    
    # write output files
    results = Writefile(results, params)
      
    return results

#===========================================================================================
# to search the queries aganist the candidate pool database
def Search(params):
    
    # save information to log file
    logfile = params.outfile + '.log'    
    logging.basicConfig(filename=logfile, level=logging.INFO, format='%(message)s')
    
    # first print the command to console/logfile
    #print '\n********************************************************************'
    print 'Starting job search_pool ' + ' '.join(sys.argv[1:]) + ' ('+str(datetime.datetime.now())[:-7] + '):\n'
    logging.info('Starting job search_pool ' + ' '.join(sys.argv[1:]) + ' ('+str(datetime.datetime.now())[:-7] + '):\n') 
    
    # start timing
    job_start = time.time()
    
	# get candidate protein sequences and their length
    # if protein name list is provided, get the sequences and their length from Uniprot (RetrieveSeq)	
    if params.candidate_pool.name.endswith('.txt') or params.candidate_pool.name.endswith('.list'):  #params.candidate_pool.name.endswith('.xlsx') or        
        retrieve_params = Namespace(candidate_pool=params.candidate_pool.name, logfile=logfile)
        seqfile_name = search_support.RetrievSeq(retrieve_params)
    # if .fasta format candidate_pool provided, use the fasta sequences file and skip length filtering & length factor in composite E-value 
    elif params.candidate_pool.name.endswith('.fasta'):
        seqfile_name = params.candidate_pool.name
        
    # get the info of sequence length
    protlen = {}
    
    for record in SeqIO.parse(seqfile_name,'fasta'):
        if len(record.id.split("|"))>1:
            protlen[record.id.split("|")[1]]=len(record.seq)
    
    
    # generate local blastp database with (degenerate) sequences and get database size
    DB_params = Namespace(seqfile=seqfile_name, logfile=logfile)
    dbsize, dbname = search_support.DBGenerate(DB_params)
    
    
    # get the mapped query sequence files for blastp search
    query_params = Namespace(queryfile=params.queryfile.name, reverse=params.reverse, verbose=params.verbose, logfile=logfile)
    query_info = search_support.QueryMap(query_params)
    # query_info list contains query name/length/file name/number of X for all the queries
    # query_info = [[query1], ... ], query1 = [name, length, file name, X_num]


    # blast query sequences against local database
    Blast_params = Namespace(query_info=query_info, dbsize=dbsize, dbname=dbname, evalue=params.evalue, gapscore=params.gapscore, logfile=logfile)
    blast_files = Blast(Blast_params)
    
    
    # interpret the blast results        
    if len(blast_files):
        # process blast results and generate the final searching results
        process_params = Namespace(blastfiles=blast_files, candidate_pool=params.candidate_pool.name, protlen=protlen, reverse=params.reverse, length=params.length, queryfile=params.queryfile.name, outfile=params.outfile, 
                                   evalue=params.evalue, verbose=params.verbose, gapscore=params.gapscore, num_candidate=params.num_candidate, query_info=query_info, logfile=logfile)
        results = Process(process_params)
        
        # print the time (in seconds) used for the job
        job_end = time.time()
        print 'Job complete in {0:.1f} seconds!'.format(job_end - job_start)
        sys.stdout.flush()
        
        # save these info to log file too
        logging.info('Job complete in {0:.1f} seconds!'.format(job_end - job_start))
        
        # if verbose plot the figure
        #if params.verbose:
        try:
            plot_params = Namespace(queryfile=params.queryfile.name, outfile=params.outfile)
            Plot(results, plot_params)
        except:
            logging.info('Failed to plot final figure.')
            print 'Failed to plot final figure.'
            sys.stdout.flush()
        
    else:
        logging.exception('No blast results found. Check your input! Exiting')        
        raise OSError('No blast results found. Check your input! Exiting')


# main function
#================================================================
if __name__ == "__main__":

    # read input arguments
    params = Parser()    

    # check if there are any errors
    params = Check(params)
    
    # search the queries aganist the candidate pool database
    Search(params)
    