import gzip
import sys
import os
import re
import pandas as pd
from skbio.parse.sequences import parse_fastq
from IPython.display import clear_output
from itertools import izip
import numpy as np


def SLdistance(s1, s2):
    """Calculates the hamming distance between s1 and s2"""
    # Initiate np array
    matrix = np.zeros((len(s1) + 1, len(s2) + 1), dtype=np.int)
    matrix[:, 0] = np.array([i for i in xrange(len(s1) + 1)])
    matrix[0, :] = np.array([i for i in xrange(len(s2) + 1)])

    # // Classical Levenshtein part
    for i in xrange(1, len(s1) + 1):
        for j in xrange(1, len(s2) + 1):
            cost = 0
            if s1[i - 1] != s2[j - 1]:
                cost = 1
            matrix[i, j] = min(matrix[(i - 1), (j - 1)] + cost,
                               matrix[i, (j - 1)] + 1,
                               matrix[(i - 1), j] + 1)
    min_distance = matrix[len(s1)][len(s2)]

    # New Sequence-Levenshtein part

    # Truncating
    for i in xrange(0, len(s1) + 1):
        min_distance = min(min_distance, matrix[i, len(s2)])
    # Elongating
    for j in xrange(0, len(s2) + 1):
        min_distance = min(min_distance, matrix[len(s1), j])
    return min_distance
    

def split_by_index(read1, read2, barcodes, bc_pos=(26,6)):
    '''
    Splits read pairs given in `read1` and `read2` according to the list of
    barcodes given in `barcode`. 

    The position and length of the barcode can be specified in `bc_pos` as a
    (start, length) tuple.
    '''
    output_files = {}
    # Read name line MUST start with @
    fastq_tpl = '@{id}\n{seq}\n+\n{q}\n' 
    cnt = 0
    assigned = 0
    error_dict = {}
    for rec1, rec2 in izip(parse_fastq(read1), parse_fastq(read2)):
        id1, seq1, q1 = rec1
        id2, seq2, q2 = rec2
        cnt += 1
        if cnt % 1000000 == 0:
            print 'Processed\t %d records...' % cnt
            sys.stdout.flush()
        istart, ilen = bc_pos
        ind = seq1[istart:istart+ilen]
        # It's an exact match for now but we really need to accomodate
        # mismatches here. Hamming distance
        
        # I have added code that corrects 1 error
        for barcode in barcodes:
            if SLdistance(ind, barcode) == 1:
                error_dict[ind] = barcode
        
        if ind in error_dict:
            assigned += 1
            mismatch = ind
            ind = error_dict[ind]
            qstr1 = ''.join([chr(val+33) for val in q1])
            qstr2 = ''.join([chr(val+33) for val in q2])
            if not(ind in output_files):
                r1 = gzip.open('../data/%s_R1.fastq.gz' % ind, 'wb')
                r2 = gzip.open('../data/%s_R2.fastq.gz' % ind, 'wb')
                print '...created output files for {0} from the correction of {1}'.format(ind, mismatch)
                sys.stdout.flush()
                output_files[ind] = (r1, r2)
            output_files[ind][0].write(fastq_tpl.format(id=id1,seq=seq1,q=qstr1))
            output_files[ind][1].write(fastq_tpl.format(id=id2,seq=seq2,q=qstr2))
        elif ind in barcodes:
            assigned += 1
            qstr1 = ''.join([chr(val+33) for val in q1])
            qstr2 = ''.join([chr(val+33) for val in q2])
            if not(ind in output_files):
                r1 = gzip.open('../data/%s_R1.fastq.gz' % ind, 'wb')
                r2 = gzip.open('../data/%s_R2.fastq.gz' % ind, 'wb')
                print '...created output files for: %s' % ind
                sys.stdout.flush()
                output_files[ind] = (r1, r2)
            output_files[ind][0].write(fastq_tpl.format(id=id1,seq=seq1,q=qstr1))
            output_files[ind][1].write(fastq_tpl.format(id=id2,seq=seq2,q=qstr2))
    print output_files.keys()
    print 'Assigned:\t%d sequences' % assigned
    for ind,files in output_files.items():
        f1, f2 = files
        f1.close()
        f2.close()

def split_pools(barcode, dirname='../data'):
    '''
    Splits the reads in R2 file of the sample specified by `barcode` into `plus`
    and `minus` pools.
    '''
    d, _, filenames = os.walk(dirname).next()
    files = [f for f in filenames if f.startswith(barcode)]
    files_R1 = [os.path.join(dirname, f) for f in files if 'R1' in f]
    files_R2 = [os.path.join(dirname, f) for f in files if 'R2' in f]
    fastq_tpl='@{id}\n{seq}\n+\n{qual}\n'
    minus_re = re.compile('^[CT][CT][CT][AG]')
    plus_re = re.compile('^[AG][AG][AG][CT]')
    for file_R1,file_R2 in zip(files_R1, files_R2):
        print "Processing files:\t{f1}\t{f2}".format(f1=file_R1, f2=file_R2)
        sys.stdout.flush()
        cnt = 0
        cnt_plus = 0
        cnt_minus = 0
        with gzip.open(file_R1, 'rb') as gzr1, gzip.open(file_R2, 'rb') as gzr2, \
             gzip.open(os.path.join(dirname, barcode+'_minus.fastq.gz'), 'wb') as gz_minus, \
             gzip.open(os.path.join(dirname, barcode+'_plus.fastq.gz'), 'wb') as gz_plus:
            for rec1,rec2 in izip(parse_fastq(gzr1), parse_fastq(gzr2)):
                cnt += 1
                if cnt % 1000000 == 0:
                    print "\t\t{0}\trecords...".format(cnt)
                    sys.stdout.flush()
                id1, seq1, qual1 = rec1
                id2, seq2, qual2 = rec2
                qual_str = ''.join([chr(33+q) for q in qual2])
                if minus_re.match(seq1):
                    gz_minus.write(fastq_tpl.format(id=id2,seq=seq2,qual=qual_str))
                    cnt_minus += 1
                elif plus_re.match(seq1):
                    gz_plus.write(fastq_tpl.format(id=id2,seq=seq2,qual=qual_str))
                    cnt_plus += 1
            print "{0}\tplus records\t{1}\tminus records".format(cnt_plus, cnt_minus)
            sys.stdout.flush()
            

