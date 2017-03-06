#!/usr/bin/env python3

""" Hash the k-mers of a given database (AMR in our case) into a set.  Screen the reads in a fastq file against
the database hash table using each k-mer in a sliding window along the read.  If any k-mer hits the hash table,
then the read passes filter; forward and backward orientations are hashed.

Credit for the compression algorithm code blocks goes to Cathal Garvey (https://github.com/cathalgarvey/dncode.git)
"""


#############
## Imports ##
#############
import argparse
import pickle
import multiprocessing as mp
import sys
import logging
import resource
import string


##########
## Vars ##
##########
db_hash = set()  # object for the hash table
chunksize = 20000000  # limit memory consumption by reading in blocks
window = 20  # k-mer size
overall = 0  # counter for stderr writing



#############
## Methods ##
#############
def worker(chunk):
    """ This code block is used in a MapReduce manner: the chunk of code is executed many times across
    the data block (each worker receives a chunk of that block).  The reads that pass filter are written to
    the logging cache (this is because writing to stdout produces thrashing).  The logging cache is then flushed
    on every iteration of the outer loop.
    :param chunk: a chunk of reads divided amongst the pool of parallel workers
    :return: void
    """
    global db_hash
    global window
    for read_name, seq in chunk:
        for i in range(len(seq) - window + 1):
            subseq = seq[i:i + window]
            if subseq in db_hash:
                logging.info('>' + seq + '\n' + seq)
                break


def fastq_parse():
    """ Parses a fastq file in chunks of 4 lines, skipping lines that don't inform fasta creation.
    This script only accepts stdin, so use cat file | script.py for correct functionality.
    :return: generator yielding tuples of (read_name, seq)
    """
    for x in range(chunksize):
        line = sys.stdin.readline()
        if line.startswith("@"):
            read_name = line.rstrip()[1:]
            seq = sys.stdin.readline().rstrip()
            sys.stdin.readline()
            sys.stdin.readline()
        if not line:
            return  # stop iteration
        yield read_name, seq


def fasta_parse(infile):
    """ Parses a fasta file in chunks of 2 lines.
    :param infile: path to the input fasta file
    :return: generator of (header, sequence) fasta tuples
    """
    with open(infile, 'r') as fasta_file:
        # Skip whitespace
        while True:
            line = fasta_file.readline()
            if line is "":
                return  # Empty file or premature end of file?
            if line[0] is ">":
                break
        while True:
            if line[0] is not ">":
                raise ValueError("Records in FASTA should begin with '>'")
            header = line[1:].rstrip()
            all_lines = []
            line = fasta_file.readline()
            while True:
                if not line:
                    break
                if line[0] is ">":
                    break
                all_lines.append(line.rstrip())
                line = fasta_file.readline()
            yield header, "".join(all_lines).replace(" ", "").replace("\r", "")
            if not line:
                return  # Stop Iteration
        assert False, "Should not reach this line"


def split(a, n):
    """ Splits an input list into n equal chunks; this works even if modulo > 0.
    :param a: list of arbitrary length
    :param n: number of groups to split into
    :return: generator of chunks
    """
    k, m = int(len(a) / n), len(a) % n
    return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))


def current_mem_usage():
    """
    :return: current memory usage in MB
    """
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024



##############
## ArgParse ##
##############
parser = argparse.ArgumentParser('amr_skewness.py')
parser.add_argument('-d', '--database', type=str, default=None, help='File path to fasta database if new hash')
parser.add_argument('-p', '--pickle', type=str, default=None, help='Optional flag: database is a pickled hash table')
parser.add_argument('-s', '--save', type=str, default=None, help='Optional: save the hash table to a pickle file')
parser.add_argument('-n', '--num_process', type=int, default=1, help='Number of processes to run in parallel')
parser.add_argument('-k', '--kmer', type=int, default=15, help='K-mer size')


##########
## Main ##
##########
if __name__ == '__main__':
    mp.freeze_support()
    ## Parse the arguments using ArgParse
    args = parser.parse_args()
    window = args.kmer

    ## Input must be on stdin; raise error if this is not the case
    if sys.stdin.isatty():
        raise IOError('Input must be on stdin.  Use stream redirect for correct functionality: cat file | script.py')

    ## Setup the logger for output of reads to stdout.  This is necessary because writing directly to stdout
    ## in parallel causes thrashing and variable results.  The logger caches reads passed to it on every loop
    ## and flushes to stdout after the reads have been tested for membership.
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    handler = logging.StreamHandler(stream=sys.stdout)
    handler.setLevel(logging.DEBUG)
    root.addHandler(handler)

    ## Either pickle or database must be defined
    if not args.pickle and not args.database:
        raise IOError('Either --database or --pickle must be defined')

    ## Unpickle if the database hash is in a pickle file, or read in the fasta if performing a new hash
    ## Construct a hash table for the fasta database for each unique k-mer.
    ## Forward and reverse k-mers are added to the table.
    translation_map = str.maketrans('ACGT', 'TGCA')
    if args.pickle:
        db_hash = pickle.load(open(args.pickle, 'rb'))
    else:
        db_fasta = [x[1] for x in fasta_parse(args.database)]
        for seq in db_fasta:
            rev = seq[::-1]
            seqtran = seq.translate(translation_map)
            revtran = rev.translate(translation_map)
            for i in range(len(seq) - window + 1):
                temp1 = seq[i:i + window]
                temp2 = rev[i:i + window]
                temp3 = seqtran[i:i + window]
                temp4 = revtran[i:i + window]
                db_hash.add(temp1)
                db_hash.add(temp2)
                db_hash.add(temp3)
                db_hash.add(temp4)


    ## Pickle the file if save option set
    if args.save:
        pickle.dump(db_hash, open(args.save, 'wb'))

    ## Read in each fastq chunk and check for membership.  Chunk size should be set such that
    ## the block size doesn't overflow memory.  Keep in mind this block size has the potential to be doubled
    ## in the logging cache.  Chunksize is set to 20 million reads for Bovine.
    pool = mp.Pool(processes=args.num_process)  # create pool of workers for parallel processing
    while True:
        chunks = [z for z in split([x for x in fastq_parse()], args.num_process)]  # divide reads into chunks
        sys.stderr.write('\n\tMemory used: {}MB'.format(current_mem_usage()))
        check = sum([len(x) for x in chunks])  # this is the break condition for the while loop (count of reads)
        overall += check  # add to overall read count for reporting to stderr
        sys.stderr.write('\n\tTotal sequences read {}, screening...'.format(overall))
        if check is 0:
            pool.close()
            pool.join()
            pool.terminate()
            del pool
            break
        res = pool.map(worker, chunks)  # pool.map is MapReduce.  All workers must finish before proceeding.
        handler.flush()  # flush the logging cache to stdout
        sys.stderr.write('\n\tFinished block.  Loading next chunk\n')
        del chunks  # remove chunks from memory.  Otherwise memory usage will be doubled.
        if check < chunksize:
            pool.close()  # ask nicely
            pool.join()  # sigterm
            pool.terminate()  # sigkill
            del pool  # make sure pool is cleared
            break
        del check
        del res
    sys.stderr.write('\n\tTotal sequences read {}'.format(overall))


