#!/usr/bin/env python3

## Take a directory of multiple sequence alignment fasta files as output by USEARCH cluster_fast.  Rename the files
## Based on the sequence headers contained within each file.  This is meant to provide a more descriptive HMM name
## when the HMMs are pressed into the master files.

## Author: Steven Lakin (Steven.Lakin@colostate.edu)


#############
## Imports ##
#############
import argparse
import os
import glob
import re
import sys
from collections import Counter
from collections import OrderedDict


###############
## Variables ##
###############
msa_table = OrderedDict()
name_table = OrderedDict()
annot_table = OrderedDict()
mydir = os.getcwd()


###############
## Functions ##
###############
def load_annots(annot_file):
    """
    Parses the production database annotation file and stores the
    gene headers as keys and the annotations as values.
    :param annot_file: production database annotation file

    :return: void
    """
    with open(annot_file, 'r') as annot:
        data = annot.read().split('\n')
        for line in data:
            temp = line.split(',')
            annot_table.setdefault(temp[0], temp[1:4])


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


def hmm_walk(top_dir, out_dir):
    """
    Walks the group directories and gives each fasta file a name based
    on the sequence members.  Produces the annotation file as output as well.

    :param top_dir: directory containing the subdirectories for model MPAs
    :param annot_file: file containing sequence annotations
    :param out_dir: output directory for renamed MPAs ready for hmmbuild
    :return: void
    """
    modelfile = open(mydir+'/mmarc_model_annotations.tsv', 'w')
    memberfile = open('/'.join(mydir.split('/')[:-1])+'/mmarc_model_members.csv', 'w')
    modelfile.write('Set\tName\tSize\tLength\tAnnotation\n')
    memberfile.write('Model_Name\tMembers\n')
    for subdir in glob.glob(top_dir+'/*'):
        modelgroup = subdir.split('/')[-1]
        for file in glob.glob(subdir+'/*'):
            filename = file.split('/')[-1]
            for entry in fasta_parse(file):
                msa_table.setdefault(filename, []).append(entry)
            annots = list(zip(*[annot_table[x[0]] for x in msa_table[filename] if x[0] in annot_table]))
            if not annots:
                continue
            classes = set(annots[0])
            if len(classes) != 1:
                classes = Counter(classes).most_common()[0][0]
            else:
                classes = str(list(classes)[0])
            mechs = '|'.join(list(set(annots[1])))
            groups = '|'.join(list(set(annots[2])))
            try:
                name_table[classes] += 1
            except KeyError:
                name_table.setdefault(classes, 1)
            modelname = 'mmarc_{}{}'.format(classes.replace(' ', '_'), name_table[classes])
            if modelgroup == 'groupIII':
                modelname += '_blast_augmented'
            model_len = set([len(x[1]) for x in msa_table[filename]])
            assert len(model_len) == 1
            model_len = int(list(model_len)[0])
            memberinfo = modelname + '\t' + ','.join([x[0] for x in msa_table[filename]])
            modelfile.write('\t'.join([modelgroup, modelname, str(len(msa_table[filename])), str(model_len), classes, mechs, groups])+'\n')
            memberfile.write(memberinfo+'\n')
            with open(out_dir+'/'+modelgroup+'/'+modelname+'.fasta', 'w') as out:
                for header, seq in msa_table[filename]:
                    out.write('>{}\n{}\n'.format(header, seq))
    modelfile.close()
    memberfile.close()


if __name__ == '__main__':
    load_annots(sys.argv[2])
    hmm_walk(sys.argv[1], sys.argv[3])







