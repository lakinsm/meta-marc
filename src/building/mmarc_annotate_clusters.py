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
from collections import OrderedDict


##########
## Vars ##
##########
translation = {'groupI': 0, 'groupII': 1, 'groupIII': 2, 'blast_augmented': 3}
revtranslation = {0: 'groupI', 1: 'groupII', 2: 'groupIII', 3: 'blast_augmented'}


#############
## Methods ##
#############
class HMMWalk(object):
    """
    This object walks the HMM FASTA build directory structure, storing information about each cluster encountered
    for each directory in the walk.  The intent is to maintain annotations across the model groups for each cluster
    containing the same sequence members.  In the case of version 1.00, the mappings are for groups I-III.  The
    BLAST augmented sequences are handled separately, since they belong to a unique model set.
    """

    def __init__(self, top_node, annot_file, _verbose=False):
        if not os.path.isdir(top_node):
            raise IOError('Top node is not a valid directory: {}'.format(top_node))
        self.top_node = top_node
        self.hmm_annot = {}
        self.node_list = [x for x in glob.glob(top_node+'/*') if os.path.isdir(x)]
        self.current_node = self.node_list.pop()
        self.clstr_map = OrderedDict()
        self.seq_map = OrderedDict()
        self.layer_map3 = OrderedDict()
        self.layer_map2 = OrderedDict()
        self.layer_map1 = OrderedDict()
        self.singletons = re.compile(r'singleton_addition')
        self.unique_mappings = False
        self.walk_completed = False
        self.verbose = _verbose
        self.outmap = OrderedDict()
        self.annot_map = OrderedDict()
        self.name_counts = OrderedDict()
        self._load_annots(annot_file)

    def __iter__(self):
        return self

    def _iterate(self):
        """
        For each node in layer0 (the model groups), iterate over the models and create a tree for the mappings:
            clusters -> sequences in each cluster
            sequences -> clusters for which this sequence is a member
        :return: next node in layer0
        """
        if self.current_node == 'msa/blast_augmented':
            return self.node_list.pop()
        for clstr in glob.glob(self.current_node+'/*'):
            with open(clstr, 'r') as model:
                data = model.read().split('\n')
                if sum([x.startswith('>') for x in data]) < 2:
                    continue
                headers = [x.replace('>', '') for x in data if x.startswith('>') and not self.singletons.search(x)]
                if not headers:
                    continue
                local_node = self.current_node.split('/')[-1]
                local_child = clstr.split('/')[-1].split('_')[-1]
                try:
                    self.clstr_map[local_node].setdefault(local_child, tuple(headers))
                except KeyError:
                    self.clstr_map.setdefault(local_node, OrderedDict({local_child: tuple(headers)}))
                for name in headers:
                    if name:
                        self.seq_map.setdefault(name, [[], [], []])
                        self.seq_map[name][translation[local_node]].append(local_child)
        return self.node_list.pop()

    def _load_annots(self, annot_file):
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
                self.hmm_annot.setdefault(temp[0], temp[1:4])

    def _interpret_set_result(self, parent_level, child_level, parent_name, set1_name):
        parent_node = revtranslation[parent_level]
        parent_set = self.clstr_map[parent_node][parent_name]
        set1 = self.clstr_map[revtranslation[child_level]][set1_name]
        if parent_level == 2:
            set2_name = self.layer_map3[parent_name][child_level]
            set2 = self.clstr_map[revtranslation[child_level]][set2_name]
        elif parent_level == 1:
            set2_name = self.layer_map2[parent_name]
            set2 = self.clstr_map[revtranslation[child_level]][set2_name]
        else:
            set2_name = self.layer_map1[parent_name]
            set2 = self.clstr_map[revtranslation[child_level]][set2_name]
        res = self._return_best_intersect(parent_set, set1, set2)
        return set1_name if res == 1 else set2_name

    @staticmethod
    def _return_best_intersect(parent, set1, set2):
        if len(set(parent) & set(set1)) < 2 and len(set(parent) & set(set2)) < 2:
            raise ValueError("Set comparison error")
        return 1 if len(set(parent) & set(set1)) > len(set(parent) & set(set2)) else 2

    @staticmethod
    def _fasta_parse(infile):
        with open(infile, 'r') as fastaFile:
            # Skip whitespace
            while True:
                line = fastaFile.readline()
                if line is "":
                    return  # Empty file or premature end of file?
                if line[0] is ">":
                    break
            while True:
                if line[0] is not ">":
                    raise ValueError("Records in FASTA should begin with '>'")
                header = line[1:].rstrip()
                allLines = []
                line = fastaFile.readline()
                while True:
                    if not line:
                        break
                    if line[0] is ">":
                        break
                    allLines.append(line.rstrip())
                    line = fastaFile.readline()
                yield header, "".join(allLines).replace(" ", "").replace("\r", "")
                if not line:
                    return  # Stop Iteration
            assert False, "Should not reach this line"

    def _calculate_mappings(self):
        """
        Verify that mappings are unique between each layer and build a graph for each layer.  This graph will be used
        later to ensure that models with shared sequences also share the same model name.
        1. Hash 3-tuples of clusters for which the same sequence was observed (the sequence appeared in those clusters).
        2. Create a dictionary where groupIII clusters from the hash table in step 1 are the unique keys and the values
            are 2-tuples of mappings from groupIII to groupsI-II.  There should be exactly one of these for each key.
        3. Repeat for groupII mappings to groupI.  There should be exactly one of these for each key.
        4. Reverse the direction and map from groupI to groupII.  There should be exactly one of these for each key.
        5. If all of these conditions are met, set the unique mappings flag to True.  If any condition is broken by
            a new set of mappings, then determine which of the two conflicting mappings is the most optimal by
            comparing their member sequences.
        6. Finally, name each model mapping set, breaking ties where necessary (first iteration deals only with clean
            mappings, second iteration completes the rest). This portion of the algorithm is performed by output_walk().
        :return: void
        """
        iterset = OrderedDict()
        for value in list(self.seq_map.values()):
            iterset.setdefault(tuple(x[0] if x else None for x in value), None)
        for entry in iterset.keys():
            try:
                if len([x for x in entry[:2] if x]) > len([x for x in self.layer_map3[entry[2]] if x]):
                    self.layer_map3[entry[2]] = entry[:2]
                elif len([x for x in entry[:2] if x]) == len([x for x in self.layer_map3[entry[2]] if x]):
                    self.layer_map3[entry[2]] = (self._interpret_set_result(2, 0, entry[2], entry[0]),
                                                 self._interpret_set_result(2, 1, entry[2], entry[1]))
                    if self.verbose:
                        sys.stdout.write('{} chosen as best set for layer3, node {}\n'.format(self.layer_map3[entry[2]], entry[2]))
            except KeyError:
                self.layer_map3.setdefault(entry[2], entry[:2])
            try:
                if entry[0] and self.layer_map2[entry[1]]:
                    self.layer_map2[entry[1]] = self._interpret_set_result(1, 0, entry[1], entry[0])
                    if self.verbose:
                        sys.stdout.write('{} chosen as initial set for layer2, node {}\n'.format(self.layer_map2[entry[1]], entry[0]))
            except KeyError:
                self.layer_map2.setdefault(entry[1], entry[0])
            try:
                if entry[0] and self.layer_map1[entry[0]]:
                    self.layer_map1[entry[0]] = self._interpret_set_result(0, 1, entry[0], entry[1])
                    if self.verbose:
                        sys.stdout.write('{} chosen as initial set for layer1, node {}\n'.format(self.layer_map2[entry[0]], entry[1]))
            except KeyError:
                if entry[1]:
                    self.layer_map1.setdefault(entry[0], entry[1])
        for l1, v1 in self.layer_map3.items():
            if v1[1]:
                try:
                    if self.layer_map2[v1[1]]:
                        if v1[0]:
                            self.layer_map2[v1[1]] = self._interpret_set_result(1, 0, v1[1], v1[0])
                            if self.verbose:
                                sys.stdout.write('{} chosen as best set for layer2, node {}\n'.format(self.layer_map2[v1[1]], v1[1]))
                except KeyError:
                    self.layer_map2.setdefault(v1[1], v1[0])
        for l2, v2 in self.layer_map2.items():
            if v2:
                try:
                    x = self.layer_map1[v2]
                    self.layer_map1[v2] = self._interpret_set_result(0, 1, v2, l2)
                    if self.verbose:
                        sys.stdout.write('{} chosen as best set for layer1, node {}\n'.format(self.layer_map1[v2], v2))
                except KeyError:
                    self.layer_map1.setdefault(v2, l2)
        self.unique_mappings = True

    def output_map(self):
        """
        Recreate the input directory structure with identical file content, but rename the files so that files with
        the same sequence members share the same name, and the file names describe the class-level annotation of the
        cluster.  Start with the core models and work outward for an intuitive naming convention.
        :param out_top_node: path to the top-level output directory
        :return: void
        """
        if not self.unique_mappings or not self.walk_completed:
            raise AttributeError("HMM mapping tree not unique or initial walk not completed for this object")
        # for clstr, v in self.clstr_map['groupIII'].items():
        #     mappings = [None, None, clstr]
        #     if self.layer_map3[clstr][1]:
        #         mappings[1] = self.layer_map3[clstr][1]
        #         if self.layer_map2[mappings[1]]:
        #             mappings[0] = self.layer_map2[mappings[1]]
        #     annots = set()
        #     for e, n in enumerate(mappings):
        #         if n:
        #             for seq in self.clstr_map[revtranslation[e]][n]:
        #                 annots.add(tuple(self.hmm_annot[seq]))
        #     class_annot = set(list(zip(*annots))[0])
        #     annot_string = ','.join(['|'.join(list(set(x))) for x in zip(*annots)])
        #     if len(class_annot) > 1:
        #         print(mappings)
        #         raise ValueError("Class-level model annotations must be identical for all linked clusters: {}".format(annots))
        #     try:
        #         model_name = 'mmarc_'+list(class_annot)[0]+str(self.name_counts[list(class_annot)[0]]+1)
        #         self.name_counts[list(class_annot)[0]] += 1
        #     except KeyError:
        #         self.name_counts.setdefault(list(class_annot)[0], 1)
        #         model_name = 'mmarc_'+list(class_annot)[0]+'1'
        #     for e, n in enumerate(mappings):
        #         if n:
        #             try:
        #
        #                 self.outmap[revtranslation[e]].setdefault(n, model_name)
        #                 self.annot_map[revtranslation[e]][n] = annot_string
        #             except KeyError:
        #                 self.outmap.setdefault(revtranslation[e], OrderedDict({n: model_name}))
        #                 self.annot_map.setdefault(revtranslation[e], OrderedDict({n: annot_string}))
        for clstr, v in self.clstr_map['groupII'].items():
            if clstr in self.outmap['groupII']:
                continue
            l1 = self.layer_map2[clstr]
            if l1 in self.outmap['groupI']:
                if self.verbose:
                    sys.stdout.write("Dual mapping from layer2 to layer1: {} -> {}, current: {}\n".format(clstr, l1, self.outmap['groupI'][l1]))
            mappings = [None, clstr]
            if l1:
                mappings[0] = l1
            annots = set()
            for e, n in enumerate(mappings):
                if n:
                    for seq in self.clstr_map[revtranslation[e]][n]:
                        annots.add(tuple(self.hmm_annot[seq]))
            class_annot = set(list(zip(*annots))[0])
            annot_string = ','.join(['|'.join(list(set(x))) for x in zip(*annots)])
            if len(class_annot) > 1:
                print(mappings)
                raise ValueError("Class-level model annotations must be identical for all linked clusters: {}".format(annots))
            try:
                model_name = 'mmarc_'+list(class_annot)[0]+str(self.name_counts[list(class_annot)[0]]+1)
                self.name_counts[list(class_annot)[0]] += 1
            except KeyError:
                self.name_counts.setdefault(list(class_annot)[0], 1)
                model_name = 'mmarc_'+list(class_annot)[0]+'1'
            for e, n in enumerate(mappings):
                if n:
                    try:
                        self.outmap[revtranslation[e]].setdefault(n, model_name)
                        self.annot_map[revtranslation[e]][n] = annot_string
                    except KeyError:
                        self.outmap.setdefault(revtranslation[e], OrderedDict({n: model_name}))
                        self.annot_map.setdefault(revtranslation[e], OrderedDict({n: annot_string}))
        for clstr, v in self.clstr_map['groupI'].items():
            if clstr in self.outmap['groupI']:
                continue
            annots = set()
            for seq in self.clstr_map['groupI'][clstr]:
                annots.add(tuple(self.hmm_annot[seq]))
            class_annot = set(list(zip(*annots))[0])
            annot_string = ','.join(['|'.join(list(set(x))) for x in zip(*annots)])
            if len(class_annot) > 1:
                raise ValueError("Class-level model annotations must be identical for all linked clusters: {}".format(annots))
            try:
                model_name = 'mmarc_'+list(class_annot)[0]+str(self.name_counts[list(class_annot)[0]]+1)
                self.name_counts[list(class_annot)[0]] += 1
            except:
                self.name_counts.setdefault(list(class_annot)[0], 1)
                model_name = 'mmarc_'+list(class_annot)[0]+'1'
            try:
                self.outmap['groupI'].setdefault(clstr, model_name)
                self.annot_map['groupI'][clstr] = annot_string
            except KeyError:
                self.outmap.setdefault('groupI', OrderedDict({clstr: model_name}))
                self.annot_map.setdefault('groupI', OrderedDict({clstr: annot_string}))

    def output_write(self, top_node, out_top_node):
        out_top_node = '/'.join(out_top_node.split('/'))
        if not os.path.isdir(out_top_node):
            os.mkdir(out_top_node)
        with open('mmarc_sequence_annotations.tsv', 'w') as annot_file:
            annot_file.write('Set\tName\tSize\tLength\tAnnotation\n')
            for name, v in self.clstr_map.items():
                if not os.path.isdir(out_top_node+'/'+name):
                    os.mkdir(out_top_node+'/'+name)
                for clstr, headers in v.items():
                    prefix = None
                    if name == 'groupI':
                        prefix = 'msaI_'
                    elif name == 'groupII':
                        prefix = 'msaII_'
                    elif name == 'groupIII':
                        prefix = 'msaIII_'
                    msa_data = OrderedDict({header: seq for header, seq in self._fasta_parse(top_node+'/'+name+'/'+prefix+clstr)})
                    if len([x for x in msa_data.values() if x]) < 2:
                        raise ValueError("Singleton cluster for {}".format(clstr))
                    if not all(x in msa_data for x in headers):
                        print(headers)
                        print(msa_data.keys())
                        raise ValueError("Not all heaaders present in clstr present in file")
                    lengths = set([len(x) for x in msa_data.values()])
                    if len(lengths) > 1:
                        raise ValueError("Not all sequences are the same length for cluster {}".format(clstr))
                    annot_file.write(name+'\t'+self.outmap[name][clstr].replace(' ', '_')+'\t'+str(len(msa_data.keys()))+'\t'+str(list(lengths)[0])+'\t'+self.annot_map[name][clstr]+'\n')
                    with open(out_top_node+'/'+name+'/'+self.outmap[name][clstr].replace(' ', '_')+'.fasta', 'w') as clstr_out:
                        for header, seq in msa_data.items():
                            clstr_out.write('>'+header+'\n'+seq+'\n')

    def augment_transfer(self, top_node, out_top_node):
        seen_headers = []
        with open('mmarc_sequence_annotations.tsv', 'a') as annot_file:
            for blast_file in glob.glob(top_node+'/blast_augmented/*'):
                msa_data = OrderedDict({header: seq for header, seq in self._fasta_parse(blast_file)})
                known_headers = [x for x in msa_data.keys() if x in self.hmm_annot]
                if not known_headers:
                    continue
                lengths = set([len(x) for x in msa_data.values()])
                if len(lengths) > 1:
                    raise ValueError("Not all sequences are the same length for cluster {}".format(known_headers))
                if any(x in seen_headers for x in known_headers):
                    raise ValueError("A header was seen twice in augmented data: {}".format(blast_file))
                annots = set(tuple(self.hmm_annot[x]) for x in known_headers)
                class_annot = set(list(zip(*annots))[0])
                annot_string = ','.join(['|'.join(list(set(x))) for x in zip(*annots)])
                if len(class_annot) > 1:
                    raise ValueError("Class-level model annotations must be identical for all linked clusters: {}".format(annots))
                try:
                    model_name = 'mmarc_'+list(class_annot)[0].replace(' ', '_')+str(self.name_counts[list(class_annot)[0]]+1)+'_blast_augmented'
                    self.name_counts[list(class_annot)[0]] += 1
                except KeyError:
                    self.name_counts.setdefault(list(class_annot)[0], 1)
                    model_name = 'mmarc_'+list(class_annot)[0].replace(' ', '_')+'1_blast_augmented'
                annot_file.write('groupIII\t'+model_name+'\t'+str(len(msa_data.keys()))+'\t'+str(list(lengths)[0])+'\t'+annot_string+'\n')
                with open(out_top_node+'/groupIII/'+model_name+'.fasta', 'w') as clstr_out:
                    for header, seq in msa_data.items():
                        clstr_out.write('>'+header+'\n'+seq+'\n')

    def __next__(self):
        try:
            next_node = self._iterate()
        except IndexError:
            self._calculate_mappings()
            self.walk_completed = True
            raise StopIteration()
        else:
            self.current_node = next_node


##############
## ArgParse ##
##############
parser = argparse.ArgumentParser('mmarc_annotate_clusters.py')
parser.add_argument('top_dir', type=str, help='Path to top level build FASTA directory')
parser.add_argument('database_annotations', type=str, help='File path to the top-level HMM annotation csv')
parser.add_argument('out_dir', type=str, help='Path to the desired top level output HMM MPA directory')


##########
## Main ##
##########
if __name__ == '__main__':
    args = parser.parse_args()
    hmm_builder = HMMWalk(args.top_dir, args.database_annotations, True)
    hmm_builder2 = HMMWalk(args.top_dir, args.database_annotations, False)
    while True:
        try:
            hmm_builder.__next__()
        except StopIteration:
            break
    while True:
        try:
            hmm_builder2.__next__()
        except StopIteration:
            break

    hmm_builder.output_map()
    hmm_builder2.output_map()
    if hmm_builder.output_map() != hmm_builder2.output_map():
        raise ValueError("Model construction changes from run to run; check the code")
    hmm_builder.output_write(args.top_dir, args.out_dir)
    hmm_builder.augment_transfer(args.top_dir, args.out_dir)
