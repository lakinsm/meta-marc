#!/usr/bin/env python3

"""Parse a HMMer tblout file line by line, storing relevant information. Optional output of coverage vectors for each
HMM with the skewness pipeline.
Note: Reads in the original fasta/fastq file MUST have a unique header for this to work properly.
"""


#############
## Imports ##
#############
import os.path
import argparse
import numpy as np
import mmarc_snpsearch as msnp
import sys
import matplotlib as mpl  # load mpl to set the output device, then load pyplot as plt
mpl.use('Agg')  # No X server running on cluster, use Agg for png generation instead
import matplotlib.pyplot as plt
import math



##########
## Vars ##
##########
total_reads_mapped = 0
plot_count = 0


#############
## Methods ##
#############
class TbloutParser:
    """This object takes as input a HMMer tblout file path and constructs an iterable that outputs
    hash-mapping of header to sequence information.
    """

    def __init__(self, filepath, outpath, filename, model_annot, model_set, evalue=10,
                 multi=False, dupfile=None, detect_snps=False, repfile=None, snpfasta=None, hmmerfasta=None,
                 cov_thresh=80):
        """
        Constructor:
        :param filepath: filepath to input hmmer tblout file
        :param outpath: path to output directory
        :param filename: basename for the output file (excluding extension)
        :param model_annot: filepath to the file containing the overall HMM annotations
        :param model_set: model set being used, maps (1,2, 3) -> (groupI, groupII, groupIII)
        :param evalue: Evalue threshold below which to keep a hit
        :param multi: optional flag to correct for multiple reads, maintaining a 1 to 1 ratio of reads to counts
        :param dupfile: optional filepath to the table of duplicate counts if used in previous pipeline steps
        :param detect_snps: boolean flag for whether or not to run the SNP detection pipeline
        :param repfile: file path to the SNP report output file if enabled
        :param snpfasta: file path to the SNP FASTA output file if enabled
        :param hmmerfasta: file path to the HMMER fasta hits output file if enabled
        :param cov_thresh: desired coverage threshold for model screening

        :return: name, len, start, stop, str of filepath
        """
        if os.path.exists(filepath):  # if file is a file, read from the file
            self.hmmer_file = str(filepath)
            self.stdin = False
        elif not sys.stdin.isatty():  # else read from standard in
            self.stdin = True
        else:
            raise ValueError("Parameter filepath must be a HMMer tblout file")
        self.outpath = outpath
        self.filename = filename
        self.reads_mapping = 0
        self.reads_total = 0
        self.ethreshold = float(evalue)
        self.coverage_threshold = int(cov_thresh)
        self.multicorrect = multi
        self.duplicate = dupfile
        self.screened_models = ()  # If a model doesn't meet coverage threshold, keep track of it here
        self.detect_snps = detect_snps  # Run the SNP screen for HMMER hits to groupII SNP models?
        self.snpcov_coeffs = {}  # Keep track of SNP correction coefficients
        self.vector_hash = {}  # Keep track of coverage over each mmarc model
        self.snp_observed = {}  # Observed SNP hits, necessary for use with duplicate table
        self.duplicate_table = {}  # Counts of duplicate reads from kmer screen
        self.hmm_observed = {}  # Mapping of HMMs to hits observed
        self.class_observed = {}  # Mapping of Classes to hits observed
        self.mech_observed = {}  # Mapping of Mechanisms to hits observed
        self.group_observed = {}  # Mapping of Groups to hits observed
        self.hmm_lengths = {}  # Length of each hmm, hmm # -> length
        self.gene_multihits = {}  # Will store HMM hits for each gene
        self.hmm_annots = {}  # HMM annotation mapping from key (hmm #) -> annotations
        self.hmmer_fasta = open(hmmerfasta, 'w') if hmmerfasta else None  # Optional HMMER fasta output file
        reldir = '/'.join(os.path.realpath(__file__).split('/')[:-1])
        if model_set == 1:
            self.set = 'groupI'
            self.detect_snps = False
        elif model_set == 2:
            self.set = 'groupII'
            self.snp_screen = msnp.MMSnpSearch(dist=15, fast=True, verbose=True, print_align=False,
                                               dna_mat=reldir+'/../src/matrices/NUC.4.4', aa_mat=reldir+'/../src/matrices/PAM30',
                                               search_enabled=self.detect_snps, report=repfile, fasta=snpfasta)
        else:
            self.set = 'groupIII'
            self.snp_screen = msnp.MMSnpSearch(dist=15, fast=True, verbose=True, print_align=False,
                                               dna_mat=reldir+'/../src/matrices/NUC.4.4', aa_mat=reldir+'/../src/matrices/PAM30',
                                               search_enabled=self.detect_snps, report=repfile, fasta=snpfasta)
        with open(model_annot, 'r') as hmm_annots:
            data = hmm_annots.read().split('\n')[1:]
            for line in data:
                line = line.split('\t')
                if line[0] != self.set:
                    continue
                if line:
                    self.hmm_lengths.setdefault(line[1], int(line[3]))
                    self.hmm_annots.setdefault(line[1], line[4:])
        if dupfile:
            with open(dupfile, 'r') as indup:
                ## Put each gene header into a key and its counts into the value.  Initialize the obs count dict
                data = indup.read().split('\n')
                for line in data:
                    if line:
                        temp = line.split()
                        if temp and int(temp[0]) > 1:
                            if temp[1] in self.duplicate_table:
                                print('err')
                            self.duplicate_table.setdefault(temp[1], int(temp[0]))

    def _add_observation(self, read_name, model, count):
        try:
            self.hmm_observed[model] += count
        except KeyError:
            self.hmm_observed.setdefault(model, count)
        try:
            self.gene_multihits[read_name][model][0] += count
        except KeyError:
            try:
                self.gene_multihits[read_name].setdefault(model, [count, count])
            except KeyError:
                self.gene_multihits.setdefault(read_name, {model: [count, count]})

    def _increment_vector(self, model, start, stop, count):
        if int(stop) < int(start):
            start, stop = stop, start
        try:
            self.vector_hash[model][(int(start) - 1):(int(stop))] += count  # increment affected region
        except KeyError:
            self.vector_hash[model] = np.zeros(self.hmm_lengths[model])  # new entry, initialize vector
            self.vector_hash[model][(int(start) - 1):(int(stop))] += count

    def _iterate(self):
        # Skip all leading whitespace
        while True:
            if self.stdin:
                hmmer_line = sys.stdin.readline()  # read from stdin
            else:
                hmmer_line = self.hmmer_file.readline()  # read from file
            if not hmmer_line:
                yield False # End of file
            if hmmer_line.startswith("#"):  # these lines contain header information
                continue
            else:  # all other lines are hits
                self.reads_total += 1
                if self.reads_total % 100000 == 0:  # update the counter on stdout every 100000 reads
                    sys.stdout.write("\rHits processed: {}".format(self.reads_total))
                    sys.stdout.flush()
                temp = hmmer_line.split()
                read_name = temp[0]
                if float(temp[12]) <= float(self.ethreshold):
                    self.reads_mapping += 1
                    ## Basic observation increment rules
                    hit = (read_name, temp[2], int(temp[4]), int(temp[5]), int(temp[6]), int(temp[7]), temp[12])
                    if self.set == 'groupII' or self.set == 'groupIII':
                        if not self.snp_screen.screen_hit(hit):
                            continue
                    if self.duplicate:
                        try:
                            multiplier = self.duplicate_table[read_name]
                        except KeyError:
                            multiplier = 1
                        self._add_observation(read_name, temp[2], multiplier)
                        self._increment_vector(temp[2], temp[4], temp[5], multiplier)
                    else:
                        self._add_observation(read_name, temp[2], 1)
                        self._increment_vector(temp[2], temp[4], temp[5], 1)
                    if self.hmmer_fasta:
                        self.hmmer_fasta.write('>'+'|'.join([str(x) for x in [self.reads_total, temp[2],
                            '|'.join(self.hmm_annots[temp[2]])]])+'\n'+read_name+'\n')
                    yield True
        self.hmmer_file.close()  # catch all in case this line is reached
        assert False, 'Should not reach this line'

    def correct_multihit(self):
        """
        For our purposes, reads need to maintain a one-to-one mapping with hits counted.  Therefore, we both need to
        correct for hits across multiple HMMs due to homology and also for different motifs within the same read
        hitting multiple times within the same HMM. Note that for this to work, reads must have a unique header in the
        original fastq/fasta file that is at the end of the read, separated by a pipe '|'.

        Correct the observed counts for each hierarchy and the HMMs, since we can't
        a priori know what is true and what is not.  **This correction should not be used for assembled contigs**, since
        we can't know from assembled contigs whether they should truly multi-map or not.  It is advised to only use
        this feature for raw reads.

        :return: void
        """
        if self.multicorrect:
            for key, subdict in self.gene_multihits.items():
                if subdict:
                    ## Remove models not passing screen filter
                    subdict = {k: v for k, v in subdict.items() if k not in self.screened_models}
                    ## Correct for HMM counts
                    for nkey, nvalue in subdict.items():
                        self.hmm_observed[nkey] -= nvalue[0]
                        self.hmm_observed[nkey] += float(nvalue[1]) * (float(nvalue[0]) / sum(list(zip(*subdict.values()))[0]))


    def correct_snpcov(self):
        """
        If the SNP screen was run for groupII models (for either groupII or groupIII HMMER runs), then models in the
        observed count dictionary need to be further corrected for coverage over the resistance SNPs.  So as to not
        underestimate reads, we take a proportion of the counts in the observed dictionary for that model and divide
        that observation by the % coverage over the most covered SNP.  For example, if there were 100 hits across
        SNP1 for Rifampin1 with 10 confirmed and 50 hits across SNP2 for Rifampin1 with 25 confirmed, we would take the
        larger of the two proportions (25 / 50 = 0.5) and use that as a correction multiplier for the observed counts
        for that model.  This calculation can be manually modified where indicated below to better fit your needs if
        applicable.

        :return: void
        """
        for model, snp_vals in self.snp_screen.snpcov_correction.items():
            try:
                ## Coeff is the multiplier, read documentation for details
                if model in self.screened_models:
                    continue
                avgcov = self.vector_hash[model][0]
                snpcovs = [float(x[0]) / x[1] if x[0] and x[1] else 0 for x in snp_vals]
                coeff = max([snp if (vals[1] > avgcov) else snp * (vals[1] / avgcov) for snp, vals in zip(snpcovs, snp_vals)])
                self.hmm_observed[model] *= coeff
                if self.hmm_observed[model] <= 0:
                    del self.hmm_observed[model]
                else:
                    self.snpcov_coeffs.setdefault(model, coeff)
            except KeyError:
                continue  # The model was deleted previously by coverage_screen

    def coverage_screen(self):
        """
        Empirically, from our groups numerous studies on antimicrobial surveillance, a filter on the percent of
        gene coverage over the model should be applied.  For example, a model has hits across only 20% of the gene,
        then the hits are likely due to homology instead of the gene being actually present.  This coverage screen
        is applied before any other screen, including the multicorrection; therefore this screen also affects the
        multicorrection in that if the coverage screen eliminates a model from consideration, the hits that were
        assigned to that eliminated model are then assigned to the ones or one that remains.

        For example: a read hits to Model1, Model2, and Model3.  Coverage is only above 80% in Model2, but 33 out of
        100 reads hit to all three models in regions of high homology.  Therefore, we get rid of Model1 and Model2
        from the gene_multihits table with this screen.  Additionally, no skewness metrics are reported for these hits,
        and these models are deleted from the hmm_observed table.  The result from all of this is that the deleted
        models do not contribute to the count aggregation, and our correct_multihit algorithm doesn't see them at all,
        thus resulting in those 33 multihit reads all going to Model2, so Model2 would end up with 100 hits.  If
        the coverage screen were NOT enabled in this same scenario, then Model1 = 11 hits, Model2 = 78 hits,
        and Model3 = 11 hits, since those 33 multihit reads get split amongst them evenly.

        :return: void
        """
        for model, vector in self.vector_hash.items():
            perc_cov = 100 * float(sum(vector > 0)) / len(vector)
            avg_cov = float(sum(vector)) / len(vector)
            if perc_cov < self.coverage_threshold:
                self.screened_models += (model,)
                del self.hmm_observed[model]
            else:
                self.vector_hash[model] = (perc_cov, avg_cov, np.array(vector))

    def aggregate_hierarchy(self):
        """
        Create the actual aggregated hierarchy counts based on the annotation file.  All counts are fundamentally
        aggregated from the HMM assignments.  We are interested in the accuracy of HMM assignment and classification,
        not necessarily in the absolute annotation of the reads themselves.  By default, reads are split evenly if more
        than one annotation is present in the HMM when aggregating up.  The user is left to decide whether to round the
        output values or leave them as floats.

        :return: void
        """
        for key, value in self.hmm_observed.items():
            annots = self.hmm_annots[key]
            class_annot = annots[0].split('|')
            mech_annot = annots[1].split('|')
            group_annot = annots[2].split('|')
            for entry in class_annot:
                try:
                    self.class_observed[entry] += value / float(len(class_annot))
                except KeyError:
                    self.class_observed.setdefault(entry, value / float(len(class_annot)))
            for entry in mech_annot:
                try:
                    self.mech_observed[entry] += value / float(len(mech_annot))
                except KeyError:
                    self.mech_observed.setdefault(entry, value / float(len(mech_annot)))
            for entry in group_annot:
                try:
                    self.group_observed[entry] += value / float(len(group_annot))
                except KeyError:
                    self.group_observed.setdefault(entry, value / float(len(group_annot)))

    def write_observed(self):
        """
        Output the counts (or corrected counts) to the output file based on the hierarchy of annotations.

        :return: void
        """
        pathname = self.outpath + '/' + self.filename + '.csv'
        with open(pathname, 'w') as observed_file:
            observed_file.write('Hierarchy,Name,Abundance\n')
            for key, values in self.hmm_observed.items():
                observed_file.write('HMM,' + key + ',' + str(values) + '\n')
            for key, values in self.class_observed.items():
                observed_file.write('Class,' + key + ',' + str(values) + '\n')
            for key, values in self.mech_observed.items():
                observed_file.write('Mechanism,' + key + ',' + str(values) + '\n')
            for key, values in self.group_observed.items():
                observed_file.write('Group,' + key + ',' + str(values) + '\n')

    def write_skew(self, skewfile, graphdir):
        with open(skewfile, 'w') as skf:
            skf.write(('Model,Model_Length,SNPcov_Coefficient,Hits,Percent_Coverage,Avg_Coverage_Depth,Shannon_Entropy,'
                            'Metric_Entropy,Min_Max,L2norm_Deviation,Vector\n'))
            if graphdir:
                plt.ioff()  # no interactive mode
                plt.hold(False)  # don't keep plot
                plot_count = 0
            for model, info in self.vector_hash.items():
                if model in self.screened_models or model not in self.hmm_observed.keys():
                    continue
                norm_vec = (info[2] ** 1) / sum(info[2])  # normalize so the vector sums to 1 (frequency)
                max_entropy = np.ones(len(norm_vec)) / len(norm_vec)  # this is used to calculate the maximum shannon
                                                                        # entropy for a given vector
                shannon = np.negative(sum([x * np.log2(x) for x in norm_vec if x > 0])) / np.negative(
                        sum([x * np.log2(x) for x in max_entropy]))  # % of max possible Shannon entropy
                l2norm = 1 - ((np.sqrt(sum(norm_vec * norm_vec)) * np.sqrt(len(norm_vec))) - 1) / (
                    np.sqrt(len(norm_vec)) - 1)  # % of max possible conformity to the L2 norm unit sphere
                metric_entropy = shannon / len(info[2])
                min_max = np.float64(np.min(info[2])) / np.max(info[2])
                try:
                    snpcov = self.snpcov_coeffs[model]
                except KeyError:
                    snpcov = 'NA'
                skf.write(
                    '{:s},{:d},{},{:.3f},{:.3f},{:.3f},{:.10f},{:.10f},{:.10f},{:.10f},{}\n'.format(
                        model,
                        len(info[2]),
                        snpcov,
                        self.hmm_observed[model],
                        info[0],
                        info[1],
                        shannon,
                        metric_entropy,
                        min_max,
                        l2norm,
                        ' '.join([str(x) for x in info[2].astype('int')])
                    )
                )
                if graphdir:
                    plot_count += 1
                    sys.stdout.write("\rPlots generated: {}".format(plot_count))  # update counter for user benefit
                    sys.stdout.flush()

                    ## Plot figure
                    plt.plot(info[2])
                    plt.xlabel('Nucleotide Position')
                    plt.ylabel('Observed Count')
                    plt.title('Coverage Plot for {}'.format(model))
                    plt.savefig(args.graphs + '/' + '{}.png'.format(model))
                    plt.close()  # make sure plot is closed
        print('\n')

    def __next__(self):
        if not self.stdin and type(self.hmmer_file) is str:  # only open file here if hmmer_file is a str and not fileIO
            self.hmmer_file = open(self.hmmer_file, 'r')
        for status in self._iterate():
            if not status:  # close file on EOF
                if not self.stdin:
                    self.hmmer_file.close()
                print('\nFinished parsing.')
                if self.detect_snps:
                    print('Results for SNP search:')
                    self.snp_screen.report_counts()
                print('Screening model coverage...')
                self.coverage_screen()
                if self.multicorrect:
                    print('Correcting multi-hits...')
                    self.correct_multihit()
                if self.detect_snps:
                    print('Correcting for SNP coverage...')
                    self.correct_snpcov()
                print('Aggregating counts...')
                self.aggregate_hierarchy()
                print('Writing results...')
                self.write_observed()
                if self.set == 'groupII' or self.set == 'groupIII':
                    self.snp_screen.close_all()
                if self.hmmer_fasta:
                    self.hmmer_fasta.close()
                print('Done.')
                raise StopIteration()
            else:
                return True

    def __iter__(self):
        return self


def sanitize_inputs():
    args.dupfile = None if args.dupfile in ('None', '') else args.dupfile
    args.evalue = 256 if args.dupfile in ('None', '') else args.evalue
    args.multicorrect = False if args.dupfile in ('None', '', 'False') else args.multicorrect
    args.skewfile = None if args.dupfile in ('None', '') else args.skewfile
    args.graphs = None if args.dupfile in ('None', '') else args.graphs
    args.hmmer_fasta_out = False if args.dupfile in ('None', '') else args.hmmer_fasta_out
    args.snp = False if args.dupfile in ('None', '') else args.snp
    args.snp_report = None if args.dupfile in ('None', '') else args.snp_report
    args.snp_fasta_out = None if args.dupfile in ('None', '') else args.snp_fasta_out


##############
## ArgParse ##
##############
parser = argparse.ArgumentParser('mmparse.py')
parser.add_argument('input', type=str, help='File path to HMMER tblout file, or "-" for stdin')
parser.add_argument('output', type=str, help='Base directory path for desired outputs')
parser.add_argument('filename', type=str, help='File name for this particular file (.csv format')
parser.add_argument('model_annots', type=str, help='Path to master model annotation file')
parser.add_argument('model_set', type=int, help='Model set to use (1, 2, 3)')
parser.add_argument('-c', '--coverage', type=str, default=80, help='Counts model hits above this coverage (0,100]')
parser.add_argument('-d', '--dupfile', type=str, default=None,
                    help='Path to duplicate count file if used in the kmer screen')
parser.add_argument('-e', '--evalue', type=float, default=10, help='Evalue under which to keep hits')
parser.add_argument('-m', '--multicorrect', action='store_true', default=False,
                    help='Flag: output counts should have a one to one mapping with unique inputs')
parser.add_argument('-s', '--skewfile', type=str, default=None, help='Optional output file for HMM skewness metrics')
parser.add_argument('-g', '--graphs', type=str, help='Path to output directory for graphs')
parser.add_argument('--hmmer_fasta_out', type=str, default=None, help='Path to optional output of all HMMER hits')
parser.add_argument('--snp', action='store_true', default=False,
                    help='Flag: hits to the groupII SNP models are screened for SNPs; details in documentation.')
parser.add_argument('--snp_report', type=str, default=None, help='Path to optional output for SNP search algorithm')
parser.add_argument('--snp_fasta_out', type=str, default=None, help='Path to optional output of SNP fasta reads')


##########
## Main ##
##########
if __name__ == '__main__':
    ## Parse the arguments using ArgParse
    args = parser.parse_args()
    sanitize_inputs()  # This makes the bash script slightly less cumbersome for passing arguments
    infile = args.input
    outpath = args.output
    if args.skewfile:
        if args.graphs:
            graph_dir = args.graphs
            if not os.path.isdir(graph_dir):
                os.mkdir(graph_dir)
    if not os.path.isdir(outpath):
        os.mkdir(outpath)
    mmarc_parse = TbloutParser(infile, outpath, args.filename, args.model_annots, args.model_set, args.evalue,
                             args.multicorrect, args.dupfile, args.snp, args.snp_report, args.snp_fasta_out,
                             args.hmmer_fasta_out, args.coverage)

    print('Parsing HMMER tblout:')
    for _ in mmarc_parse:
        pass

    if args.skewfile:
        mmarc_parse.write_skew(args.skewfile, args.graphs)

