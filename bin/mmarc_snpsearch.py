"""
If a HMMER output string hits to the SNP models AND overlaps a known SNP, then:

    If the lengths of the query string and the model string are the same and fast option is enabled (default):
        1. Do local alignment around the SNP (30bp max) using Smith-Waterman
        2. Return the SW local alignment with score
    Else if the lengths do not match or the fast option is disabled:
        1. Do a global alignment on both strings using Needleman-Wunsch
        2. Return the NW global alignment with score

    1. Check walk up the aligned model string, ignore gaps, locate the SNP codon
    2. If the SNP is present in the same location on the query string, then return the location, else False
    3. If pass step 2, then score the protein context in the alignment query string,
        if there is a stop codon return False, else return the protein substitution score and protein alignment
    4. Report all hits passing this point according to user specified level of detail

Scoring matrices:
1. NUC.4.4 for DNA alignments
2. PAM30 for protein scoring

"""

import numpy as np
import mmarc_needle as mn
import os


def expand_regex(reglist):
    expanded = []
    while reglist:
        x = reglist.pop()
        if x == ']':
            multi = ''
            x = reglist.pop()
            while x != '[':
                multi += x
                x = reglist.pop()
            expanded.append(multi)
        else:
            expanded.append(x)
    return expanded[::-1]


def load_scoring_matrix(mat):
        with open(mat, 'r') as mb:
            count = 0
            submat = None
            col = None
            for line in mb:
                if line.startswith("#"):
                    continue
                elif line.startswith(" "):
                    col = line.split()
                    n = len(col)
                    submat = np.empty((n, n), dtype=np.int32)
                else:
                    row = line.split()
                    if not row:
                        continue
                    submat[:, count] = row[1:]
                    count += 1
            return submat, col


class MMSnpSearch(object):
    def __init__(self, dist, fast, verbose, print_align, dna_mat, aa_mat, search_enabled, report, fasta):
        self.dist = dist
        self.fast = fast
        self.verbose = verbose
        self.print_align = print_align
        self.search_enabled = search_enabled
        self.report = report
        self.fasta = open(fasta, 'w') if fasta else None
        self.report = open(report, 'w') if report else None
        self.snp_models = {}
        self.snpcov_correction = {}
        self.mbmat, self.mbcol = load_scoring_matrix(aa_mat)
        self.dnamat, self.dnacol = load_scoring_matrix(dna_mat)
        reldir = '/'.join(os.path.realpath(__file__).split('/')[:-1])
        self.load_snp_metadata(reldir+'/../src/mmarc_snpsearch_metadata2.txt')
        self.snp_modelhit = 0
        self.snp_overlap = 0
        self.snp_found = 0
        self.dna2aa = {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': 'Stop', 'TAG': 'Stop',
            'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'TGT': 'C', 'TGC': 'C', 'TGA': 'Stop', 'TGG': 'W',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }
        self.ambig = {
            'U': 'T', 'M': 'AC', 'R': 'AG', 'W': 'AT', 'S': 'CG',
            'Y': 'CT', 'K': 'GT', 'V': 'ACG', 'H': 'ACT',
            'D': 'AGT', 'B': 'CGT', 'N': 'ACGT'
        }

    def load_snp_metadata(self, meta_file):
        with open(meta_file, 'r') as meta:
            data = meta.read().strip().split('\n')
            for line in data[1:]:
                temp = line.split(',')
                if temp[1] == '1':
                    try:
                        self.snp_models[temp[0]].append(temp[2:])
                    except KeyError:
                        self.snp_models.setdefault(temp[0], [temp[2:]])

    def expand_ambig(self, ambig_string):
        expanded = ['']
        for x in ambig_string:
            if x in 'ACGT-':
                expanded = [y + x for y in expanded]
            else:
                expanded *= len(self.ambig[x])
                for e, i in enumerate(self.ambig[x]):
                    for z in range(0, len(expanded), len(self.ambig[x])):
                        expanded[z + e] += i
        return expanded

    def codon_translate(self, codon):
        if any(x == '-' for x in codon):
            return '*'
        else:
            return self.dna2aa[codon]

    def get_mbcost(self, a1, a2):
        if len(a1) > 1:
            return max([self.mbmat[self.mbcol.index(x), self.mbcol.index(a2)] for x in a1])
        else:
            return self.mbmat[self.mbcol.index(a1), self.mbcol.index(a2)]

    def glbaln(self, s1, s2, alg):
        aln = mn.MMarcNW(s1, s2, -12, -4, self.dnamat, self.dnacol)
        aln.fill_matrix(alg)
        return aln.traceback(alg)

    def snp_search(self, s1, s2, snp, count):
        i = 0
        while count > 1:
            if s2[i] == '-':
                i += 1
            else:
                count -= 1
                i += 1
        if any(x == '-' for x in s2[i:i + 3]) or any(x == '-' for x in s1[i:i + 3]):
            return False
        snpq = s1[i:i + 3]
        snpq_aa = [self.codon_translate(x) for x in self.expand_ambig(snpq)]
        if any(x in snp for x in snpq_aa):
            return i

    def score_context(self, qstr, i, lcont, rcont):
        laa = expand_regex(list(lcont))
        raa = expand_regex(list(rcont))
        if any(x in self.ambig for x in qstr):
            return (None,) * 2
        snp_aa = self.codon_translate(qstr[i:i + 3])
        lind = 3 * (len(qstr[0:i]) // 3) if 3 * (len(qstr[0:i]) // 3) < 15 else 15
        rind = 3 * (len(qstr[i + 3:]) // 3) if 3 * (len(qstr[i + 3:]) // 3) < 15 else 15
        tldna = qstr[i - lind:i]
        trdna = qstr[i + 3:i + 3 + rind]
        max_val = -256
        max_prot = None
        if any(x in tldna for x in self.ambig.keys()) or any(x in trdna for x in self.ambig.keys()):
            ldna = self.expand_ambig(tldna)
            rdna = self.expand_ambig(trdna)
        else:
            ldna = [tldna]
            rdna = [trdna]
        for l in ldna:
            for r in rdna:
                if any(self.codon_translate(l[i:i + 3]) == 'Stop' for i in range(0, len(l), 3)):
                    return (None,) * 2
                elif any(self.codon_translate(r[i:i + 3]) == 'Stop' for i in range(0, len(r), 3)):
                    return (None,) * 2
                try:
                    cost = sum([self.get_mbcost(laa[5 - len(l) // 3:][e], self.codon_translate(l[i:i + 3])) for e, i in enumerate(range(0, len(l), 3))]) +\
                        sum([self.get_mbcost(raa[e], self.codon_translate(r[i:i + 3])) for e, i in enumerate(range(0, len(r), 3))])
                except IndexError:  # Limitated by context annotation length
                    left_length = int(min(len(l) / 3, len(laa)))
                    right_length = int(min(len(r) / 3, len(raa)))
                    cost = sum([self.get_mbcost(laa[5 - (3 * left_length) // 3:][e], self.codon_translate(l[i:i + 3])) for e, i in enumerate(range(0, 3 * left_length, 3))]) +\
                        sum([self.get_mbcost(raa[e], self.codon_translate(r[i:i + 3])) for e, i in enumerate(range(0, 3 * right_length, 3))])
                if cost > max_val:
                    max_val = cost
                    max_prot = ''.join([self.codon_translate(l[i:i + 3]) for e, i in enumerate(range(0, len(l), 3))]) + '|' + \
                        snp_aa + '|' + ''.join([self.codon_translate(r[i:i + 3]) for e, i in enumerate(range(0, len(r), 3))])
        return max_val, max_prot if max_prot else (None,) * 2

    def screen_hit(self, hit):
        if hit[1] in self.snp_models:
            self.snp_modelhit += 1
            if self.search_enabled:
                for m, snp in enumerate(self.snp_models[hit[1]]):
                    if hit[2] <= int(snp[3]) and hit[3] >= int(snp[4]):
                        self.snp_overlap += 1
                        try:
                            self.snpcov_correction[hit[1]][m][1] += 1
                        except KeyError:
                            self.snpcov_correction.setdefault(hit[1], [[0, 0] for x in range(len(self.snp_models[hit[1]]))])
                        if hit[4] > hit[5]:
                            hit = list(hit)
                            hit[0] = hit[0][::-1]
                            hit[4], hit[5] = hit[5], hit[4]
                        hmmdiff = hit[3] - hit[2]
                        alndiff = hit[5] - hit[4]
                        if abs(hmmdiff - alndiff) == 0 and self.fast:
                            lindex = hit[5] - (hit[3] - int(snp[3])) - self.dist \
                                if hit[5] - (hit[3] - int(snp[3])) - self.dist > hit[4] else hit[4]
                            rindex = hit[5] - (hit[3] - int(snp[4])) + self.dist \
                                if hit[5] - (hit[3] - int(snp[4])) + self.dist < hit[5] else hit[5]
                            offset = np.ceil((rindex - lindex - 2) / 2)
                            lsnp = int(int(snp[3]) - offset)
                            rsnp = int(int(snp[4]) + offset)
                            aln1, aln2, cost, g1, g2 = self.glbaln(hit[0][lindex-1:rindex], snp[-1][lsnp-1:rsnp], 'SW')
                        else:
                            lindex, rindex = hit[4], hit[5]
                            lsnp, rsnp = hit[2], hit[3]
                            offset = int(snp[3]) - lsnp
                            aln1, aln2, cost, g1, g2 = self.glbaln(hit[0][lindex-1:rindex], snp[-1][lsnp-1:rsnp], 'NW')
                        codon = self.snp_search(aln1, aln2, snp[2], int(offset)+1)
                        if codon:
                            score, prot = self.score_context(aln1, codon, snp[5], snp[6])
                            if score:
                                self.snp_found += 1
                                self.snpcov_correction[hit[1]][m][0] += 1
                                refstring = '|'.join([snp[5], snp[2], snp[6]])
                                if self.report:
                                    if self.verbose:
                                        self.report.write('\t'.join([str(x) for x in [hit[1], snp[0], snp[7], snp[8],
                                            hit[-1], cost, score, g1, g2, abs(hmmdiff-alndiff)]])+'\n')
                                    else:
                                        self.report.write('\t'.join([str(x) for x in [hit[1], snp[0], snp[7],
                                            snp[8]]])+'\n')
                                if self.print_align:  # For debugging purposes
                                    print(prot)
                                    print(refstring)
                                    print(aln1)
                                    print(aln2+'\n')
                                if self.fasta:
                                    self.fasta.write('>'+'|'.join([str(x) for x in [hit[1], snp[0], snp[7], snp[8],
                                            hit[-1], cost, score, g1, g2, abs(hmmdiff-alndiff)]])+'\n'+hit[0]+'\n')
            else:
                return False  # If SNP search is disabled, don't count any hits to the SNP models
        else:
            return True  # Count all hits not hitting to SNP models
        return True  # If SNP search is enabled, count all hits and correct with snpcov_correction at the end

    def report_counts(self):
        print('Hits to SNP models: {}'.format(self.snp_modelhit))
        if self.search_enabled:
            print('Hits overlapping known SNP codons: {}'.format(self.snp_overlap))
            print('Potential resistance SNPs detected: {}\n'.format(self.snp_found))
        else:
            print('SNP search not enabled\n')

    def close_all(self):
        if self.fasta:
            self.fasta.close()
        if self.report:
            self.report.close()

