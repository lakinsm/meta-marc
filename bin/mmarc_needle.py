## Module for small-scale Needlman-Wunch pairwise alignments for the Meta-MARC package

import numpy as np


class MMarcNW(object):
    def __init__(self, seq1, seq2, gap, extend, submat, subcol):
        self.s1 = seq1
        self.s2 = seq2
        self.as1 = ''
        self.as2 = ''
        self.rows = len(seq1) + 1
        self.cols = len(seq2) + 1
        self.f = np.zeros((self.rows, self.cols), dtype=np.int32)
        self.g = int(gap)
        self.e = int(extend)
        self.submat = submat
        self.subcol = subcol
        self.total_cost = None
        self.gap1 = 0
        self.gap2 = 0

    def _edit_cost(self, a, b):
        return self.submat[self.subcol.index(a), self.subcol.index(b)]

    def fill_matrix(self, algorithm):
        if algorithm == 'NW':
            for i in range(self.rows):
                self.f[i, 0] = self.g * i + 1
            for j in range(self.cols):
                self.f[0, j] = self.g * j + 1
            for i in range(1, self.rows):
                for j in range(1, self.cols):
                    diag = self.f[i-1, j-1] + self._edit_cost(self.s1[i-1], self.s2[j-1])
                    left = self.f[i-1, j] + self.g
                    up = self.f[i, j-1] + self.g
                    self.f[i, j] = max(diag, left, up)
        elif algorithm == 'SW':
            for i in range(self.rows):
                self.f[i, 0] = self.g * i + 1
            for j in range(self.cols):
                self.f[0, j] = self.g * j + 1
            for i in range(1, self.rows):
                for j in range(1, self.cols):
                    diag = self.f[i-1, j-1] + self._edit_cost(self.s1[i-1], self.s2[j-1])
                    left = self.f[i-1, j] + self.g
                    up = self.f[i, j-1] + self.g
                    self.f[i, j] = max(0, diag, left, up)
        else:
            raise ValueError('Choice of algorithm must be NW (Needleman-Wunsch) or SW (Smith-Waterman)')

    def traceback(self, algorithm):
        i = len(self.s1)
        j = len(self.s2)
        self.total_cost = self.f[i, j]
        while i > 0 and j > 0:
            current_val = self.f[i, j]
            diag_val = self.f[i-1, j-1]
            up_val = self.f[i, j-1]
            left_val = self.f[i-1, j]
            if algorithm == 'NW':
                if current_val == diag_val + self._edit_cost(self.s1[i-1], self.s2[j-1]):
                    self.as1 += self.s1[i-1]
                    self.as2 += self.s2[j-1]
                    i -= 1
                    j -= 1
                elif current_val == up_val + self.g:
                    self.as1 += '-'
                    self.as2 += self.s2[j-1]
                    self.gap1 += 1
                    j -= 1
                elif current_val == left_val + self.g:
                    self.as1 += self.s1[i-1]
                    self.as2 += '-'
                    self.gap2 += 1
                    i -= 1
                else:
                    raise ValueError('NW matrix filled incorrectly')
            elif algorithm == 'SW':
                opt = max(diag_val, up_val, left_val)
                if opt == diag_val:
                    self.as1 += self.s1[i-1]
                    self.as2 += self.s2[j-1]
                    i -= 1
                    j -= 1
                elif opt == up_val:
                    self.as1 += '-'
                    self.as2 += self.s2[j-1]
                    self.gap1 += 1
                    j -= 1
                elif opt == left_val:
                    self.as1 += self.s1[i-1]
                    self.as2 += '-'
                    self.gap2 += 1
                    i -= 1
                else:
                    raise ValueError('SW matrix filled incorrectly')
            else:
                raise ValueError('Choice of algorithm must be NW (Needleman-Wunsch) or SW (Smith-Waterman)')
        while i > 0:
            self.as1 += self.s1[i-1]
            self.as2 += '-'
            i -= 1
        while j > 0:
            self.as1 += '-'
            self.as2 += self.s2[j-1]
            j -= 1
        return self.as1[::-1], self.as2[::-1], self.total_cost, self.gap1, self.gap2


if __name__ == '__main__':
    import sys
    seq1 = str(sys.argv[1])
    seq2 = str(sys.argv[2])
    print(seq1)
    print(seq2)
    nw = MMarcNW(seq1, seq2, -12, -4)
    nw.fill_matrix('NW')
    as1, as2, cost, g1, g2 = nw.traceback('NW')
    print(as1)
    print(as2)
    print(cost, g1, g2)




