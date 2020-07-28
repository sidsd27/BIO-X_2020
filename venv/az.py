import numpy as np
import pandas as pd
import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature


def CCTop(guide, off_target): #higher = bad
    score = 0
    for x in np.where(off_target != guide)[0]:
        score = score + 1.2 ** (x + 1)
    return score/224

def MITscore(guide, target): # high = good, low = bad

   # the MIT score as described in Hsu et al. and their website. The MIT Specificity score summarizes all off-targets into a single number from 0-100. The higher the number, the fewer off-target effects are expected.

    M = [0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804,
         0.685, 0.583]
    mm_pos: None = np.where(guide != target)[0]

    n_mm = len(mm_pos)

    if n_mm > 0:
        dtot = 0
        for pos1 in range(0, n_mm):
            for pos2 in range(pos1 + 1, n_mm):
                dtot += abs(mm_pos[pos1] - mm_pos[pos2])
        if n_mm > 1:
            davg = dtot / (n_mm * (n_mm - 1) / 2.)
        else:
            davg = 0.

        factor_distance = 1 / ((19 - davg) / 19 * 4 + 1) * 1. / n_mm ** 2

        factor_positions = 1
        for pos in mm_pos:
            factor_positions *= (1 - M[pos])

        score = factor_distance * factor_positions
    else:
        score = 1.0
    return score


def CFD(guide, target, guide_comp): #0 to 100 #high is bad
    mm_pos = np.where(guide != target)[0]
    guide_comp = np.where(guide_comp == 'T', 'U', guide_comp)
    guide = np.where(guide == 'T', 'U', guide)
    cfd_table = pd.read_excel('/Users/sidsd27/Desktop/TABLE19.xlsx')
    score = 1

    for MMpos in mm_pos:


        MMtype = 'r' + guide_comp[MMpos] + ":d" + target[MMpos].upper()

        score*=cfd_table[(cfd_table['Mismatch Type'] == MMtype) & (cfd_table['Position'] == MMpos + 1)]['Percent-Active'].iloc[0]

    return score


