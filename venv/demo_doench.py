import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '/Users/sidsd27/Documents/GitHub/Azimuth')
from azimuth import  model_comparison
import numpy as np
import os
import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature


def main():

    doench_predictions()




def doench_predictions():
   f_new = open("doench_scores.txt", 'a+')
   f_new.truncate(0)
   num_lines = 0
   for line in open('just_guides.txt'):
       if "+" in line or  "-" in line:
            num_lines +=1

   set = np.empty(num_lines , dtype = 'object')
   with open('just_guides.txt') as fp:
        cnt = 0
        line = fp.readline()
        while (not line.isspace()):
            if not "+" in line and not "-" in line:
                break
            else:
                print(line)
                if "+" in line:
                    line = line[:line.find("+")]

                if "-" in line:
                    line = str(Seq(line[:line.find("-")]).reverse_complement())

                set[cnt] = line

                line = fp.readline()
                cnt +=1



   predictions = model_comparison.predict(set, None, None)

   for i, prediction in enumerate(predictions):
       f_new.write(str(predictions[i]) + '\n')

main()