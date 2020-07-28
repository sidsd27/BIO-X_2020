import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
import numpy as np
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.keys import Keys
import pandas as pd
from selenium.webdriver.common.desired_capabilities import DesiredCapabilities
import os
from bs4 import BeautifulSoup
import lxml
from webdriver_manager.chrome import ChromeDriverManager
import time
from subprocess import call

import requests
from az import CFD, MITscore, CCTop
from gRNA_finder import find_gRNAs
import subprocess
global gDNA
global set_gRNAs
import sys
import os
np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(precision=3)
np.set_printoptions(formatter={'float': '{: 0.2f}'.format})
global feature_name
def main():
    feature_name = "sig_peptide"
    mRNA = (str(sys.argv[1]))
    gDNA = (str(sys.argv[2]))
    name = (str(sys.argv[3]))
    if mRNA.startswith("M"):
        feature_name = "Aligned"
    print(mRNA, gDNA)
    result = find_gRNAs(feature_name, mRNA, gDNA, name)
    print(result)
    feature_name = feature_name
    set_gRNAs = result[0]
    if set_gRNAs == None:
        print("failed")
    print(set_gRNAs)


    f = open("just_guides.txt", 'a+')
    f.truncate(0)
    gDNA = result[1]
    fp = open("gDNA.txt", 'a+')
    fp.truncate(0)
    fp.write(str(gDNA))
    fp.close()
    fp = open("arms.txt", 'a+')
    fp.truncate(0)
    feature_name = result[2]
    for x in set_gRNAs:
        if len(x.crRNA)== 20:
            if gDNA.find(x.crRNA) > 0 :
                f.write(str(gDNA[gDNA.find(x.crRNA) - 4: gDNA.find(x.crRNA) + 26]) + "+"+ "\n")
                fp.write(str(x.hUS) + "  " +  str(x.hDS) + "  " +  str(x.d_insert) + " "+ feature_name  + "\n")
            else:
                crRNA = (x.crRNA).reverse_complement()
                if crRNA in gDNA:

                     f.write(str(gDNA[gDNA.find(crRNA) - 6: gDNA.find(crRNA) + 24]) + "-"+ "\n")
                     fp.write(str(x.hUS) + "  " + str(x.hDS) + "  " + str(x.d_insert) + " " + feature_name + "\n" )

                else:
                    print("Something went wrong, guide RNA not found:", x.crRNA)
    f.close()

def get_genomic_sequence(mRNA):
    Entrez.email = 'siddhant.suridhawan@gmail.com'
    handle_signal_sequence = Entrez.efetch(db="sequences", id="NM_000949.7", rettype="gb", retmode="fasta")
    lines = ((handle_signal_sequence.readlines()))

    indices = {}
    for i in range(0, len(lines)):
        line = lines[i]
        if "PRIMARY" in line:
            i = i + 1
            while (i <= len(lines)):
                if "FEATURES" in lines[i]:
                    break
                line = lines[i].split(" ")
                while ("" in line):
                    line.remove("")
                if len(line) >= 2:
                    if "\n" in line[2]:
                        line[2] = (line[2])[: line[2].find("\n")]
                    indices[line[1]] = (line[2])
                i = i + 1
    print(indices)
    gDNA = ""
    for key in indices:
        value = indices[key]
        Entrez.email = 'siddhant.suridhawan@gmail.com'
        handle_signal_sequence = Entrez.efetch(db="sequences", id=key, rettype="gb", retmode="fasta")
        for seq_record in SeqIO.parse(handle_signal_sequence, "genbank"):
            lower = value[:value.find("-")]
            upper = value[value.find("-") + 1:]
            print(lower, upper)
            gDNA = gDNA + seq_record.seq[(int)(lower):(int)(upper)]
    print(gDNA)

main()


