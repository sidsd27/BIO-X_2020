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
from prediction_util import *
from training_util import *
import requests
from az import CFD, MITscore, CCTop
from gRNA_finder import find_gRNAs
import subprocess
global gDNA
import sys
import os
np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(precision=3)
np.set_printoptions(formatter={'float': '{: 0.2f}'.format})

def complement(Nucleotide):

    comp = []
    for i in Nucleotide:
        if i == "T" or i=='t':
            comp.append("A")
        if i == "A" or i=='a':
            comp.append("T")
        if i == "G" or i=='g':
            comp.append("C")
        if i == "C" or i=='c':
            comp.append("G")

    return ''.join(comp)

def Actual_genomic_hit(index):
    global output
    targets = []
    rows = BeautifulSoup(driver.page_source, 'lxml').find(id='myTable0').tbody.find_all('tr')[0:5]
    for i in range(5):
        text = ((rows[i].find("td", {"style": "white-space: nowrap;font-family:Courier New;"}).text))
        targets.append(text[:text.find("-")])
    output[index] = targets





def get_off_targets():
    global driver
    options = webdriver.ChromeOptions()
    options.add_argument('headless')
    options.add_argument('window-size=1920x1080')
    options.add_argument("disable-gpu")
    driver = webdriver.Chrome(ChromeDriverManager().install(), options=options)
    driver.set_page_load_timeout(20)
    global output
    output = {}
    with open('just_guides.txt', 'r') as inputs_file:
        lines = inputs_file.readlines()
        index = 0
        for line in lines:
            line = line[:-1]
            if line != '\n' and line != '':
                driver.get('https://cm.jefferson.edu/Off-Spotter/')
                time.sleep(5)
                driver.find_element_by_id('input20mers').send_keys(line[4:len(line) - 7])
                time.sleep(5)
                driver.find_element_by_id('submitmerSearch').click()
                time.sleep(5)
                Actual_genomic_hit(index)
                index = index + 1
    driver.close()
    return output


def calculateScore(guide, off_target):
    CFD_score = []
    MIT_score = []
    CCTop_score = []
    for x in off_target:
        x = np.array(list(x))
        guide_comp = complement((guide))
        print(guide_comp, "guide_comp")
        guide = np.array(list(guide))
        guide_comp = np.array(list(guide_comp))
        CFD_score.append(CFD(guide, x, guide_comp))
        MIT_score.append(MITscore(guide, x))
        CCTop_score.append(CCTop(guide, x))
    return [max(MIT_score), min(CCTop_score), max(CFD_score)]

def main():
    get_off_targets()
    results = get_data()
    matrix = results[0]
    guides = results[1]
    print(matrix)
    scores = []
    for i in range(0, (matrix.shape)[1]):
        scores.append(matrix[6][i])
    f = open("REPEAT_list.txt", 'a+')
    for i in range(0,15):
        index = scores.index(max(scores))
        f.write(str("gRNA:" + guides[index] + '\n'))
        fp = open("arms.txt", 'r')
        f.write(str(str(i) + "   " + guides[scores.index(max(scores))] + "  "+ (fp.readlines())[scores.index(max(scores))] + '\n'))
        scores[index] = 0

    f.close()




def print_file(filename):
    with open(filename) as fp:
        line = fp.readline()
        cnt = 1
        while line:

            line = fp.readline()
            cnt += 1

def get_data():
    off_target_dict = output
    print(off_target_dict)
    with open('just_guides.txt', 'r') as inputs_file:
        guides = inputs_file.readlines()

    inputs_file.close()
    with open('doench_scores.txt', 'r') as inputs_file:
        doench_score = inputs_file.readlines()

    inputs_file.close()
    matrix = np.zeros((7, len(doench_score)), dtype=float)
    print(len(guides))
    for i, guide in enumerate(guides):
        print(guide, i)
        seq_guide = guide
        guide = str(guide)
        matrix[0][i] = i

        scores = calculateScore(guide[4:len(guide)-8], off_target_dict[i])


        matrix[1][i] = (float)(doench_score[i]) #doench
        matrix[2][i] = (float)((scores[0])) #CCTop
        matrix[3][i] = (float)((scores[2])) #CFD
        matrix[4][i] = (float)((scores[1])) #MIT


        enzyme = 'esp'
        file = open('gDNA.txt', 'a+')
        seq = file.readline()

        if "+" in guide:
            guide = guide[4:len(guide)-5]

        if "-" in guide:
            guide = complement(str(guide[:guide.find("-")]))
            guide = guide[::-1]
            guide = guide[4:len(guide) - 2]

        DeepHF = effciency_predict(guide, 'esp')
        matrix[5][i] = (float)(DeepHF)
        matrix[6][i] = float((float)(matrix[1][i]) * 1.2 +  (float)(matrix[5][i]) * 1.2 - (float)(matrix[3][i]) * 0.25 + (float)(matrix[4][i]*0.15) - (float)(matrix[2][i]*0.15)) #CFD)

    return matrix, guides




main()
