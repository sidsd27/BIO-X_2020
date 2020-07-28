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
from classes import guide_RNA, guide_data
import os
from bs4 import BeautifulSoup
import lxml
from webdriver_manager.chrome import ChromeDriverManager
import time
import timeit



def get_data(set_gRNAs, gDNA, f):
    lines = []
    options = webdriver.ChromeOptions()
    options.add_argument('headless')
    options.add_argument('window-size=1920x1080')
    options.add_argument('--no-sandbox')
    cnt = 1
    for i in set_gRNAs:
        cnt = cnt +1
        get_off_targets()
        lines.append(i)
        data = []
    matrix = np.zeros((6, cnt -1), dtype=float)
    for i, guide in enumerate(lines):
        seq_guide = guide
        guide = str(guide.crRNA)
        matrix[0][i] = i
        f.write(f'{i + 1}/{len(lines)}')
        driver = webdriver.Chrome(ChromeDriverManager().install(), options=options)
        f.write("\n"+ "crRNA sequence:" + guide + "\n")
        driver.get(
                'https://design.synthego.com/#/validate/results?genome=homo_sapiens_gencode_26_primary&nuclease=cas9&guide=' + guide.replace(
                    '\n', ''))
        time.sleep(15)
        page_source = BeautifulSoup(driver.page_source, 'lxml')
        if(page_source ==None):
            break
        #text = BeautifulSoup(driver.page_source, 'lxml').find(class_='results-text-summary')
        if(text!=None):
            f.write(str(seq_guide.d_insert) + ": " + "Distance from target site \n")
            f.write(text.h2.text + "\n")
            temp = page_source.tbody.tr.find_all('td')
            radio_buttons = []
            #for t in temp:
               # radio_buttons.append('1' if t.i.text == 'check_circle' else '0')
            #doench_score = temp[2]['uib-tooltip']

            temp = page_source.find(
                    class_='table table-hover table-striped table-vertical table-offtargets ng-scope').find_all(
                    class_='mono ng-binding')
            off_target = []
            for t in temp:
                if guide != t.text:
                   off_target.append(t.text)
            f.write("Upstream homology arm:" + str(seq_guide.hUS) + "\n" + "Length:" + str(len(str(seq_guide.hUS))) + "\n")
            f.write("Downstream homology arm:" + str(seq_guide.hDS) + "\n" + "Length:" + str(len(str(seq_guide.hDS)))+ "\n")
            f.write(doench_score + "\n")
            scores = calculateScore(guide, off_target)
            #matrix[1][i] = (float)(doench_score.split(" ")[2]) #doench
            matrix[2][i] = (float)((scores[0])) #CCTop
            matrix[3][i] = (float)((scores[2])) #CFD
            matrix[4][i] = (float)((scores[1])) #MIT

            f.write(str(scores[0]) + ": MIT score\n")
            f.write(str(scores[2]) + ": CFD score\n")
            f.write(str(scores[1])  + ": CCTop score\n")
            enzyme = 'esp'
            seq = str(gDNA)
            if seq_guide.strand == 1:
               DeepHF = effciency_predict(str(seq_guide.crRNA) + str(seq_guide.PAM), 'esp')
               f.write(str(DeepHF) + ": DeepHF score")
            else:
              DeepHF = effciency_predict((str(seq_guide.PAM) + str(seq_guide.crRNA)), 'esp')
              f.write(str(DeepHF) + ": DeepHF score")
            matrix[5][i] = (float)(DeepHF)
        else:
            print("aborted", guide)
    driver.close()
    f.close()
    return [f, matrix]