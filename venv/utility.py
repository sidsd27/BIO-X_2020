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

from az import CCTop, MITscore, CFD
from classes import guide_RNA, guide_data



def introduce_silent_mutations(sequence, num_mutations):
    seq = str(sequence)
    codons = [('TTT', 'TTC'), ('TTA', 'TTG'), ('TCU', 'TCA', 'TCT', 'TCG'), ('TAT', 'TAC'), ('TGC', 'TGT'),
              ('CTT', 'CTC', 'CTG', 'CTA'), ('CCT', 'CCG', 'CCA', 'CCG'), ('CAT', 'CAC'), ('CAA', 'CAG'),
              ('CGT', 'CGA', 'CGC' ,'CGG'), ('ATT', 'ATC', 'ATA') , ('ACU', 'ACC', 'ACA', 'ACG'), ('AAT', 'AAC'), ('AAA', 'AAG'),
              ('AGT', 'AGC'), ('AGA', 'AGG'), ('GTT', 'GTC', 'GTA', 'GTG'), ('GCT', 'GCC', 'GCA', 'GCG'), ('AAT', 'AAC'), ('AAA', 'AAG'),
              ('AGU', 'AGC'), ('AGA', 'AGG'), ('GTA', 'GTC', 'GTC', 'GTG'), ('GCT', 'GCC', 'GCA', 'GCG'), ('GAT', 'GAC'), ('GAA', 'GAG'),
              ('GGT', 'GGA', 'GGC', 'GGG')]

    while (num_mutations >= 0):
        length = len(seq)
        for i in range(-3, length):
            if (i%3 == 0):
                window = seq[i:i + 3]
                for x in codons:
                    if window in x:
                        for j in x:
                            if j != window:
                                seq = seq[0:i] +  str(j) + seq[i+3:]
                                num_mutations = num_mutations - 1
    return seq

def print_file(filename):
    with open(filename) as fp:
        line = fp.readline()
        cnt = 1
        while line:
            print(line)
            line = fp.readline()
            cnt += 1

def complement(Nucleotide):
    comp = []
    for i in Nucleotide:
        if i == "T":
            comp.append("A")
        if i == "A":
            comp.append("T")
        if i == "G":
            comp.append("C")
        if i == "C":
            comp.append("G")

    return ''.join(comp)

def translate(seq):
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    protein = ""
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        protein += table[codon]
    return protein

def align_prot(input):
    handle_signal_sequence = Entrez.efetch(db="sequences", id="NM_000949.7", rettype="gb", retmode="fasta")
    for seq_record in SeqIO.parse(handle_signal_sequence, "genbank"):
        gDNA = seq_record.seq
        for (index, feature) in enumerate(seq_record.features):
            if feature.type == "CDS":
                translation = (feature.qualifiers['translation'])
                print(translation)
                if input in translation:
                    return[gDNA[:translation.find(input)*3]]



