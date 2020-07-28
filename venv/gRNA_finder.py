import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from classes import guide_RNA
from classes import guide_data



def find_gRNAs(feature, mRNA, gDNA, name):  # need to reconfigure a bit
    if feature == "Aligned":
        feature_sequence = align_prot(mRNA, gDNA)
        print(str(feature_sequence), gDNA)
    else:
        feature_sequence = find_feature_annotated(feature, mRNA) #Finds feature in mRNA file
    alignment = align_gDNA(feature_sequence, feature, gDNA, name)

    if alignment == "None":

        result = find_gRNAs("exon2", mRNA, gDNA, name)
        if result == None:
            result = find_gRNAs("N_terminus", mRNA, gDNA, name)
        print(result)
        return result
    gDNA = alignment[1] #gDNA sequence
    search_window = alignment[0] #window to search for guides
    final_list1 = idenitfy_elements(search_window, gDNA, feature_sequence, feature, 5, 57) #identify gRNA + features
    return final_list1, gDNA, feature

def find_feature_annotated(feature_name, mRNA_identifier):
        print(feature_name)
        Entrez.email = 'siddhant.suridhawan@gmail.com'
        handle_signal_sequence = Entrez.efetch(db="nucleotide", id=mRNA_identifier, rettype="gb", retmode="text")
        for seq_record in SeqIO.parse(handle_signal_sequence, "genbank"):
            if feature_name == "N_Terminus":
                return  seq_record.seq[:]
            if feature_name == "C_Terminus":
                return seq_record.seq[:]
            if feature_name == "CDS":
                return seq_record.seq[:]
            if feature_name == "Aligned":
                return seq_record.seq[:]

            else:
                f = (seq_record.features[index_genbank_features(seq_record, feature_name)])
                start = f.location.start
                end = f.location.end
                feature = (seq_record.seq[start:end])
        return feature

def index_genbank_features(gb_record, feature_type):


    if "N_Terminus" in feature_type:
        return "N_terminus"
    if "C_Terminus" in feature_type:
        return "C_Terminus"
    if "exon" in feature_type:
        exon_number = feature_type[len(feature_type) - 1]
        exon_index =  index_exon(gb_record, exon_number)
        return exon_index
    for (index, feature) in enumerate(gb_record.features):
        if feature.type == feature_type:
           return index

def index_exon(gb_record, exon_number): #need to adjust homology arm search for this
    i = 0
    for (index, feature) in enumerate(gb_record.features):
        if feature.type == "exon":
                i = i + 1
                if i == (int)(exon_number):
                    print("Finding gRNAs in exon", exon_number)
                    if index != None:
                       return index

def align_gDNA(feature, feature_name, gDNA_file, name):
    f = open("REPEAT_list.txt", 'a+')

    f.write(str(gDNA_file) + "  " + name)
    handle_signal_sequence = Entrez.efetch(db="nucleotide", id=gDNA_file, rettype="gbwithparts", retmode="text")
    for seq_record in SeqIO.parse(handle_signal_sequence, "genbank"):
            gDNA = seq_record.seq
            if feature_name == "Aligned":

                return feature, gDNA
            if feature_name == "N_Terminus":
                for (index, feature) in enumerate(seq_record.features):
                    if feature.type == "CDS":
                        return gDNA[feature.location.start: feature.location.start + 20], gDNA
            if feature_name == "C_Terminus":
                for (index, feature) in enumerate(seq_record.features):
                    if feature.type == "CDS":
                       return gDNA[feature.location.end - 200: feature.location.end], gDNA

            if feature.reverse_complement() in gDNA:
                gDNA = gDNA.reverse_complement()

            if feature in gDNA:
                print(feature)
                starts = gDNA.find(feature)
                ends = starts + len(feature)
                search_window = gDNA[ends - 50: ends+50]

                if "exon" in feature_name:
                    search_window = gDNA[starts - 20: starts + 100]
                    return search_window, gDNA

                return search_window, gDNA
            else:
                print("whoopsie")
                return "None"



def idenitfy_elements(search_window, gDNA, feature_sequence, feature_name, index=0, homology_length=57):
    gRNAs = []
    PAM = ["AGG", "TGG", "CGG", "GGG"]
    silent_mut_inserts = {"CGG": "CGA", "AGG": "AGA", "TGG": "TGC", "GGG": "GGA"}
    crRNAs = []
    for i in range(0, 90):
        sliding_window = search_window[20 + i: 23 + i]
        for x in PAM:
            if str(sliding_window) == x:
                gRNAs.append(i)

                mut = silent_mut_inserts[str(x)]
                hUS = gDNA[gDNA.find(search_window[i: i + 20]) + 20 - homology_length: gDNA.find(search_window[i: i + 20]) + 20]
                hDS = mut + gDNA[gDNA.find(search_window[i: i + 20]) + 23: gDNA.find(search_window[i: i + 20]) + 23 + homology_length -3]
                distance = gDNA.find(feature_sequence) + len(feature_sequence) - gDNA.find(
                    search_window[i: i + 20]) - 17
                if feature_name == "N_Terminus":
                    distance = gDNA.find(search_window)  - gDNA.find(
                    search_window[i: i + 20]) - 17
                if feature_name == "C_Terminus":
                    distance = gDNA.find(search_window)  - gDNA.find(
                    search_window[i: i + 20]) - 17
                crRNA = guide_RNA(search_window[i: i + 20], search_window[i + 20: i + 23], distance, hUS, hDS, 1)
                crRNAs.append(crRNA)

        i = i + 1

    PAMs = {"CCA", "CCT", "CCC", "CCG"}
    for i in range(0, 90):
        sliding_window = search_window[i: i + 3]
        for x in PAMs:
            if str(sliding_window) == x:
                gRNAs.append(i)
                # hUS = (search_window[i + 3- 57: i + 3])
                # if the gRNA binds the antisense strand, the PAM is a glycine on the sense strand, so
                # I mutate it to alanine (CGA (GCT) and add CCC (glycine) right after
                hUS = str(gDNA[gDNA.find(search_window[i + 3 : i +23]) - homology_length + 3: gDNA.find(search_window[i + 3 : i +23])]) + "CGA"
                hDS = "CCC" + gDNA[gDNA.find(search_window[i + 3: i + 23]):gDNA.find(search_window[i + 3: i + 23]) + homology_length -3]
                distance = gDNA.find(feature_sequence) + len(feature_sequence) - gDNA.find(
                    search_window[i+3: i + 23]) - 3
                if feature_name == "N_Terminus":
                    distance = gDNA.find(search_window)  - gDNA.find(
                    search_window[i: i + 20]) + 17
                if feature_name == "C_Terminus":
                    distance = gDNA.find(search_window)  - gDNA.find(
                    search_window[i: i + 20]) + 17
                crRNA = guide_RNA(search_window[i + 3: i + 23].reverse_complement(), search_window[i: i + 3], distance,
                                  hUS, hDS, -1)
                crRNAs.append(crRNA)
        i = i + 1
    final_list = []
    for x in crRNAs:
        if x.d_insert in range(-30, 30):
            final_list.append(x)

    return final_list

def get_genomic_sequence(mRNA):
    Entrez.email = 'siddhant.suridhawan@gmail.com'
    handle_signal_sequence = Entrez.efetch(db="sequences", id="NM_000949.7", rettype="gbwithparts", retmode="fasta")
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

    gDNA = ""
    for key in indices:
        value = indices[key]
        lower = value[:value.find("-")]
        upper = value[value.find("-") + 1:]
        Entrez.email = 'siddhant.suridhawan@gmail.com'
        handle_signal_sequence = Entrez.efetch(db="sequences", id=key, rettype="gb", retmode="fasta")
        for seq_record in SeqIO.parse(handle_signal_sequence, "genbank"):
            gDNA += str(seq_record.seq[(int)(lower):(int)(upper)])
            if len(str(seq_record.seq[(int)(lower):(int)(upper)])) == (int)(upper) - (int)(lower):
                print(seq_record.seq)

def align_prot(input, accession):
    handle_signal_sequence = Entrez.efetch(db="sequences", id=accession, rettype="gb", retmode="fasta")
    print(input, accession)
    for seq_record in SeqIO.parse(handle_signal_sequence, "genbank"):
        gDNA = seq_record.seq
        for (index, feature) in enumerate(seq_record.features):
            if feature.type == "CDS":
                translation = (feature.qualifiers['translation'])
                translation = translation[0]
                print(translation)

                if input in translation:
                    return gDNA[feature.location.start + len(input)*3 - 10: feature.location.start+ len(input)*3 + 20]