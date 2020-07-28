import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature

Entrez.email = 'siddhant.suridhawan@gmail.com'
handle_signal_sequence = Entrez.efetch(db="sequences", id="NM_000949.7", rettype="gb", retmode="fasta")
lines = ((handle_signal_sequence.readlines()))
indices = {}
for i in range (0, len(lines)):
    line = lines[i]
    if "PRIMARY" in line:
        i = i +1
        while(i <= len(lines)):
            if  "FEATURES" in lines[i]:
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

