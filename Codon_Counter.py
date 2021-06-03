import Bio
from Bio import SeqIO
from matplotlib import pyplot as plt

Query = ('TCA', 'TCG', 'TAG')            # ---> Codons which will be counted    #
FileName = ''                            # ---> Insert name of the GenBank file #

colors = {'TAG':'black',
          'TCA':'orange',
          'TCG':'red'}

bars_width = {'TAG': 3,
              'TCA': 3,
              'TCG': 3}

TargetCodons = {}

for codons in Query:
    TargetCodons[codons] = []

def ReverseComplement(x:str):
    revcomp={"A" : "T", "C" : "G", "G" : "C", "T" : "A",
             "R" : "Y", "Y" : "R", "S" : "S", "W" : "W", "K" : "M", "M" : "K",
             "B" : "V", "D" : "H", "H" : "D", "V" : "B",
             "N" : "N"}
    rc = str()
    for m in x[::-1]:
        rc = rc+revcomp[m]
    return(rc)
    
def Analyse_CDS(item, location):
    start, end = location._start, location._end
    length = end - start
    sequence = record.seq[start:end]
    strand = location._strand
    if strand == -1:
        sequence = ReverseComplement(sequence)
    if item.qualifiers['codon_start'] == ['1'] and length % 3 == 0:
        for index in range(length // 3):
            codon = sequence[3*index : 3*index+3]
            if codon in Query:
                if strand == 1:
                    TargetCodons[codon].append(start + 3*index + strand)
                else:
                    TargetCodons[codon].append(end - 3*index + 2*strand)
    else:
        print('Odd CDS found:', item)

for record in Bio.SeqIO.parse(FileName+'.gb', "genbank"):
    genome_length = len(record.seq)
    for item in record.features:
        if item.type == 'CDS':
            if isinstance(item.location, Bio.SeqFeature.FeatureLocation):
                Analyse_CDS(item, item.location)
            elif isinstance(item.location, Bio.SeqFeature.CompoundLocation):
                for component in item.location.parts:
                    Analyse_CDS(item, component)
                    
graph_width = genome_length/1000
graph_hegth = 5

for codons in Query:    
    fig_bars, ax_bars = plt.subplots(figsize = (graph_width , graph_hegth))
    ax_bars.set_xlim(0, genome_length)
    ax_bars.yaxis.set_visible(False)
    for item in TargetCodons[codons]:
        ax_bars.axvline(x = item, c = colors[codons], linewidth = bars_width[codons])
    fig_bars.savefig(FileName + '_' + codons + '_Barcode.svg', format = "svg")
    
for codon in TargetCodons:
    print(codon + ': ', len(TargetCodons[codon]))