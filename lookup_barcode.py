import pickle
from Bio import SeqIO
from Bio.Seq import Seq

bc_list = []
counter = 0
def lookup_barcode(ref):
    #import fastaq file with the help of biopython
    for record in SeqIO.parse( 'LAM-2-24_S46_R1_001.fq', 'fastq' ):
        if str(record.seq) in ref:
            values = ref[str(record.seq)]
            keys = str(record.seq)
            data = [keys, values]
            bc_list.append(data)
    return(bc_list)

with (open("../pubs.pickle", "rb")) as openfile:
    while True:
        try:
            reference = pickle.load(openfile)
            pickle_out = open( 'LAM-2-24_S46_R1_001_v2_list.pickle', 'wb' )
            pickle.dump(lookup_barcode(reference), pickle_out)
            pickle_out.close()
            #print(lookup_barcode(reference))
        except EOFError:
            break
