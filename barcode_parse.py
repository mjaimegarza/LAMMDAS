from Bio import SeqIO
from Bio.Seq import Seq
import pysam
import pickle
import time
start = time.time()

#assert len(set(A)) == len(A)/2, 'not pairwise' #make sure file is pairwise
#the assertion passed

#importing SAM file with the help of pysam
samfile = pysam.AlignmentFile( '/Users/matt/OneDrive/UCSF/PUBS/outv3_vfl.sam', 'rb' )
A = []
def SAM_parse(samfile):
	for counter, line in enumerate(samfile.fetch()): #start a loop through SAM file
		if counter % 2 == 0: #Skip every other line
			if line.reference_name is not None: #Make sure the read found a target reference
				ref = line.reference_name
				mut_pos = int((int(ref[-10:].split(':')[0]) - 5514) / 3) #extact the location from the name
				t = ref[-1:] #target mutation
				s = ref[-5:-4] #source variant
				if t == s: #if target and source mutation are the same set target mutation to WT
					t = 'WT'
				A.append((line.query_name, mut_pos, t)) #return a list A: (Header, position of mutation, mutation)
	return(A)

B = {}
def barcode_parse():
	#import fastaq file with the help of biopython
	for record in SeqIO.parse( '/Users/matt/OneDrive/UCSF/PUBS/180831_M02564_0268_000000000-BHYPT_PUBS/Undetermined_S0_L001_I1_001.fastq', 'fastq' ):
		B[record.id] = record.seq #create a dictionary called B: {barcode: (pos, mut)}
	return(B)

C = {}
def barcode_lookup(A, B):
	for line in A: #each line in list A
		#Create a dictionary called C that uses the header from list A (above) to lookup
		#the keys (pos, mut) in list B and return the final dictionary with {sequence: (pos, mut)}
		C[ str(B[line[0]]) ] = (line[1], line[2])
	return(C)

#run all of the functions and save final dictionary as pickle file
if __name__ == '__main__':
	SAM_parse(samfile)
	barcode_parse()
	pickle_out = open( 'VFL.pickle', 'wb' )
	pickle.dump(barcode_lookup(A, B), pickle_out)
	pickle_out.close()
end = time.time()
print(end - start)
