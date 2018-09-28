#Opens R1 and R2 files and reads them per line
R1 = open("Undetermined_S0_L001_R1_001.fastq", "r").readlines()
R2 = open("Undetermined_S0_L001_R2_001.fastq", "r").readlines()

def assertion(R1, R2):
    assert (len(R1) == len(R2)), 'R1 and R2 are not the same size'
    #step across each 4th line which corresponds to each sequence header
    for i in range(0,len(R1),4):
        #This gets rid of the unique string at the end of the header and keeps the non-unique part
        R1_kept, partition, removed = R1[i].partition(' ')
        R2_kept, partition, removed = R2[i].partition(' ')
        assert (R1_kept == R2_kept), 'R1 and R2 are not in the same order'

assertion(R1, R2)
