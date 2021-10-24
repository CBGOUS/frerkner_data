from Bio import SeqIO

names=[]
records=[]
sequence_count=0
skip_count=0
duplicate_count=0
keep_count=0
for record in SeqIO.parse("/home/sr/research/ensembl/ensembl_hsa_3utrs.fa", "fasta"):
    sequence_count+=1
    if(sequence_count%10000==0):
        print(str(sequence_count)+"..")
        
    if not "unavailable" in record.seq:
        
        try:
            i= names.index(record.id.split("|")[2])
            #print(record.id.split("|")[2] + " already exists at index <" + str(i) + ">")
            duplicate_count+=1
            
        except ValueError:
            #print("add <" + str(record.id.split("|")[2]) + ">")
            names.append(record.id.split("|")[2])
            records.append(record)
            keep_count+=1
            
    else:
        skip_count+=1
        #print("missing sequence <" + record.id + "> - skipping")
#print names
    
with open("/home/sr/research/ensembl/ensembl_hsa_3utrs.uniq.fa", 'w') as f_out:

    r=SeqIO.write(records, f_out, 'fasta')    
    
print("read " + str(sequence_count) + " records:")
print("removed " + str(duplicate_count) + "duplicate entries")
print("removed " + str(skip_count) + " entries without sequence")
print("retained " + str(keep_count) + "records")
