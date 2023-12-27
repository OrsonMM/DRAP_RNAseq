#usr/bin/pyhton3

#  Autor: Orson Mestanza 
#  The script is the first version 0.1

from time import time
import pandas as pd
import numpy as np
import sys


def header_fasta(open_fasta):

	header = []

	with open(open_fasta) as f:
		for line in f:
			if line[0] == ">":
				line = line.rstrip('\n')
				header.append(line[1:])

	return(header)

def tpm_fpkm(open_count):

	header_dict = {}
	# dict = { header: (TOTAL,FPKM,TPM) } structure of data

	with open(open_count) as f:
		for line in f:
			line = line.rstrip('\n')
			line = line.split('\t')
			if line[0] != "bundle_id":

				header_dict[line[1]] = int(line[4]),float(line[10]),float(line[14])

			else:
				pass

	return(header_dict)

def kegg_annotation(open_kegg):

	kegg_tax = {}

	kegg_other = {}

	with open(open_kegg) as f:
		for line in f:
			line = line.rstrip("\n")
			line = line.split("\t")

			if len(line[1]) != 0: 
				#print(line)
				if line[2] == 'Animals' and line[3] == 'Arthropods':

					new_1 = line[0].replace("user:","")[:-2]
					#print(new_1)
					if new_1 not in kegg_tax:
						kegg_tax[new_1] = (line[1],line[5],line[6])
						
					else:
						
						if float(kegg_tax[new_1][2]) > float(line[6]):
							pass
						elif float(kegg_tax[new_1][2]) < float(line[6]):
							kegg_tax[new_1] = (line[1],line[5],line[6])
						else:
							pass

				elif line[2] == 'Bacteria':
					
					new_2 = line[0].replace("user:","")[:-2]
					#print(new_2)
					if new_2 not in kegg_tax:
						kegg_tax[new_2] = (line[1],line[5],line[6])
					else:
						if float(kegg_tax[new_1][2]) > float(line[6]):
							pass
						elif float(kegg_tax[new_1][2]) < float(line[6]):
							kegg_tax[new_1] = (line[1],line[5],line[6])
						else:
							pass

				else:

					new_3 = line[0].replace("user:","")[:-2]
					#print(new_3)
					if new_3 not in kegg_other:
						kegg_other[new_3] = (line[2]+"_"+line[3],line[1],line[5],line[6])
					else:
						if float(kegg_other[new_3][3]) > float(line[6]):
							pass
						elif float(kegg_other[new_3][3]) < float(line[6]):
							kegg_other[new_3] = (line[2]+"_"+line[3],line[1],line[5],line[6])
						else:
							pass


	#print(len(kegg_tax),len(kegg_other))

	return(kegg_tax,kegg_other)

def resume_results(header,count,kegg):

	num_kegg = 0
	num_not = 0
	out = 0

	uniq_kegg = {}

	for i in header:       # This section is for transcripts annotated kegg[0] var. 
		if i in count:
			if i in kegg[0]:
				#print(kegg[0][i][0],count[i][0])
				if kegg[0][i][0] not in uniq_kegg:
					val_1 = int(count[i][0])
					uniq_kegg[kegg[0][i][0]] = val_1
				else:
					val_2 = int(count[i][0])
					uniq_kegg[kegg[0][i][0]] = uniq_kegg[kegg[0][i][0]]+val_2


				num_kegg += 1

			elif i in kegg[1]:
				
				#print(kegg[1][i])
				#print(str("This ID is not Arthropods or Bacteria"+"	"+i+str(kegg[1][i])))
				num_not += 1
			else:
				#print(str(i+"	"+"header dont have annotation"))
				out += 1
		else:
			#print(str("WARM ##### This ID is not in count"+i))
			pass

	for key, value in uniq_kegg.items():
		print(str(key+"\t"+str(value)))
	#print(uniq_kegg)
	#print(len(header),num_kegg,num_not,out)
	#print(str("TRANSCRIPTS WITH KEGG ANOTATIONS :"+" "+str(num_kegg)+"|"),str("TRANSCRIPTS WITH NO KEGG Bacteria OR Arthropods :"+" "+str(num_not)),str("TRANSCRIPTS WITH NO KEGG ANOTATIONS :"+" "+str(out)))

def main():
	t1 = time()
	null_0 = sys.argv[0]

	open_fasta = sys.argv[1]
	open_count = sys.argv[2]
	open_kegg = sys.argv[3]
	out_file = sys.argv[4]

	header = header_fasta(open_fasta)
	count = tpm_fpkm(open_count)
	kegg = kegg_annotation(open_kegg)

	resume_results(header,count,kegg)

if __name__ == "__main__":
	main()
