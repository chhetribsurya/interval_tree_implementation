#master_TF_tree_dict[key]["custom_overlap_list"].append(overlap_print)
						#master_TF_tree_dict[key]["Bs_overlap_list"].append(each.data)
						#print overlap_tuples
						# with open(outfile_path, "a") as appendfile:
						# 	appendfile.write(overlap_print + "\n")
import os
import scipy.stats
import json
import time
import intervaltree


start_time = time.time()

####Input files for script to run are attached in the dir input_file_for script:
cluster_TF_file = "/home/surya/Desktop/scripts/data/chrom_impute_chromHMM_data/wgEncodeRegTfbsClusteredWithCellsV3_tail.bed"  
dmr_file = "/home/surya/Desktop/scripts/data/metilene_dmr_data/DMR_data.bed"
background_dmr_file = "/home/surya/Desktop/scripts/data/metilene_dmr_data/Background_DMR.bed"

####Output files from the script run (provide the path location to save the file):
#Provide the file path location to save the output of ENCODE_TF_cluster file as JSON format (could be useful for future usage):
JSON_dict_file = "/home/surya/Desktop/scripts/data/metilene_dmr_data/ENCODE_TF_JSON_1.txt"
TF_fishers_enrichment = "/home/surya/Desktop/scripts/data/metilene_dmr_data/TF_fishers_enrichment_result_1.txt"
out_file = "/home/surya/Desktop/scripts/data/metilene_dmr_data/significant_dmr_overlaps_with_TF_1.bed"
out_file_2 = "/home/surya/Desktop/scripts/data/metilene_dmr_data/background_dmr_overlaps_with_TF_1.bed"



master_TF_dict_return = {}
header = ["chrom", "start", "end", "TF_name", "cell_line"]
head_print = "\t".join(header)


def generate_interval_tree(TF_file, Json_output_file):
	with open(cluster_TF_file,"r") as file:
		with open(JSON_dict_file, 'w') as outfile:
			for line in file:
				splitted = line.strip("\n").split("\t",4)
				chrom, start, end, TF = splitted[0], splitted[1], splitted[2], splitted[3]
				cell_line = splitted[4].split()[1]
				line_info = [chrom, str(start), str(end), str(TF)]

				#Generate a interval tree to implement binary search tree algorithm:
				#Data structure is Intervaltree here, with intervals as its elements 
				#(analogous to elements of the list, and elements of the dictionary, characters of the string array)
				#my_tree[start:end]:
				if TF in master_TF_dict_return.keys():
					if chrom in master_TF_dict_return[TF].keys():
						master_TF_dict_return[TF][chrom].appendi(int(start), int(end), "\t".join(line_info))
					else:
						master_TF_dict_return[TF][chrom] = intervaltree.IntervalTree()
						master_TF_dict_return[TF][chrom].appendi(int(start), int(end), "\t".join(line_info))
						#master_TF_dict_return[TF].update({chrom = intervaltree.IntervalTree()}
				else:
					master_TF_dict_return[TF] = {chrom: intervaltree.IntervalTree()}				
					master_TF_dict_return[TF][chrom].appendi(int(start), int(end), "\t".join(line_info))

					#unique_dmr_hits:
					master_TF_dict_return[TF].update({"significant_dmr_hits":0})
					#no. of TF covered:
					master_TF_dict_return[TF].update({"significant_dmr_hits_1":0})

					master_TF_dict_return[TF].update({"background_dmr_hits":0})
					master_TF_dict_return[TF].update({"custom_overlap_list":[]})
					master_TF_dict_return[TF].update({"As_overlap_list":[]})
					master_TF_dict_return[TF].update({"Bs_overlap_list":[]})


			#json.dump(master_TF_dict_return, outfile)
			return(master_TF_dict_return)



master_TF_tree_dict = generate_interval_tree(cluster_TF_file, JSON_dict_file)
print "\n\nTime to parse the ENCODE_TF_cluster = ", time.time()-start_time

		

##################

def genome_overlap_1(a, b):
	if (a[2] < b[1]) or (a[1] > b[2]):
		return False 
	else:
		return True 

#################




start_time = time.time()
print "\n\nStep 1.... significant dmr analysis..."


#Generate a interval tree to implement binary search tree algorithm:
#Data structure is Intervaltree here, with intervals as its elements (analogous to elements of the list)
count = 0
custom_overlap_list = [] 
#with open(out_file, "w") as out_file:
with open(dmr_file, "r") as file:
	for line in file:
		splitted = line.strip().split()
		chrom = splitted[0]
		start = splitted[1]
		end = splitted[2]
		qval = splitted[3]
		meth_perc = splitted[4]

		line_coords = [chrom, start, end, qval, meth_perc]
		A_tuple = (chrom, int(start),int(end), "\t".join(line_coords))


		for key,value in master_TF_tree_dict.iteritems():
			each_tf_tree_dict = master_TF_tree_dict[key]
			#name the output file with its keyword:
			# file_name = key + ".txt"
			# outfile_path = directory + "/" + file_name

			B_intervaltree_list = each_tf_tree_dict.get(chrom)
			if B_intervaltree_list:	
				if B_intervaltree_list.overlaps_range(int(start), int(end)):
						#unique_dmr_hits:
						master_TF_tree_dict[key]["significant_dmr_hits_1"] += 1
						#for single cpgs:
						#mytree.overlaps_point(single_cpg_coord i.e int(start))

				#B_overlap_list = B_intervaltree_list.search(153770108, 154027427)
				B_overlap_list = B_intervaltree_list.search(int(start), int(end))
				#For single cpg overlaps:
				# mytree.search(point)
				# mytree.search(5)
				# Out[27]: {Interval(4, 7, (4, 7)), Interval(5, 9, {5: 9})}
				# or, 

				if B_overlap_list:
					master_TF_tree_dict[key]["As_overlap_list"].append(A_tuple[3])					
					#master_TF_tree_dict[key]["significant_dmr_hits"] += 1
					for each in B_overlap_list:
						B_tuple_3 = each.data
						#unique_dmr_hits:
						master_TF_tree_dict[key]["significant_dmr_hits"] += 1
						#print key, A_tuple[3], each.data
						overlap_print = "%s\t%s\t%s" %(key, A_tuple[3], B_tuple_3)
						print overlap_print
						#total_no. of TF covered:
						count +=1

						#master_TF_tree_dict[key]["custom_overlap_list"].append(overlap_print)
						#master_TF_tree_dict[key]["Bs_overlap_list"].append(B_tuple_3)
						#print overlap_tuples
						# with open(outfile_path, "a") as appendfile:
						# 	appendfile.write(overlap_print + "\n")

	print "count:", count
 

	print "Total sig count for last TF in dict constructed ::: ", each_tf_tree_dict["significant_dmr_hits"]
	print "But, total background count for last TF in dict ::: ", each_tf_tree_dict["background_dmr_hits"]

print "\nsignificant_dmr_overlap completed!!!!!"
print "Time for significant dmr analysis = ", time.time()-start_time



#####################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#####################


				#B_overlap_list = B_intervaltree_list.search(153770108, 154027427)
				B_overlap_list = B_intervaltree_list.search(start, end)
				if B_overlap_list:
					for each in B_overlap_list:
						#print A_tuple[3], each.data
						overlap_print = "%s\t%s" %(A_tuple[3], each.data)
						print overlap_print
 
						#overlap_tuples = "%s\t%s\t%s" %(key, A_tuple[3], B_tuple[3])
						#custom_overlap_list.append(overlap_tuples)

						#master_TF_tree_dict[key]["custom_overlap_list"].append(overlap_tuples)
						#master_TF_tree_dict[key]["As_overlap_list"].append(A_tuple[3])
						#master_TF_tree_dict[key]["Bs_overlap_list"].append(B_tuple[3])
						#print overlap_tuples

						# with open(outfile_path, "a") as appendfile:
						# 	appendfile.write(overlap_tuples + "\n")

	print "Total sig count for last TF in dict constructed ::: ", each_tf_tree_dict["significant_dmr_hits"]
	print "But, total background count for last TF in dict ::: ", each_tf_tree_dict["background_dmr_hits"]

print "\nsignificant_dmr_overlap completed!!!!!"
print "Time for significant dmr analysis = ", time.time()-start_time


				if B_overlap_list:
					master_TF_tree_dict[key]["As_overlap_list"].append(A_tuple[3])					
					#master_TF_tree_dict[key]["significant_dmr_hits"] += 1
					for each in B_overlap_list:
						#unique_dmr_hits:
						master_TF_tree_dict[key]["significant_dmr_hits"] += 1
						print A_tuple[3], each.data


Example:					
chr5	118664517	118666699	6.4e-11	0.623359	chr5	118666360	118666916	HDAC8
chr5	118664517	118666699	6.4e-11	0.623359	chr5	118664444	118665000	HDAC8
chr5	118664517	118666699	6.4e-11	0.623359	chr5	118666193	118666589	HDAC8



#####################
#####################


print "\n\nStep 2.... background dmr analysis....."
start_time = time.time()

##Following piece of code prevents the append file mode to append the outputs to be appended
#to the same file:
directory = "/home/surya/Desktop/scripts/data/metilene_dmr_data/python_dir_background_dmr"
if not os.path.exists(directory):
    os.makedirs(directory)
    #os.rmdir(prior_dir)
    #shutil.rmtree()


#Check if file exists to prevent the append mode for each keys/TF value being written to the file:
file_name_list = [key+".txt" for key in master_TF_dict.iterkeys()]
for file_name in file_name_list:
	prior_file = directory + "/" +file_name
	if os.path.exists(prior_file):
		os.remove(prior_file) 


custom_overlap_list_2 = []
with open(out_file_2, "w") as outfile:
	with open(background_dmr_file, "r") as file:
		for line in file:
			splitted = line.strip().split()
			chrom = splitted[0]
			start = splitted[1]
			end = splitted[2]

			line_coords = [chrom, start, end]
			A_tuple = (chrom, int(start),int(end), "\t".join(line_coords))


			for key,value in master_TF_dict.iteritems():
				each_tf_dict = master_TF_dict[key]
				
				#name the output file with its keyword:
				file_name = key + ".txt"
				outfile_path = directory + "/" + file_name

				B_tuple_list = each_tf_dict.get(chrom)       
				if B_tuple_list:
					for B_tuple in B_tuple_list:               
						if genome_overlap_1(A_tuple, B_tuple):
								master_TF_dict[key]["background_dmr_hits"] += 1
								overlap_tuples = "%s\t%s\t%s" %(key, A_tuple[3], B_tuple[3])
								custom_overlap_list_2.append(overlap_tuples)

								master_TF_dict[key]["custom_overlap_list"].append(overlap_tuples)
								master_TF_dict[key]["As_overlap_list"].append(A_tuple[3])
								master_TF_dict[key]["Bs_overlap_list"].append(B_tuple[3])
								#print overlap_tuples

								with open(outfile_path, "a") as appendfile:
									appendfile.write(overlap_tuples + "\n")

		print "Total background count for last TF in dict constructed ::: ", each_tf_dict["background_dmr_hits"]
		print "Also, total significant count for last TF in dictionary::: ", each_tf_dict["significant_dmr_hits"]


print "\nbackground_dmr_overlap completed!!!!!"
print "Time for background dmr analysis = ", time.time()-start_time
print("\n\n")



#Append significant dmr and background dmr hits to list_array and
#calculate fisher exact test & bonferroni correction for multiple test correction:

with open(TF_fishers_enrichment, "w") as out_file:	
	FishTable = []
	for i in range(5):
		FishTable.append([])

	for dict_key,dict_value in master_TF_dict.iteritems():
		FishTable[0].append(dict_key)
		FishTable[1].append(master_TF_dict[dict_key]["significant_dmr_hits"])
		FishTable[2].append(master_TF_dict[dict_key]["background_dmr_hits"])

	# FishTable[3] = [items for i, items in zip(FishTable[0], FishTable[1])]
	# FishTable[4] = [items for i, items in zip(FishTable[0], FishTable[1])]

	#Calculating the enrichment pvalue using scipy fishers exact test:
	Significant_hits = len(custom_overlap_list)
	Background_hits = len(custom_overlap_list_2)
	header = ("TF_name", "Significant_DMR_hits", "Background_dmr_Hits", "Enrichment_pval", "Enrichment_Bonferonni_Adj")
	out_file.write("\t".join(header) + "\n")
	print "\t".join(header)
	for i in range(len(FishTable[0])):
		fishers_pvalue = scipy.stats.fisher_exact([[FishTable[1][i], Significant_hits - FishTable[1][i]], [FishTable[2][i], Background_hits - FishTable[2][i]]], 'greater')[1]
		FishTable[3].append(fishers_pvalue)
		bonferroni_correction = min(scipy.stats.fisher_exact([[FishTable[1][i], Significant_hits - FishTable[1][i]], [FishTable[2][i], Background_hits - FishTable[2][i]]], 'greater')[1] * len(FishTable[0]), 1)
		FishTable[4].append(bonferroni_correction)
		fishers_test_result = (str(FishTable[0][i]), str(FishTable[1][i]), str(FishTable[2][i]), str(FishTable[3][i]), str(FishTable[4][i]))
		fisher_print = "\t".join(fishers_test_result)
		print fisher_print
		out_file.write(fisher_print + "\n")

print "\nEnd of job!!\n"




In [2]: mytree = intervaltree.IntervalTree()

In [3]: type(mytree)
Out[3]: intervaltree.intervaltree.IntervalTree

In [4]: mytree
Out[4]: IntervalTree()

In [5]: mytree[1:2] = "55-89"

In [6]: mytree
Out[6]: IntervalTree([Interval(1, 2, '55-89')])

In [7]: mytree[4:7] = (4,7)

In [8]: mytree[5:9] = {5:9}

In [9]: mytree
Out[9]: IntervalTree([Interval(1, 2, '55-89'), Interval(4, 7, (4, 7)), Interval(5, 9, {5: 9})])



In [10]: sorted(t[6])
---------------------------------------------------------------------------
NameError                                 Traceback (most recent call last)
<ipython-input-10-228f05cfff21> in <module>()
----> 1 sorted(t[6])

NameError: name 't' is not defined

In [11]: sorted(mytree[6])
Out[11]: [Interval(4, 7, (4, 7)), Interval(5, 9, {5: 9})]

In [12]: sorted(mytree[67])
Out[12]: []

In [13]: sorted(mytree[5])
Out[13]: [Interval(4, 7, (4, 7)), Interval(5, 9, {5: 9})]

In [14]: sorted(mytree[4])
Out[14]: [Interval(4, 7, (4, 7))]







In [24]: mytree.overlaps_point(5)
Out[24]: True

In [25]: mytree.overlaps_point(111)
Out[25]: False

In [26]: mytree.overlaps_range(2,5)
Out[26]: True

In [27]: mytree.search(5)
Out[27]: {Interval(4, 7, (4, 7)), Interval(5, 9, {5: 9})}

In [28]: mytree.search(5, 10)
Out[28]: {Interval(4, 7, (4, 7)), Interval(5, 9, {5: 9})}

In [29]: mytree.search(1, 5)
Out[29]: {Interval(1, 2, '55-89'), Interval(4, 7, (4, 7))}


#Shows the search is exclusive/not inclusive (since 5(end coord) wont be included in the search):
In [63]: mytree.appendi(3,5, {"dict": "cool"})

In [64]: mytree
Out[64]: IntervalTree([Interval(1, 2, '55-89'), Interval(1, 5, {'ram': 'cool'}), Interval(3, 5, {'dict': 'cool'}), Interval(3, 9, {'dict': 'cool'}), Interval(4, 7, (4, 7)), Interval(5, 9, {5: 9})])

In [65]: mytree.search(5)
Out[65]: 
{Interval(3, 9, {'dict': 'cool'}),
 Interval(4, 7, (4, 7)),
 Interval(5, 9, {5: 9})}

In [66]: mytree.appendi(3,6, {"dict": "cool"})

In [67]: mytree.search(5)
Out[67]: 
{Interval(3, 6, {'dict': 'cool'}),
 Interval(3, 9, {'dict': 'cool'}),
 Interval(4, 7, (4, 7)),
 Interval(5, 9, {5: 9})}



In [30]: for my in mytree:
   ....:     print my
   ....:     
Interval(1, 2, '55-89')
Interval(5, 9, {5: 9})
Interval(4, 7, (4, 7))

In [31]: for my in mytree:
    print my.begin
   ....:     
1
5
4

In [32]: for my in mytree:
    print my.data
   ....:     
55-89
{5: 9}
(4, 7)
