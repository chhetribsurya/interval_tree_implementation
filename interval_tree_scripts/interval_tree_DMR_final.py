import os
import scipy.stats
import json
import time
import intervaltree
import cPickle as pickle 

start_time = time.time()


####Input files for script to run are attached in the dir input_file_for script:
cluster_TF_file = "/home/surya/Desktop/scripts/data/chrom_impute_chromHMM_data/wgEncodeRegTfbsClusteredWithCellsV3.bed"  
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
				#(analogous to elements of the list, and elements of the dictionary, 
				#characters of the string array) my_tree[start:end]:
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

					master_TF_dict_return[TF].update({"significant_unique_dmr_hits":0})
					master_TF_dict_return[TF].update({"background_unique_dmr_hits":0})
					master_TF_dict_return[TF].update({"TotalTF_coveredBy_sig_unique_dmr_hits":0})
					master_TF_dict_return[TF].update({"TotalTF_coveredBy_bg_unique_dmr_hits":0})

					master_TF_dict_return[TF].update({"custom_overlap_list_sig":[]})
					master_TF_dict_return[TF].update({"custom_overlap_list_bg":[]})
					master_TF_dict_return[TF].update({"As_overlap_list_sig":[]})
					master_TF_dict_return[TF].update({"As_overlap_list_bg":[]})
					master_TF_dict_return[TF].update({"Bs_overlap_list_sig":[]})
					master_TF_dict_return[TF].update({"Bs_overlap_list_bg":[]})


			#json.dump(master_TF_dict_return, outfile)
			return(master_TF_dict_return)


master_TF_tree_dict = generate_interval_tree(cluster_TF_file, JSON_dict_file)
print "\n\nTime to parse the ENCODE_TF_cluster = ", time.time()-start_time

my_data = master_TF_tree_dict
# Store data (serialize)
with open('interval_tree_DMR_TF.pkl', 'wb') as handle:
    pickle.dump(my_data, handle, protocol=pickle.HIGHEST_PROTOCOL)

# Load data (deserialize)
# with open('filename.pickle', 'rb') as handle:
#     master_TF_tree_dict = pickle.load(handle)


##################

def genome_overlap_1(a, b):
	if (a[2] < b[1]) or (a[1] > b[2]):
		return False 
	else:
		return True 

#################





start_time = time.time()
print "\n\nStep 1.... significant dmr analysis..."


directory = "/home/surya/Desktop/scripts/data/metilene_dmr_data/sig_dmr"
if not os.path.exists(directory):
    os.makedirs(directory)
    #os.rmdir(prior_dir)
    #shutil.rmtree()


##Following piece of code prevents the append file mode to append the outputs to be appended
#to the same file. Check if file exists to prevent the append mode for each keys/TF value 
#being written to the file:
file_name_list = [key+".txt" for key in master_TF_tree_dict.iterkeys()]
for file_name in file_name_list:
	prior_file = directory + "/" +file_name
	if os.path.exists(prior_file):
		os.remove(prior_file)

Sig_count = 0
count = 0
with open(out_file, "w") as out_file:
	with open(dmr_file, "r") as file:
		for line in file:
			splitted = line.strip().split()
			chrom = splitted[0]
			start = splitted[1]
			end = splitted[2]
			qval = splitted[3]
			meth_perc = splitted[4]
			Sig_count +=1

			line_coords = [chrom, start, end, qval, meth_perc]
			A_tuple = (chrom, int(start),int(end), "\t".join(line_coords))


			for key,value in master_TF_tree_dict.iteritems():
				each_tf_tree_dict = master_TF_tree_dict[key]
				#name the output file with its keyword:
				file_name = key + ".txt"
				outfile_path = directory + "/" + file_name

				B_intervaltree_list = each_tf_tree_dict.get(chrom)

				if B_intervaltree_list:		
					B_overlap_list = B_intervaltree_list.search(int(start), int(end))	

					if B_overlap_list:
						master_TF_tree_dict[key]["significant_unique_dmr_hits"] += 1
						master_TF_tree_dict[key]["As_overlap_list_sig"].append(A_tuple[3])	

						for each in B_overlap_list:
							B_tuple_3 = each.data
							#Total TF covered by unique hits:
							master_TF_tree_dict[key]["TotalTF_coveredBy_sig_unique_dmr_hits"] += 1
							master_TF_tree_dict[key]["Bs_overlap_list_sig"].append(each.data)
							overlap_tuple_print = "%s\t%s\t%s" %(key, A_tuple[3], each.data)
							master_TF_tree_dict[key]["custom_overlap_list_sig"].append(overlap_tuple_print)
							#print overlap_tuple_print

							#Count for Total TF covered by unique hits | wc -l of unix:
							count +=1											
							with open(outfile_path, "a") as appendfile:
								appendfile.write(overlap_tuple_print + "\n")

		print "\n\nTotal line count(~wc -l of unix):", count

	print "Total unique dmr hit count for last TF in dict ::: ", each_tf_tree_dict["significant_unique_dmr_hits"]
	print "But, total background count for last TF in dict ::: ", each_tf_tree_dict["background_unique_dmr_hits"]

print "\nsignificant_dmr_overlap completed!!!!!"
print "Time for significant dmr analysis = ", time.time()-start_time




#####################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#####################



print "\n\nStep 2.... background dmr analysis....."
start_time = time.time()


directory = "/home/surya/Desktop/scripts/data/metilene_dmr_data/background_dmr"
if not os.path.exists(directory):
    os.makedirs(directory)
    #os.rmdir(prior_dir)
    #shutil.rmtree()

#Check if file exists to prevent the append mode for each keys/TF value being written to the file:
file_name_list = [key+".txt" for key in master_TF_tree_dict.iterkeys()]
for file_name in file_name_list:
	prior_file = directory + "/" +file_name
	if os.path.exists(prior_file):
		os.remove(prior_file)

Bg_count = 0
count = 0
custom_overlap_list_2 = []
with open(out_file_2, "w") as outfile:
	with open(background_dmr_file, "r") as file:
		for line in file:
			splitted = line.strip().split()
			chrom = splitted[0]
			start = splitted[1]
			end = splitted[2]
			Bg_count +=1
		
			line_coords = [chrom, start, end, qval, meth_perc]
			A_tuple = (chrom, int(start),int(end), "\t".join(line_coords))

			for key,value in master_TF_tree_dict.iteritems():
				each_tf_tree_dict = master_TF_tree_dict[key]
				#name the output file with its keyword:
				file_name = key + ".txt"
				outfile_path = directory + "/" + file_name

				B_intervaltree_list = each_tf_tree_dict.get(chrom)

				if B_intervaltree_list:		
					B_overlap_list = B_intervaltree_list.search(int(start), int(end))	

					if B_overlap_list:
						master_TF_tree_dict[key]["background_unique_dmr_hits"] += 1
						master_TF_tree_dict[key]["As_overlap_list_bg"].append(A_tuple[3])	

						for each in B_overlap_list:
							B_tuple_3 = each.data
							#Total TF covered by unique hits:
							master_TF_tree_dict[key]["TotalTF_coveredBy_bg_unique_dmr_hits"] += 1
							master_TF_tree_dict[key]["Bs_overlap_list_bg"].append(each.data)
							overlap_tuple_print = "%s\t%s\t%s" %(key, A_tuple[3], each.data)
							master_TF_tree_dict[key]["custom_overlap_list_bg"].append(overlap_tuple_print)
							#print overlap_tuple_print

							#Count for Total TF covered by unique hits | wc -l of unix:
							count +=1
							with open(outfile_path, "a") as appendfile:
								appendfile.write(overlap_tuple_print + "\n")
			
		print "\n\nTotal line count(~wc -l of unix):", count

	print "Total background count for last TF in dict constructed ::: ", each_tf_tree_dict["background_unique_dmr_hits"]
	print "Also, total significant count for last TF in dictionary::: ", each_tf_tree_dict["significant_unique_dmr_hits"]

print "\nbackground_dmr_overlap completed!!!!!"
print "Time for background dmr analysis = ", time.time()-start_time
print("\n\n")




#Append significant dmr and background dmr hits to list_array and
#calculate fisher exact test & bonferroni correction for multiple test correction:

with open(TF_fishers_enrichment, "w") as out_file:	
	FishTable = []
	for i in range(5):
		FishTable.append([])

	for dict_key,dict_value in master_TF_tree_dict.iteritems():
		FishTable[0].append(dict_key)
		FishTable[1].append(master_TF_tree_dict[dict_key]["significant_unique_dmr_hits"])
		FishTable[2].append(master_TF_tree_dict[dict_key]["background_unique_dmr_hits"])

	# FishTable[3] = [items for i, items in zip(FishTable[0], FishTable[1])]
	# FishTable[4] = [items for i, items in zip(FishTable[0], FishTable[1])]

	#Calculating the enrichment pvalue using scipy fishers exact test:
	original_sig_count = Sig_count 
	original_bg_count = Bg_count
	header = ("TF_name", "Significant_DMR_hits", "Background_dmr_Hits", "Enrichment_pval", "Enrichment_Bonferonni_Adj")
	out_file.write("\t".join(header) + "\n")
	print "\t".join(header)
	for i in range(len(FishTable[0])):
		fishers_pvalue = scipy.stats.fisher_exact([[FishTable[1][i], original_sig_count - FishTable[1][i]], [FishTable[2][i], original_bg_count - FishTable[2][i]]], 'greater')[1]
		FishTable[3].append(fishers_pvalue)
		bonferroni_correction = min(scipy.stats.fisher_exact([[FishTable[1][i], original_sig_count - FishTable[1][i]], [FishTable[2][i], original_bg_count - FishTable[2][i]]], 'greater')[1] * len(FishTable[0]), 1)
		FishTable[4].append(bonferroni_correction)
		fishers_test_result = (str(FishTable[0][i]), str(FishTable[1][i]), str(FishTable[2][i]), str(FishTable[3][i]), str(FishTable[4][i]))
		fisher_print = "\t".join(fishers_test_result)
		print fisher_print
		out_file.write(fisher_print + "\n")

print "\nEnd of job!!\n"



