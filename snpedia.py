import sys,pickle,snp
from wikitools import wiki,category,page
from glob import glob


def snpedia_items():										# get all snps from snpedia
	try:
		snpedia_infile = open("snpedia.data","rb")			# try getting cached snps 
	except IOError:
		snpedia = snpedia_getter()							# if none avail get from snpedia
	else:	
		snpedia = pickle.load(snpedia_infile)				# else load from cache
		snpedia_infile.close()
		
	return snpedia	

		
def snpedia_getter():
	site = wiki.Wiki("http://snpedia.com/api.php")			# open snpedia
	snps = category.Category(site, "Is_a_snp")
	snpedia = {}
	
	for article in snps.getAllMembersGen(namespaces=[0]):	# get all snp-names
		snpedia[article.title.lower()] = "in snpedia"
		print article.title
		
	snpedia_outfile = open("snpedia.data","wb")				# save all snps to cache
	pickle.dump(snpedia,snpedia_outfile)
	snpedia_outfile.close()
	return snpedia


def snp_compare(snpedia,my_snps):							# compare my snps to snpedia 
	filtered_snps = []
	for single_snp in my_snps:
		if len(single_snp.genotype) >= 2:					# only use homozygotous snps
			if single_snp.genotype[0] == single_snp.genotype[1]:
				if snpedia.has_key(single_snp.name):
					filtered_snps.append(single_snp)
	return filtered_snps 


def snpedia_genotypes(my_filtered_snps):					# get genotypes for snps from snpedia
	try:
		genotype_infile = open("genotypes.data","rb")		# again, try to cache data to save load
	except IOError:
		genotypes = genotype_getter(my_filtered_snps)
	else:
		genotypes = pickle.load(genotype_infile)
		genotype_infile.close()
		
	return genotypes 


def genotype_getter(my_filtered_snps):
	site = wiki.Wiki("http://snpedia.com/api.php")
	genotypes = {}
	
	for single_snp in my_filtered_snps:
		type_counter = 1
		wikipage = page.Page(site,single_snp.name)
		snp_page = wikipage.getWikiText()
		
		while snp_page.find("geno"+str(type_counter)) != -1:
			
			if genotypes.has_key(single_snp.name):
				current_genotypes = genotypes[single_snp.name]
				type_start = snp_page.find("geno"+str(type_counter))
				type_start = snp_page.find("(",type_start)
				type_stop = snp_page.find(")",type_start)
				current_genotypes.append(str(snp_page[type_start:type_stop+1])) 
				genotypes[single_snp.name] = current_genotypes
				
			else:
				type_start = snp_page.find("geno"+str(type_counter))
				type_start = snp_page.find("(",type_start)
				type_stop = snp_page.find(")",type_start)
				genotypes[single_snp.name] = [str(snp_page[type_start:type_stop+1])]
				
			type_counter +=1
			
		print "Got genotypes for " + str(single_snp.name)
	genotype_outfile = open("genotypes.data","wb")
	pickle.dump(genotypes,genotype_outfile)
	genotype_outfile.close()
	return genotypes


def genotype_descriptions(genotypes):					# get descriptions for relevant genotypes from snpedia
	try:
		description_infile = open("description.data","rb")		# again, try to cache data to save load...
	except IOError:
		descriptions = description_getter(genotypes)
	else:
		descriptions = pickle.load(description_infile)
		description_infile.close()
		
	return descriptions


def description_getter(genotypes):
	site = wiki.Wiki("http://snpedia.com/api.php")
	genotype_descriptions = {}
	for single_type in genotypes:
		for variant in genotypes[single_type]:
			genotype_name = str(single_type)+str(variant)
			print genotype_name
			wikipage = page.Page(site,genotype_name)
			if wikipage.exists == True:
				genotype_page = wikipage.getWikiText()
				if genotype_page.find("summary=") != -1:
					summary_start = genotype_page.find("summary=") + 8
					summary_stop = genotype_page.find("\n",summary_start)
					print genotype_page[summary_start:summary_stop]
					genotype_descriptions[genotype_name] = genotype_page[summary_start:summary_stop]
	description_outfile = open("description.data","wb")
	pickle.dump(genotype_descriptions,description_outfile)
	description_outfile.close()
	return genotype_descriptions
		


def genotype_comparer(my_filtered_snps,genotypes,descriptions):
	counter = 1
	bases = {'A': "T", "T" : "A", "G" : "C", "C" :"G", "I" : "I", "D" : "D"}
	candidates = {}	
	for single_snp in my_filtered_snps:	
																							# iterate over all SNPs from user 
		my_genotype = ""
		if genotypes.has_key(single_snp.name) and single_snp.genotype != "--":													# is this SNP called and are there genotypes for this SNP?  
			genotype = genotypes[single_snp.name]
			for single_genotype in genotype:																					# iterate over all known genotypes of a SNP
				if "("+single_snp.genotype[0]+";"+single_snp.genotype[1]+")" == single_genotype:						
					my_genotype = "("+single_snp.genotype[0]+";"+single_snp.genotype[1]+")"										# check for genotype of user (or its complementary bases)
				elif "("+bases[single_snp.genotype[0]]+";"+bases[single_snp.genotype[1]]+")" == single_genotype:
					my_genotype = "("+bases[single_snp.genotype[0]]+";"+bases[single_snp.genotype[1]]+")"
	
			if my_genotype != "":																								# if genotype of user is known: go on!
				output = "-----------\nMy genotype at " + single_snp.name+": "+single_snp.genotype+" or "+my_genotype+"\n"
				
				for single_genotype in genotype:	
					genotype_name = single_snp.name + single_genotype
					if descriptions.has_key(genotype_name):
						output = output + genotype_name +" is associated with: "+ descriptions[genotype_name]+"\n"
				if output.find("associated") != -1:
					print counter
					print output
					candidates[single_snp.name] = single_snp
					
					counter += 1
	return candidates


def freq_counter(candidates):
	other_genomes = glob("other_genomes/*.txt")
	variance = {}
	
	for single_genome in other_genomes:
		second_genome = {}
		second_genome = snp.reader_dict(single_genome)
		for my_snp in candidates:
			if second_genome.has_key(my_snp):
				if candidates[my_snp].other_genotypes.has_key(second_genome[my_snp].genotype):
					candidates[my_snp].other_genotypes[second_genome[my_snp].genotype] += 1
				else:
					candidates[my_snp].other_genotypes[second_genome[my_snp].genotype] = 1
				print candidates[my_snp].other_genotypes
	for single_snp in candidates:
		if variance.has_key(len(candidates[single_snp].other_genotypes)):
			variance[len(candidates[single_snp].other_genotypes)] += 1
		else:
			variance[len(candidates[single_snp].other_genotypes)] = 1
		if len(candidates[single_snp].other_genotypes) > 3:
			print single_snp + ": " + str(candidates[single_snp].other_genotypes)
	print variance
		


def main():
	print "Start getting SNPs from snpedia.com"
	snpedia_snps = snpedia_items()																				# get snps as hash: the snp-names are saved as keys
	print "Got " +str(len(snpedia_snps)) + " SNPs from snpedia.com\n"
	print "Start getting SNPs from user"
	my_snps = snp.reader("genome_Bastian_Greshake_Full_20110503120911.txt")										# get snps of user, saved as a list
	print "Got " + str(len(my_snps)) + " SNPs from user\n"
	print "Start filtering SNPs"
	my_filtered_snps = snp_compare(snpedia_snps,my_snps)														# get homozygotous snps of user which are in snpedia, saved as list
	print "Got " +str(len(my_filtered_snps)) + " SNPs that are available on SNPedia & are homozygotous\n"
	print "Start getting genotypes"
	genotypes = snpedia_genotypes(my_filtered_snps)																# get different genotypes for each snp, saved as hash. key = snp, value = list of genotypes
	print "Got " +str(len(genotypes)) + " SNPs from SNPedia which have known genotypes\n"
	print "Start getting genotype descriptions"
	descriptions = genotype_descriptions(genotypes) 															# get different descriptions for each genotype, saved as hash. key = genotype name, value = description
	print "Got "+ str(len(descriptions)) +" genotypes with a sufficient description\n"
	print "Start getting relevant genotypes\n"
	candidates = genotype_comparer(my_filtered_snps,genotypes,descriptions)
	print "Got relevant genotypes\n"
	print "Start counting frequency of candidates\n"
	freq_counter(candidates)

		
main()