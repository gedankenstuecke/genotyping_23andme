import sys

def reader(infile):
	handle = open(infile,"r")
	data = []
	for x in handle:
		if x[0] != "#":
			single = x.split("\t")
			data.append(snp(single[0],single[1],single[2],single[3].rstrip()))
	handle.close()
	return data	

def reader_dict(infile):
	handle = open(infile,"r")
	data = {}
	for x in handle:
		if x[0] != "#":
			single = x.split("\t")
			data[single[0]] = snp(single[0],single[1],single[2],single[3].rstrip())
	handle.close()
	return data


class snp():
	def __init__(self,name,chromosome,position,genotype):
		self.name = name
		self.chromosome = chromosome
		self.position = position
		self.genotype = genotype
		self.other_genotypes = {}

#reader(sys.argv[1])