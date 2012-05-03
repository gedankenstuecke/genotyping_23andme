import sys,os,pickle,snp,datetime, PyRSS2Gen
from Bio import Entrez
from glob import glob

# This script crawls PubMed for all SNPs included in a 23andme-raw-data-file 
# and publishes the resulting hits as a RSS-feed. The script is meant to run
# once a month using cron etc. 
# New items will be added to the RSS-feed while old items are kept. 
# To reduce load on the NCBI-servers, old items are stored locally

# Dependencies: 
# BioPython http://www.biopython.org/
# PyRSS2Gen http://www.dalkescientific.com/Python/PyRSS2Gen.html
# The SNP class found at https://github.com/gedankenstuecke/genotyping_23andme


# CONFIG, EDIT FOR YOUR NEEDS! # 
Entrez.email  = "bgreshake@googlemail.com"	                        # Your email, so NCBI can contact you, if needed
raw_data      = "genome_Bastian_Greshake_Full_20110503120911.txt"   # The genotyping-raw-data-file of interest
last_update   = "last.update"					                    # The file where the last-update-time will be saved
old_items     = "pubmed_rss.backup"				                    # Were old feed-times shall be stored
feed_location = "23andme_feed.xml"			                        # Were the RSS-Feed shall be published (/var/www/younameit)
pickleFile    = 'pubmed_crawler.pickle'

def store(data, file):
	f = open(file, "wb")
	try:
		pickle.dump(data, f)
	except:
		print "error saving"
		sys.exit(1)
	f.close()

def load(file):
	if os.path.exists(file):
		f = open(file, 'rb')
		data = pickle.load(f)
		f.close()
		return data
	return None

def paper_getter(snps,last_update, startAt=1):
	new_hits = []
	if glob(last_update) == []:										#check for time of last update
		today      = datetime.date.today().strftime("%Y %b")
		updated_on = datetime.datetime.strptime(today, "%Y %b")
		pickle_out = open(last_update,"wb")
		pickle.dump(updated_on,pickle_out)
		pickle_out.close()										#create new "last update" if needed
		updated_on = ""
	else:													#else load time of last update
		pickle_in  = open(last_update,"rb")
		updated_on = pickle.load(pickle_in)
		
	
	if os.path.exists(pickleFile):
		data     = load(pickleFile)
		if startAt == 1: startAt  = data[0]
		new_hits = data[1]
		print "RE STARTING AT POSITION", startAt, "WITH",len(new_hits), "HITS"
	
	snp_counter = 1	
	for single_snp in snps:											#iterate over all SNPs
		try:
			snp_counter += 1
			if snp_counter <= startAt:
				print "skipping snp", snp_counter, 'till',startAt
				continue

			handle  = Entrez.esearch(db="pubmed", term=single_snp.name)	
			results = Entrez.read(handle)								#Get all papers for the SNP
			print snp_counter,

			if results["IdList"] != []:								#If papers are found for SNP:
				for single_paper in results["IdList"]:
					try:
						summary_handle = Entrez.esummary(db="pubmed", id=single_paper)	 #Get summary
						summary = Entrez.read(summary_handle)
						if updated_on == "" or updated_on < datetime.datetime.strptime(summary[0]["PubDate"][:8], "%Y %b"): #If new: Get further information
							paper         = Papers(single_paper)					# Create new paper with basic information about it
							paper.snp     = single_snp.name
							paper.title   = summary[0]["Title"]
							paper.author  = summary[0]["AuthorList"][0]
							paper.journal = summary[0]["Source"]
							try:
								detail_handle = Entrez.efetch(db="pubmed",id=single_paper,retmode="xml")
								details = Entrez.read(detail_handle)
								try:
									paper.abstract = str(details[0]["MedlineCitation"]["Article"]["Abstract"]["AbstractText"])
								except:
									paper.abstract = "no abstract available"
								print paper.snp + ": " + paper.title
								print "Written by "    + paper.author + " et al"
								print "Published in "  + paper.journal
								print paper.abstract
								new_hits.append(paper)
							except:
								print "NO UPDATE"
								pass
					except:
						print "NO SUMMARY"
						pass
			else:
				print "NONE"
			data     = [ snp_counter, new_hits ]
			store(data, pickleFile)
		except:
			pass
	if os.path.exists(pickleFile):
		os.remove(pickleFile)
	return new_hits


def item_creator(new_paper):				# Create RSS-items out of all new papers
	items = []
	for paper in new_paper:
		heading = paper.snp + ": " + paper.title
		item = PyRSS2Gen.RSSItem(title = heading)
		description = "Written by " + paper.author + "\n"
		description = description + "Published in " + paper.journal + "\n"
		description = description + "Abstract: " + paper.abstract + "\n"
		description = description + "PMID: " + paper.pmid
		link = "http://www.ncbi.nlm.nih.gov/pubmed?term=" + paper.pmid
		item.description = description
		item.link = link
		item.guid = PyRSS2Gen.Guid(link)
		items.append(item)
	return items


def item_import(rss_items,backup):			# Import old RSS-Items
	if glob(backup) != []:
		pickle_in = open(backup,"rb")
		old_rss_items = pickle.load(pickle_in)
		for old_item in old_rss_items:
			rss_items.append(old_item)
		pickle_in.close()
		
		pickle_out = open(backup,"wb")
		pickle.dump(rss_items,pickle_out)
		pickle_out.close()
		
	else:
		pickle_out = open(backup,"wb")
		pickle.dump(rss_items,pickle_out)
		pickle_out.close()
	return rss_items


def create_feed(all_rss_items,feed_location):	# Create RSS-Feed
	rss = PyRSS2Gen.RSS2(
		title = "Papers on SNPs analyzed by 23andme",
		link = "https://github.com/gedankenstuecke/genotyping_23andme",
		description = "Monthly updated feed on all new papers listed in PubMed about SNPs which are analyzed by 23andme",
		lastBuildDate = datetime.datetime.utcnow(),
		items = all_rss_items,
		)
	rss.write_xml(open(feed_location,"w"))


class Papers():								# simple Class for storing basic information about papers
	"""docstring for paper"""
	def __init__(self, pmid):
		self.pmid = pmid
		self.title = ""
		self.abstract = ""
		self.author = ""
		self.snp =  ""
		self.journal = ""


if sys.argv[1] is not None:
	if os.path.exists(sys.argv[1]):
		raw_data      = sys.argv[1] # The genotyping-raw-data-file of interest

startPos = 1

if len(sys.argv) > 2 and sys.argv[2] is not None:
	if sys.argv[2].isdigit() and sys.argv[2] > 0:
		startPos = int(sys.argv[2])


sys.stderr.write("INPUT FILE " + raw_data + " STARTING AT " + str(startPos) + "\n")

sys.stderr.write("READING SNPS\n")
snps          = snp.reader(   raw_data)

sys.stderr.write("GETTING PAPER\n")
new_paper     = paper_getter( snps,         last_update, startPos)

sys.stderr.write("CREATING RSS\n")
rss_items     = item_creator( new_paper)

sys.stderr.write("IMPORT RSS\n")
all_rss_items = item_import(  rss_items,    old_items)

sys.stderr.write("CREATING FEED\n")
create_feed(                  all_rss_items,feed_location)

sys.stderr.write("FINISHED\n")

# todo: backup new items, import old items, create feed
