== Introduction ==

This repository includes the raw genotyping data, as obtained by 23andme, using their latest chip (as used in april 2011). The data was downloaded on May, 3rd 2011. Updates will be released as soon as the chips used by 23andme provide further information. 

Feel free to use this data, out of curiosity or for your research. I'd love to hear from people who use this data and are abled to get some nice results from it.
If you need some phenotypic data just let me know. You can reach me via email (bgreshake@googlemail.com) or phone (+49 176 213 044 66). 

Cheers,
Bastian Greshake
May, 3rd 2011

Last Revision of README: May, 26th 2011


== Content == 

* ''genome_Bastian_Greshake_Full_20110503120911.txt'' My genotyping raw data, as provided by 23andme
* ''genome_Bastian_Greshake_Full_20110503120911.html'' The Promethease-report, generated from my raw data
* ''snp.py'' A simple parser for the raw data, including a simple SNP-class
* ''snpedia.py'' A script which screens all homozygotous SNPs in the raw data and gets additional information about their effects from SNPedia 
* ''pubmed_crawler.py'' A script which crawls PubMed for all papers which concern SNPs which are analyzed by 23andme and publishes it's result as a RSS-feed

== snp.py ==
Contains a method to parse the genotyping-raw-data, SNPs are parsed into a class which contains information about the SNPs name, chromosome, position and genotype

== snpedia.py ==
The script gets information about each homozygotous SNP from SNPedia.com (as far as they have information) and gives back all homozygotous SNPs, which are annotated in the SNPedia. Including all known genotypes and their function. 

I wrote this script to get information about homozygotous, dominant SNPs i carry, to infer information about the phenotypes of my parents. The script returns candidate SNPs which have known genotypes with a sufficient description. 

Up to know these candidate SNPs have to be screened by hand, to see if I carry the dominant genotype (otherwise no information about parental phenotype can be won). To reduce load on SNPedia data on genotypes and descriptions are stored locally in *.data-files

=== Dependencies ===

* snp.py
* [http://code.google.com/p/python-wikitools/ Python-Wiki-Tools]

== pubmed_crawler.py ==
The script crawls the PubMed-database for scientific papers on the SNPs which are analyzed by 23and me. And republishes the results as a RSS-feed. Already crawled papers are stored locally to reduce load on the NCBI-servers. 

The script is written to be run once a month (using cron etc.) to look for new papers, which are then added to the existing RSS-feed.

Last-update-time is stored in last.update, old RSS-items are stored in pubmed_rss.backup, the feed is published to 23andme_feed.xml    

Whats missing: An option to split the job into chunks, otherwise a full run on 1 Mio SNPs will take ~4 days to complete.   

An example RSS-Feed can be found in 23andme_feed.xml

=== Dependencies ===

* snp.py
* [http://biopython.org/wiki/Main_Page BioPython]
* [http://www.dalkescientific.com/Python/PyRSS2Gen.html PyRSS2Gen] 
