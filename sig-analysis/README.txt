####Categorising Signals ####

# 1. low-freq and rare variants in known IBD regions from 3 sources:
#	a). 163 loci from Jostins et al
#	b). 40 genes from Holm (3 overlaps with Jostins et al.)
#	c). 40 genes from TEAs

	## To find the gene locations, I applied the following procedure:
		# (1). For a given list of rsIDs (hits), using VEP to get the most related 
			# http://grch37.ensembl.org/Homo_sapiens/Tools/VEP/
			# out-vep.txt
					
		# (2). Then use BioMart to determine Gene locations:
			# http://grch37.ensembl.org/biomart/martview/
			# out-gene.txt
									
		# (3). After which using sig-locations.Rscript to get overlapping regions. For intergenic variants ±5kbp are taken for the given BP
