#####locating variants that are in the known regions####
#R CMD BATCH --no-restore '--args sigfile="gene-file.txt" window=5000 output="known-sig.txt"' search-known-sigs.R

args=(commandArgs(TRUE)) 
 if(length(args)==0){
	     print("No filename arguments supplied.")
 }else{
	      for(i in 1:length(args)){
		                eval(parse(text=args[[i]]))
      }
  }

sigs<-NULL 
regions<-read.table(sigfile,h=T,stringsAsFactors=FALSE)

for ( chrom in unique(regions$Chr) ){
		
		print(chrom)

	load(paste("chrom",chrom,"/chrom",chrom,"-snptest-assoc.RData",sep="")) #resource file which has BP info
		
		region<-regions[regions$Chr==chrom,]
		sigs<-c(sigs,as.character(data[sapply(data$BP,function(x) { any( (region$Start-window) <x & x< (region$End+window) )}),]$SNP))	
			rm(data)
}

write.table(sigs,file=output,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)	
