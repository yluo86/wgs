### This R script outputs a given signal region, assuming it has been annotated using VEP and has gene locations using Biosmart

	genes<-read.table("out-gene.txt",h=T,stringsAsFactors=FALSE)
	vep<-read.table("out-vep.txt",h=T,stringsAsFactors=FALSE)
	vep$Start<-(genes[match(vep$Gene,genes$GENE),]$START)
	vep$End<-(genes[match(vep$Gene,genes$GENE),]$END)
	vep$Chr<-(genes[match(vep$Gene,genes$GENE),]$CHR)

	rsids<-as.character(unique(vep$Uploaded_variation))
	regions<-function(x,window=5000){
		if (sum(!is.na(vep[vep[,1]==x,]$Start))==0 & sum(!is.na(vep[vep[,1]==x,]$End))==0){
			bp<-unlist(strsplit(as.character(vep[vep[,1]==x,]$Location),":"))[2]
			chr<-unlist(strsplit(as.character(vep[vep[,1]==x,]$Location),":"))[1]
			
			if(length(unlist(strsplit(bp,"-")))==2){
				start<-as.numeric(unlist(strsplit(bp,"-"))[1])-window
				end<-as.numeric(unlist(strsplit(bp,"-"))[2])+window
			}
			else{
				start<-as.numeric(bp)-window
				end<-as.numeric(bp)+window
			}
		}
		else{
			chr<-unique(vep[vep[,1]==rsids[3] & !is.na(vep$Chr),]$Chr)
			start<-min(vep[vep[,1]==x,]$Start,na.rm=T)
			end<-min(vep[vep[,1]==x,]$End,na.rm=T)
		}
		return(c(x,chr,start,end))
	} 
	
	output<-matrix(sapply(rsids,regions),ncol=4,byrow=T)

	test<-subset(vep,select=c("Uploaded_variation","SYMBOL"))
	genename<-aggregate(test[2],test[-2],FUN = function(X) paste(unique(X), collapse=","))
	output<-cbind(output,genename[match(output[,1],genename$Uploaded_variation),2])
	write.table(output,file="out-signal-regions.txt",quote=FALSE,sep="\t",col.names=c("rsID","Chr","Start","End","Gene"),row.names=FALSE)

