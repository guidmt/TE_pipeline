library(data.table)
# Rscript filtering_5prime_and_3prime_portions_of_modern_specific_empty_sites.R args[1] args[2] args[3] args[4]
# args[1] = empty sites file(in fasta format)
# args[2] = 5prime empty sites file (in fasta format)
# args[3] = 3prime empty sites file (in fasta format)
# args[4] = choose "sense" or "antisense"

args <- commandArgs(trailingOnly = TRUE)
library(data.table)
print(paste(Sys.time()," -> filtering_5prime_and_3prime_portions_of_modern_specific_empty_sites.R",sep=""))

empty_sites=fread(args[1],header=F,comment.char="",colClasses="character",blank.lines.skip=F)
sites_5p=fread(args[2],header=F,comment.char="",colClasses="character",blank.lines.skip=F)
sites_3p=fread(args[3],header=F,comment.char="",colClasses="character",blank.lines.skip=F)

sense<-as.character(args[4])
print(paste(Sys.time()," -> finished importing all the arguments",sep=""))

es_tab<-empty_sites[,c(1,4)]
tab_5p<-sites_5p[,c(1,4)]
tab_3p<-sites_3p[,c(1,4)]

es_tab=es_tab[(es_tab[,2]!=""),]
tab_5p=tab_5p[(tab_5p[,2]!=""),]
tab_3p=tab_3p[(tab_3p[,2]!=""),]

print(paste(Sys.time()," -> finished tabularizing the fasta files",sep=""))

# consider to use parApply for big data
splitting_table<-apply(DF_RES,1,FUN=function(X){
   
   splitone=unlist(strsplit(as.character(X[1]),split=":"))
   splittwo=unlist(strsplit(as.character(splitone[2]),split="-"))
   
   riga=data.frame(splitone[1],splittwo[1],splittwo[2],as.character(X[2]))
   
})

tab_3p_sep<-do.call(rbind,splitting_table(tab_3p))
print(paste(Sys.time()," -> finished splitting the 3prime table",sep=""))
colnames(tab_3p_sep)=c("V1","V2","V3","V4")

tab_5p_sep<-do.call(rbind,splitting_table(tab_5p))
print(paste(Sys.time()," -> finished splitting the 5prime table",sep=""))
colnames(tab_5p_sep)=c("V1","V2","V3","V4")

es_tab_sep<-do.call(rbind,splitting_table(es_tab))
print(paste(Sys.time()," -> finished splitting the empty sites table",sep=""))
colnames(es_tab_sep)=c("V1","V2","V3","V4")

livels_analysis=levels(as.factor(es_tab_sep[,1]))

es_tab_sep=data.table(es_tab_sep)
tab_5p_sep=data.table(tab_5p_sep)
tab_3p_sep=data.table(tab_3p_sep)

save.image(paste("test_filtering_ref_es_",sense,".RData",sep=""))

setkey(es_tab_sep,V1)
setkey(tab_5p_sep,V1)
setkey(tab_3p_sep,V1)

tab_5p_sep_filt=NULL
tab_3p_sep_filt=NULL
es_tab_sep_filt=NULL

for (j in 1:length(livels_analysis)) {

   temp_es=es_tab_sep[.(as.character(livels_analysis[j]))]
   temp_5p=tab_5p_sep[.(as.character(livels_analysis[j]))]
   temp_3p=tab_3p_sep[.(as.character(livels_analysis[j]))]
   
   temp_es=as.data.frame(temp_es,stringsAsFactors=FALSE)
   temp_5p=as.data.frame(temp_5p,stringsAsFactors=FALSE)
   temp_3p=as.data.frame(temp_3p,stringsAsFactors=FALSE)
   
   if (sense=="sense") {
      
      temp_es_filt=temp_es[((is.element(temp_es[,3],temp_3p[,3]))&(is.element(temp_es[,2],temp_5p[,2]))),]
      temp_3p_filt=temp_3p[is.element(temp_3p[,3],temp_es_filt[,3]),]
      temp_5p_filt=temp_5p[is.element(temp_5p[,2],temp_es_filt[,2]),]
      
   }
   if (sense=="antisense") {
      
      temp_es_filt=temp_es[((is.element(temp_es[,2],temp_3p[,2]))&(is.element(temp_es[,3],temp_5p[,3]))),]
      temp_3p_filt=temp_3p[is.element(temp_3p[,2],temp_es_filt[,2]),]
      temp_5p_filt=temp_5p[is.element(temp_5p[,3],temp_es_filt[,3]),]
      
   }
   
   temp_es_filt=temp_es_filt[order(as.numeric(temp_es_filt[,2]),as.numeric(temp_es_filt[,3]),decreasing=F),]
   temp_3p_filt=temp_3p_filt[order(as.numeric(temp_3p_filt[,2]),as.numeric(temp_3p_filt[,3]),decreasing=F),]
   temp_5p_filt=temp_5p_filt[order(as.numeric(temp_5p_filt[,2]),as.numeric(temp_5p_filt[,3]),decreasing=F),]
   es_tab_sep_filt=rbind(es_tab_sep_filt,temp_es_filt)
   tab_3p_sep_filt=rbind(tab_3p_sep_filt,temp_3p_filt)
   tab_5p_sep_filt=rbind(tab_5p_sep_filt,temp_5p_filt)
   
}
print(paste(Sys.time()," -> finished splitting, sorting and filtering the tables, now producing fasta files",sep=""))
fasta_es_filt=file(paste("reference_sorted_empty_sites_for_archaic_insertions_",sense,".fa",sep=""),"w")
fasta_5p_filt=file(paste("reference_sorted_empty_sites_5prime_portion_",sense,".fa",sep=""),"w")
fasta_3p_filt=file(paste("reference_sorted_empty_sites_3prime_portion_",sense,".fa",sep=""),"w")

for (i in 1:nrow(es_tab_sep_filt)) {
   writeLines(paste(">",as.character(es_tab_sep_filt[i,1]),":",as.character(es_tab_sep_filt[i,2]),"-",as.character(es_tab_sep_filt[i,3]),"_emptysite",i,sep=""),con=fasta_es_filt)
   writeLines(as.character(es_tab_sep_filt[i,4]),con=fasta_es_filt)
}

for (i in 1:nrow(tab_5p_sep_filt)) {
     writeLines(paste(">",as.character(tab_5p_sep_filt[i,1]),":",as.character(tab_5p_sep_filt[i,2]),"-",as.character(tab_5p_sep_filt[i,3]),"_emptysite",i,sep=""),con=fasta_5p_filt)
     writeLines(as.character(tab_5p_sep_filt[i,4]),con=fasta_5p_filt)
}

for (i in 1:nrow(tab_3p_sep_filt)) {
     writeLines(paste(">",as.character(tab_3p_sep_filt[i,1]),":",as.character(tab_3p_sep_filt[i,2]),"-",as.character(tab_3p_sep_filt[i,3]),"_emptysite",i,sep=""),con=fasta_3p_filt)
     writeLines(as.character(tab_3p_sep_filt[i,4]),con=fasta_3p_filt)
}

close(fasta_es_filt)
close(fasta_5p_filt)
close(fasta_3p_filt)

print(paste(Sys.time()," -> DONE!",sep=""))

