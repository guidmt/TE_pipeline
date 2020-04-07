library(readr)

# Rscript filtering_3prime_from_blast_archaic.R args[1] args[2] args[2] args[3] args[4] args[5]
# args[1] = blastn output table (with no header and default table columns + "qlen" & "slen")
# args[2] = minimum length of the match on the 3' portion of the transposable element ("length" column; integer value)
# args[3] = maximum distance between "qlen" and "qend" values (integer value)
# args[4] = minimum number of non-matching bases on the 3' flanking portion (sense:"slen"-"send"; antisense:"sstart"-"send"; integer value)
# args[5] = output fasta files prefix (prefix + "_3prime(_flankings)_sense/antisense.fa")

args <- commandArgs(trailingOnly = TRUE)
library(data.table)
print(paste(Sys.time()," -> running script filtering_3prime_from_blast_archaic_definitive.R on file ",args[2],sep=""))

#
# Get the parameters for filtering
#

mmele=as.numeric(args[2])
maxdist=as.numeric(args[3])
mmfla=as.numeric(args[4])

#
# Read blastn file
#

blastn_table_name=as.character(args[1])
number_of_rows_blastn<-system(paste(paste('wc -l', blastn_table_name),'|cut -f 1'),intern=T)
number_of_rows_blastn<-as.numeric(sapply(strsplit(number_of_rows_blastn,split=' '),'[[',1))
size_chunks_to_create<-round(number_of_rows_blastn/50000)

idx<-split(data.frame(1:number_of_rows_blastn), rep(1:ceiling(number_of_rows_blastn/size_chunks_to_create), each=size_chunks_to_create, length.out=number_of_rows_blastn))

fastasense=paste(args[5],"_3prime_sense.fa",sep="")
fastantisense=paste(args[5],"_3prime_antisense.fa",sep="")
fastaflankingsense=paste(args[5],"_3prime_flanking_sense.fa",sep="")
fastaflankingantisense=paste(args[5],"_3prime_flanking_antisense.fa",sep="")

for(i in 1:length(idx)){

   print(i) 
   
   as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
   
   if(i==1){
   blastn_table<-read_lines(blastn_table_name,skip=min(idx[[i]])-1,n_max=max(idx[[i]])-min(idx[[i]]))
   }else{
   blastn_table<-read_lines(blastn_table_name,skip=min(idx[[i]]),n_max=max(idx[[i]])-min(idx[[i]]))
   }
   
   table_blastn_df<- do.call(rbind,lapply(X=strsplit(blastn_table,split='\t'),FUN=function(X){as.data.frame(t(X))}))
   table_blastn_df<-data.frame(table_blastn_df[,2],table_blastn_df[,-2])
   
   print(nrow(table_blastn_df))
   
   table_blastn_df[,-c(1:2,ncol(table_blastn_df))]<-sapply( table_blastn_df[,-c(1:2,ncol(table_blastn_df))],as.numeric.factor)
   table_blastn_df[,c(1:2,ncol(table_blastn_df))]<-sapply( table_blastn_df[,c(1:2,ncol(table_blastn_df))],as.character.factor)
   
   table_blastn_df[,2]=paste(">",table_blastn_df[,2],sep="")
   
   tab_2=table_blastn_df[table_blastn_df[,4]>=mmele,]
   
   tab_sense=tab_2[which((tab_2[,10]-tab_2[,9])>0),]
   
   if(nrow(tab_sense)!=0){
      
   tab_sense_2=tab_sense[which((tab_sense[,13]-tab_sense[,8])<=maxdist),]
   tab_sense_filt=tab_sense_2[which(c(tab_sense_2[,14]-tab_sense_2[,10])>=mmfla),]
   
   }
   
   tab_antisense=tab_2[which((tab_2[,10]-tab_2[,9])<0),]
   
   if(nrow(tab_antisense)!=0){
      
   tab_antisense_2=tab_antisense[which((tab_sense[,13]-tab_sense[,8])<=maxdist),]
   tab_antisense_3=tab_antisense_2[which(abs(tab_antisense_2[,14]-tab_antisense_2[,9])<=1),]
   tab_antisense_filt=tab_antisense_3[which(tab_antisense_3[,10]>=mmfla),]
   
   }
   
   #
   # If the blast table with the sense reads is not empty
   #
   
   if(exists("tab_sense_filt")){
      if(nrow(tab_sense_filt)!=0){
      
   print("true sense")
      
   vector_sense<-rep("fill",nrow(tab_sense_filt)*2)
   
   vector_sense[seq(1,length(vector_sense)-1,by=2)]<-tab_sense_filt[,2]
   vector_sense[seq(2,length(vector_sense),by=2)]<-tab_sense_filt[,ncol(tab_sense_filt)]
      
   write.table(vector_sense,file=fastasense,append=T,quote=F,col.names=F,row.names=F)
   
   asd<-apply(tab_sense_filt,1,FUN=function(X){
      
      seqstring=unlist(strsplit(as.character(X[length(X)]),split=""))
      paste(seqstring[(as.numeric(X[10])+1):as.numeric(X[14])],collapse="")
      
   })
   
   vector_sense_flank<-rep("fill",nrow(tab_sense_filt)*2)
   
   vector_sense_flank[seq(1,length(vector_sense_flank)-1,by=2)]<-tab_sense_filt[,2]
   vector_sense_flank[seq(2,length(vector_sense_flank),by=2)]<-asd
   
   
   write.table(vector_sense_flank,file=fastaflankingsense,append=T,quote=F,col.names=F,row.names=F)
   
      }
   }
   
   #
   # If the blast table with the antisense reads is not empty
   #
   
   if(exists("tab_antisense_filt")) {
      if(nrow(tab_antisense_filt)!=0){
      
      print("true antisense")
      
      vector_antisense<-rep("fill",nrow(tab_antisense_filt)*2)
      
      vector_antisense[seq(1,length(vector_antisense)-1,by=2)]<-tab_antisense_filt[,2]
      vector_antisense[seq(2,length(vector_antisense),by=2)]<-tab_antisense_filt[,ncol(tab_antisense_filt)]
      
      write.table(vector_antisense,file=fastantisense,append=T,quote=F,col.names=F,row.names=F)
      
      asd<-apply(tab_antisense_filt,1,FUN=function(X){
         
         seqstring=unlist(strsplit(as.character(X[length(X)]),split=""))
         paste(seqstring[1:(as.numeric(X[10])-1)],collapse="")
         
      })
      
      vector_antisense_flank<-rep("fill",nrow(tab_antisense_filt)*2)
      
      vector_antisense_flank[seq(1,length(vector_antisense_flank)-1,by=2)]<-tab_antisense_filt[,2]
      vector_antisense_flank[seq(2,length(vector_antisense_flank),by=2)]<-asd
      
      print(asd)
      
      write.table(vector_antisense_flank,file=fastaflankingantisense,append=T,quote=F,col.names=F,row.names=F)
   
   }#end if antisense
   }
   

}
