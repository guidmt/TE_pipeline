# Rscript filtering_3prime_from_blast_modern.R args[1] args[2] args[3] args[4] args[5] args[6]
# args[1] = blastn output table (with no header and default table columns + "qlen" & "slen")
# args[2] = minimum length of the match on the 3' portion of the transposable element ("length" column; integer value)
# args[3] = maximum distance between "qlen" and "qend" values (integer value)
# args[4] = minimum number of non-matching bases on the 3' flanking portion (sense:"slen"-"send"; antisense:"sstart"-"send"; integer value)
# args[5] = flanking length for future processing (with perl scripts)
# args[6] = output fasta files prefix (prefix + "_3prime_sense/antisense_pos.txt")
# da modificare il ciclo di produzione delle posizioni nel caso reale!
library(readr)

args <- commandArgs(trailingOnly = TRUE)
length_flanking<-as.numeric(args[5])

blastn_table_name=as.character(args[1])
number_of_rows_blastn<-system(paste(paste('wc -l', blastn_table_name),'|cut -f 1'),intern=T)
number_of_rows_blastn<-as.numeric(sapply(strsplit(number_of_rows_blastn,split=' '),'[[',1))
size_chunks_to_create<-round(number_of_rows_blastn/1000)

idx<-split(data.frame(1:number_of_rows_blastn), rep(1:ceiling(number_of_rows_blastn/size_chunks_to_create), each=size_chunks_to_create, length.out=number_of_rows_blastn))

for(i in 1:length(idx)){

print(i)

as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
   
if(i==1){
   table_blastn<-read_lines(blastn_table_name,skip=min(idx[[i]])-1,n_max=max(idx[[i]])-min(idx[[i]]))
}else{
   table_blastn<-read_lines(blastn_table_name,skip=min(idx[[i]]),n_max=max(idx[[i]])-min(idx[[i]]))
}

table_blastn<-data.frame(do.call(rbind,strsplit(table_blastn,split='\t')))
table_blastn[,-c(1:2)]<-sapply( table_blastn[,-c(1:2)],as.numeric.factor)
table_blastn[,c(1:2)]<-sapply( table_blastn[,(1:2)],as.character.factor)
   
mmele=as.numeric(args[2])
maxdist=as.numeric(args[3])
mmfla=as.numeric(args[4])
flalen=as.numeric(args[5])
print(paste(Sys.time()," -> finished reading the arguments",sep=""))

tab_2=table_blastn[table_blastn[,4]>=mmele,]
tab_sense=tab_2[(tab_2[,10]-tab_2[,9])>0,]
tab_sense_2=tab_sense[(tab_sense[,13]-tab_sense[,8])<=maxdist,]
tab_sense_filt=tab_sense_2[(tab_sense_2[,14]-tab_sense_2[,10])>=mmfla,]
tab_antisense=tab_2[(tab_2[,10]-tab_2[,9])<0,]
tab_antisense_2=tab_antisense[(tab_antisense[,13]-tab_antisense[,8])<=maxdist,]
tab_antisense_filt=tab_antisense_2[tab_antisense_2[,10]>=mmfla,]
tab_sense_filt=tab_sense_filt[order(tab_sense_filt[,3],decreasing=TRUE),]
tab_antisense_filt=tab_antisense_filt[order(tab_antisense_filt[,3],decreasing=TRUE),]

print(paste(Sys.time()," -> finished filtering the blastn output tab",sep=""))

apply(tab_sense_filt,1,FUN=function(X,...){

   pos_sense=NULL

   riga=paste(X[2],(as.numeric(X[10])+1),length_flanking,sep="\t")

   print("start selection")

   output_string<-args[6]

   if (!is.element(riga,pos_sense[,1])) {

      print(!is.element(riga,pos_sense[,1]))

      write(riga,file=paste(output_string,"_3prime_sense_pos.txt",sep=""),append=T)

      elen=as.numeric(X[4])
      numeroterzo=(as.numeric(X[10])+1)-elen-flalen
      numeroquarto=2*(flalen)
	
      riga_fs=paste(X[2],numeroterzo,numeroquarto,sep="\t")
      write(riga_fs,file=paste(output_string,"_filled_sites_sense_pos.txt",sep=""),append=T)

      riga_5p=paste(tab_sense_filt[i,2],numeroterzo,length_flanking,sep="\t")
      write(riga_5p,file=paste(output_string,"_5prime_sense_pos.txt",sep=""),append=T)

   }
}
)

apply(tab_antisense_filt,1,FUN=function(X,...){
   
   pos_antisense=NULL
   
   riga=paste(X[2],(as.numeric(X[10])-1),length_flanking,sep="\t")
   
   print("start selection")
   
   output_string<-args[6]
   
   if (!is.element(riga,pos_antisense[,1])) {
      
      write(riga,file=paste(output_string,"_3prime_antisense_pos.txt",sep=""),append=T)
      
      elen=as.numeric(X[4])
      numeroterzo=(as.numeric(X[10])-1)+elen+flalen
      numeroquarto=2*(flalen)
      
      riga_fs=paste(X[2],numeroterzo,numeroquarto,sep="\t")
      write(riga_fs,file=paste(output_string,"_filled_sites_antisense_pos.txt",sep=""),append=T)
      
      riga_5p=paste(tab_sense_filt[i,2],numeroterzo,length_flanking,sep="\t")
      write(riga_5p,file=paste(output_string,"_5prime_antisense_pos.txt",sep=""),append=T)
   
   }
}
)

print(paste(Sys.time()," -> finished producing all the pos tabs!",sep=""))

}
