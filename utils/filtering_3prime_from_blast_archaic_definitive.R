# Rscript filtering_3prime_from_blast_archaic.R args[1] args[2] args[3] args[4] args[5] args[6]
# args[1] = blastn output table (with no header and default table columns + "qlen" & "slen")
# args[2] = reads database (in fasta format)
# args[3] = minimum length of the match on the 3' portion of the transposable element ("length" column; integer value)
# args[4] = maximum distance between "qlen" and "qend" values (integer value)
# args[5] = minimum number of non-matching bases on the 3' flanking portion (sense:"slen"-"send"; antisense:"sstart"-"send"; integer value)
# args[6] = output fasta files prefix (prefix + "_3prime(_flankings)_sense/antisense.fa")

args <- commandArgs(trailingOnly = TRUE)
library(data.table)
print(paste(Sys.time()," -> running script filtering_3prime_from_blast_archaic_definitive.R on file ",args[2],sep=""))

#
# Get the parameters for filtering
#

mmele=as.numeric(args[3])
maxdist=as.numeric(args[4])
mmfla=as.numeric(args[5])

#
# Read blastn file
#

blastn_table_name=as.character(args[1])
number_of_rows_blastn<-system(paste(paste('wc -l', blastn_table_name),'|cut -f 1'),intern=T)
number_of_rows_blastn<-as.numeric(sapply(strsplit(number_of_rows_blastn,split=' '),'[[',1))
size_chunks_to_create<-round(number_of_rows_blastn/10000)

idx<-split(data.frame(1:number_of_rows_blastn), rep(1:ceiling(number_of_rows_blastn/size_chunks_to_create), each=size_chunks_to_create, length.out=number_of_rows_blastn))

for(i in 1:length(idx)){

   print(i) 
   
   as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
   
   blastn_table<-read_lines(blastn_table_name,skip=min(idx[[i]]),n_max=max(idx[[i]]))
   blastn_table<- do.call(rbind,lapply(X=strsplit(blastn_table,split='\t'),FUN=function(X){as.data.frame(t(X))}))
   blastn_table[,-c(1:2)]<-sapply( blastn_table[,-c(1:2)],as.numeric.factor)
   blastn_table[,c(1:2)]<-sapply( blastn_table[,(1:2)],as.character.factor)
   
   
   # Filtering using the quality control
   paste(Sys.time()," -> finihed reading table ",args[1],sep="")
   
   blastn_table[,2]=paste(">",blastn_table[,2],sep="")
   
   tab_2=blastn_table[blastn_table[,4]>=mmele,]
   
   tab_sense=tab_2[which((tab_2[,10]-tab_2[,9])>0),]
   tab_sense_2=tab_sense[which((tab_sense[,13]-tab_sense[,8])<=maxdist),]
   tab_sense_filt=tab_sense_2[which(c(tab_sense_2[,14]-tab_sense_2[,10])>=mmfla),]
   
   tab_antisense=tab_2[which((tab_2[,10]-tab_2[,9])<0),]
   tab_antisense_2=tab_antisense[which((tab_sense[,13]-tab_sense[,8])<=maxdist),]
   tab_antisense_3=tab_antisense_2[abs(tab_antisense_2[,14]-tab_antisense_2[,9])<=1,]
   tab_antisense_filt=tab_antisense_3[tab_antisense_3[,10]>=mmfla,]
   

   if(nrow(tab_sense_filt)!=0 | nrow(tab_antisense_filt)!=0){
   
   reads_file_name=as.character(args[2])
   number_of_rows_reads<-system(paste(paste('wc -l', reads_file_name),'|cut -f 1'),intern=T)
   number_of_rows_reads<-as.numeric(sapply(strsplit(number_of_rows_reads,split=' '),'[[',1))

   n <- 10000
   nr <- number_of_rows_reads
   idx<-split(data.frame(1:number_of_rows_reads), rep(1:ceiling(nr/n), each=n, length.out=nr))
   
   idx2<-data.frame(do.call(rbind,lapply(idx,range)))
   
   fastasense=file(paste(args[6],"_3prime_sense.fa",sep=""),"w")
   fastantisense=file(paste(args[6],"_3prime_antisense.fa",sep=""),"w")
   fastaflankingsense=file(paste(args[6],"_3prime_flanking_sense.fa",sep=""),"w")
   fastaflankingantisense=file(paste(args[6],"_3prime_flanking_antisense.fa",sep=""),"w")


   }

   
   #
   # If the blast table with the sense reads is not empty
   #
   
   if(nrow(tab_sense_filt)!=0){
      
   match_sense<-apply(idx2,1,FUN=function(X,...){
      
      
   temp_reads<-read_lines(reads_file_name,skip=idx2[i,1],n_max=idx2[i,2])
   
   length_reads=length(temp_reads)
   
   reads_int=temp_reads[seq(1,length_reads-1,by=2)]
   reads_seq=temp_reads[seq(2,length_reads,by=2)]
   
   reads_tab=data.frame(reads_int,reads_seq,stringsAsFactors=FALSE,comment.char="")
   
   reads_tab_sense_filt=reads_tab[is.element(reads_tab[,2],tab_sense_filt[,2]),]
   
   tab_sense_filt=tab_sense_filt[is.element(tab_sense_filt[,2],reads_tab_sense_filt[,2]),]
   
   reads_tab_sense_filt=reads_tab_sense_filt[match(tab_sense_filt[,2],reads_tab_sense_filt[,2]),]
   
   writeLines(as.character(reads_tab_sense_filt[,1]),con=fastasense)
   
   # fasta_flanking_sense[j,1]=as.character(paste(">",reads_tab_sense_filt[i,1],sep=""))
   writeLines(reads_tab_sense_filt[,1],con=fastaflankingsense)
   
   writeLines(as.character(reads_tab[,2]),con=fastasense)
   
   asd=unlist(strsplit(as.character(reads_tab[,2]),split=""))[(as.numeric(temp_reads[i,10])+1):as.numeric(temp_reads[i,14])]
   flanking_seq=paste(asd,sep="",collapse="")
   
   # fasta_flanking_sense[j,1]=flanking_seq
   writeLines(as.character(flanking_seq),con=fastaflankingsense)
   
         }
      )
   
   
   } #end if sense

   
   #
   # If the blast table with the antisense reads is not empty
   #
   
   if(nrow(tab_antisense_filt)!=0){
      
   reads_antisense<-tab_antisense_filt[,2]

   match_antisense<-apply(idx2,1,FUN=function(X,...){
      
      temp_reads<-read_lines(reads_file_name,skip=idx2[i,1],n_max=idx2[i,2])
      
      length_reads=length(temp_reads)
      
      reads_int=temp_reads[seq(1,length_reads-1,by=2)]
      
      reads_seq=temp_reads[seq(2,length_reads,by=2)]
      
      reads_tab=data.frame(reads_int,reads_seq,stringsAsFactors=FALSE,comment.char="")
      
      reads_tab_antisense_filt=reads_tab[is.element(reads_tab[,1],tab_antisense_filt[,2]),]
      
      tab_antisense_filt=tab_antisense_filt[is.element(tab_antisense_filt[,2],reads_tab_antisense_filt[,2]),]
      
      reads_tab_antisense_filt=reads_tab_antisense_filt[match(tab_antisense_filt[,2],reads_tab_antisense_filt[,2]),]
      
      writeLines(as.character(reads_tab_antisense_filt[,1]),con=fastantisense)

      writeLines(as.character(reads_tab_antisense_filt[,1]),con=fastaflankingantisense)

      writeLines(as.character(reads_tab_antisense_filt[,2]),con=fastantisense)
      asd=unlist(strsplit(as.character(reads_tab_antisense_filt[,2]),split=""))[1:(as.numeric(tab_antisense_filt[i,10])-1)]
      flanking_seq=paste(asd,sep="",collapse="")

      writeLines(as.character(flanking_seq),con=fastaflankingantisense)
      
   }
   )
   
   }#end if antisense


}
