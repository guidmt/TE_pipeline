# Rscript filtering_5prime_from_blast_archaic_v2.R args[1] args[2] args[3] args[4] args[5]
# args[1] = suffix of the filtered blastn output table (with no header and default table columns + "qlen" & "slen"; blasted_archaic_5prime_flanking_args[5]_args[1].txt for example)
# args[2] = reads database (in fasta format)
# args[3] = minimum number of non-matching bases
# args[4] = minimum lenght of the match
# args[5] = choose "sense" or "antisense", for the direction of the non-matching bases

library(data.table)

args <- commandArgs(trailingOnly = TRUE)

# 
# Read the parameters from the command line
# 
sense=as.character(args[5])
mmele=as.numeric(args[3])
mmfla=as.numeric(args[4])

print(paste(Sys.time()," -> executing filtering_5prime_from_blast_archaic_v2.R script on file ",args[2]," ",args[5],sep=""))

inputfile=paste("blasted_archaic_5prime_flanking_",sense,"_",as.character(args[1]),".txt",sep="")


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
  
    tab_filt_2<-blastn_table

    # 
    # Filtering
    # 

    tab_filt_2=tab_filt_2[tab_filt_2[,4]>=mmfla,]
    tab_filt_2=tab_filt_2[(tab_filt_2[,14])-(tab_filt_2[,4])>=mmele,]
    tab_filt_2[,2]=paste(">",tab_filt_2[,2],sep="")
    
    if(exists("tab_filt2")){
      
      if(nrow(tab_filt2)!=0){
        
      vector_to_use<-rep("fill",nrow(tab_sense_filt)*2)
        
      vector_to_use[seq(1,length(vector_to_use)-1,by=2)]<-tab_filt2[,2]
      vector_to_use[seq(2,length(vector_to_use),by=2)]<-tab_filt2[,ncol(tab_filt2)]
      
      fasta_out=file(paste("putative_archaic_specific_5prime_reads_",sense,"_",as.character(args[1]),".fa",sep=""),"w")
      
      writeLines(vector_to_use,con=fasta_out)

      }
    }
  }
        
        


