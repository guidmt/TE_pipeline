# Rscript generate_longer_reads_from_blasted_db_archaic.R args[1] args[2] args[2] args[3] args[4] args[5] args[6] args[8]
# args[1] = blastn output table prefix (with no header and default table columns + "qlen" and "slen", .txt format)
# args[2] = choose "3prime" or "5prime"
# args[3] = choose "sense" or "antisense"
# args[4] = choose "flankings" or "reads" (whatever argument that is not "flankings" is recognized as "reads")
# args[5] = minimum match to consider in the blastn output table
# args[6] = output file in fasta format
# args[7] = choose "partial" or "final"

args <- commandArgs(trailingOnly = TRUE)

print(paste(Sys.time()," -> running script generate_longer_reads_from_blasted_db_archaic.R on file ",args[2],sep=""))

side=as.character(args[2])

sense=as.character(args[3])

type=as.character(args[4])

mmatch=as.numeric(args[5])

operation="final"

number_of_rows_blastn<-system(paste(paste('wc -l', paste(args[1],".txt",sep="")),'|cut -f 1'),intern=T)
number_of_rows_blastn<-as.numeric(sapply(strsplit(number_of_rows_blastn,split=' '),'[[',1))
size_chunks_to_create<-round(number_of_rows_blastn/50000)

idx<-split(data.frame(1:number_of_rows_blastn), rep(1:ceiling(number_of_rows_blastn/size_chunks_to_create), each=size_chunks_to_create, length.out=number_of_rows_blastn))

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
   

   fasta_out=file(args[6],"w")
   
   temp_tab=NULL
   empty_sites=NULL
   fl5p=NULL

   if (type=="flankings") {
      
      blast_tab=table_blastn_df[!duplicated(table_blastn_df[,1]),]
      
      vec=as.factor(paste(">",blast_tab[,2],":",blast_tab[,9],sep=""))
      
      blast_tab=cbind(blast_tab,vec)
      
      lev=levels(blast_tab[,15])
   } 
   
   
   if (operation=="final") {
   
   for (i in 1:length(lev)) {
   
         temp_tab=blast_tab[blast_tab[,15]==lev[i],]
                     
         temp_tab=temp_tab[temp_tab[,4]>=mmatch,]
         
         diff=as.numeric(temp_tab[,10])-as.numeric(temp_tab[,9])+1
         
         temp_tab=cbind(temp_tab,diff)
         
         diff_sort=sort(temp_tab$diff,decreasing=TRUE)
         
         temp_tab_high=temp_tab[temp_tab$diff==diff_sort[1],]
         
         temp_tab_low=temp_tab[temp_tab$diff==diff_sort[length(diff_sort)],]
         
         temp_tab_high=t(temp_tab_high[1,])
         
         temp_tab_low=t(temp_tab_low[1,])
         
         if (side=="3prime"){ reads_tab_high=temp_tab[temp_tab[,ncol(temp_tab)-1]==temp_tab_high[1],] }
         
         if (side=="3prime"){ reads_tab_low=temp_tab[temp_tab[,ncol(temp_tab)-1]==temp_tab_low[1],] }
         
         reads_tab_high=t(reads_tab_high[1,])
         
         reads_tab_low=t(reads_tab_low[1,])
         
         #
         # Qui e' dove vengono prodotti gli output che ci interessano
         #
         
         if (type=="flankings") {
      
            writeLines(paste(">site",i,"_",side,"_",sense,sep=""),con=fasta_out)
            writeLines(reads_tab_high[2],con=fasta_out)
            if (sense=="sense") posizione=as.numeric(temp_tab_high[9])
            if (sense=="antisense") posizione=as.numeric(temp_tab_high[10])
            riga_5pf=paste(temp_tab_high[2],posizione,sep="\t")
            fl5p=rbind(fl5p,riga_5pf)
            posizione_es=posizione-100
            riga_es=paste(temp_tab_high[2],posizione_es,sep="\t")
            empty_sites=rbind(empty_sites,riga_es)
         }
      
      }
         
      }
      
      close(fasta_out)
      if (type=="flankings") {
         if (sense=="sense") {
            write.table(empty_sites,"empy_sites_sense_pos_list.txt",col.names=F,row.names=F,quote=F)
            write.table(fl5p,"5prime_pos_sense_list.txt",col.names=F,row.names=F,quote=F)
         }
         if (sense=="antisense") {
            write.table(empty_sites,"empy_sites_antisense_pos_list.txt",col.names=F,row.names=F,quote=F)
            write.table(fl5p,"5prime_pos_antisense_list.txt",col.names=F,row.names=F,quote=F)
         }
      }

      
}

print(paste(Sys.time()," -> ",args[2]," is DONE!",sep=""))