# Rscript generate_longer_reads_from_blasted_db_archaic.R args[1] args[2] args[3] args[4] args[5] args[6] args[7] args[8]
# args[1] = blastn output table prefix (with no header and default table columns + "qlen" and "slen", .txt format)
# args[2] = filtered reads database (in fasta format)
# args[3] = choose "3prime" or "5prime"
# args[4] = choose "sense" or "antisense"
# args[5] = choose "flankings" or "reads" (whatever argument that is not "flankings" is recognized as "reads")
# args[6] = minimum match to consider in the blastn output table
# args[7] = output file in fasta format
# args[8] = choose "partial" or "final"

args <- commandArgs(trailingOnly = TRUE)

print(paste(Sys.time()," -> running script generate_longer_reads_from_blasted_db_archaic.R on file ",args[2],sep=""))

#
# Reads the parameters for analysis
#

side=as.character(args[3])
senso=as.character(args[4])
tipo=as.character(args[5])
mmatch=as.numeric(args[6])
operazione=as.character(args[8])

blast_tab=read.table(paste(args[1],".txt",sep=""),header=F,comment.char="",colClasses=c("character","character","numeric","integer","integer","integer","integer","integer","integer","integer","numeric","numeric","integer","integer"))

reads_tab=read.table(args[2],header=F,comment.char="",colClasses="character")

print(paste(Sys.time()," -> finished importing and tabularizing all the arguments for file ",args[2],sep=""))

fasta_out=file(args[7],"w")

temp_tab=NULL

empty_sites=NULL

fl5p=NULL

if (tipo=="flankings") {
   
   blast_tab=blast_tab[!duplicated(blast_tab[,1]),]
   
   vec=as.factor(paste(">",blast_tab[,2],":",blast_tab[,9],sep=""))
   
   blast_tab=cbind(blast_tab,vec)
   
   lev=levels(blast_tab[,15])
   
} else {
   
   if (side=="3prime") lev=levels(as.factor(blast_tab[,2]))
   
   if (side=="5prime") lev=levels(as.factor(blast_tab[,2]))
   
}

if (tipo=="flankings") print(paste(Sys.time()," -> starting the 'flankings' procedure for file ",args[2],": producing assembled flankings and positions for 5prime flankings and empty sites...",sep=""))
if (tipo=="reads") print(paste(Sys.time()," -> starting the 'reads' procedure for file ",args[2],": producing assembled reads...",sep=""))

   if (operazione=="final") {
      
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
         
         if (side=="3prime"){ reads_tab_high=reads_tab[reads_tab$reads_int==temp_tab_high[1],] }
         
         if (side=="3prime"){ reads_tab_low=reads_tab[reads_tab$reads_int==temp_tab_low[1],] }
         
         reads_tab_high=t(reads_tab_high[1,])
         
         reads_tab_low=t(reads_tab_low[1,])
         
         #
         # Qui e' dove vengono prodotti gli output che ci interessano
         #
         
         if (tipo=="flankings") {
            
            writeLines(paste(">site",i,"_",side,"_",senso,sep=""),con=fasta_out)
            writeLines(reads_tab_high[2],con=fasta_out)
            if (senso=="sense") posizione=as.numeric(temp_tab_high[9])
            if (senso=="antisense") posizione=as.numeric(temp_tab_high[10])
            riga_5pf=paste(temp_tab_high[2],posizione,sep="\t")
            fl5p=rbind(fl5p,riga_5pf)
            posizione_es=posizione-100
            riga_es=paste(temp_tab_high[2],posizione_es,sep="\t")
            empty_sites=rbind(empty_sites,riga_es)
         }
         
      }
      
   }

close(fasta_out)

if (tipo=="flankings") {
   
   if (senso=="sense") {
      
      write.table(empty_sites,"empy_sites_sense_pos_list.txt",col.names=F,row.names=F,quote=F)
      
      write.table(fl5p,"5prime_pos_sense_list.txt",col.names=F,row.names=F,quote=F)
      
   }
   
   if (senso=="antisense") {
      write.table(empty_sites,"empy_sites_antisense_pos_list.txt",col.names=F,row.names=F,quote=F)
      
      write.table(fl5p,"5prime_pos_antisense_list.txt",col.names=F,row.names=F,quote=F)
   }
}

print(paste(Sys.time()," -> ",args[2]," is DONE!",sep=""))
