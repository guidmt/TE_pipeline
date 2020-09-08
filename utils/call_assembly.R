# Rscript generate_longer_reads_from_blasted_db_modern_definitive.R args[1] args[2] args[3] args[4] args[5] args[6]
# args[1] = blastn output table of 3prime modern flankings against all the archaic reads (with no header and default table columns + qlen + slen)
# args[2] = blastn output table of 5prime modern flankings against all the archaic reads (with no header and default table columns + qlen + slen)
# args[3] = complete reads database (in fasta format)
# args[4] = modern specific filled sites file (in fasta format)
# args[5] = choose "sense" or "antisense"
# args[6] = minimum length of the match

args <- commandArgs(trailingOnly = TRUE)
print(paste(Sys.time()," -> running script generate_longer_reads_from_blasted_db_modern_definitive.R on file ",args[3]," and its associates...",sep=""))
blast_tab_3p=fread(args[1],header=F,data.table=F)
blast_tab_5p=fread(args[2],header=F,data.table=F)

reads_tab=fread(args[3],header=F,data.table=F)
fs_tab=fread(args[4],header=F,data.table=F)

sense=as.character(args[5])
mmatch=as.numeric(args[6])

print(paste(Sys.time()," -> finished importing all the arguments for file ",args[3]," and its associates...",sep=""))

blast_tab_3p_filt=blast_tab_3p[as.numeric(blast_tab_3p[,4])>=mmatch,]
blast_tab_5p_filt=blast_tab_5p[as.numeric(blast_tab_5p[,4])>=mmatch,]

lungh=nrow(filled_sites)

nominativi_fs=NULL
nominativi_3p=NULL
nominativi_5p=NULL

temp=unlist(strsplit(fs_tab[,1],split="_"))
nominativi_fs=temp[seq(2,length(temp),2)]

fs_tab=cbind(fs_tab,nominativi_fs)
lev=levels(as.factor(nominativi_fs))

print(paste(Sys.time()," -> finished filtering and/or tabularizing the arguments for file ",args[3]," and its associates...",sep=""))
print(paste(Sys.time()," -> starting operation FINAL on file ",args[3]," and its associates...",sep=""))

fasta_out_def=NULL
fasta_fs_filt_def=NULL
fasta_out_part=NULL
fasta_fs_filt_part=NULL

for (i in 1:length(lev)) {
   
   temp_tab_3p=blast_tab_3p_filt[blast_tab_3p_filt[,15]==as.character(lev[i]),]
   
   temp_tab_5p=blast_tab_5p_filt[blast_tab_5p_filt[,15]==as.character(lev[i]),]  
   
   rtemp=reads_tab[is.element(reads_tab[,1],temp_tab_3p[,2])|is.element(reads_tab[,1],temp_tab_5p[,2]),]
   
   rtemp=rtemp[!duplicated(rtemp[,1]),]
   
   fasta_temp=file("temp_to_assemble.fa","w")
   for (k in 1:nrow(rtemp)) {
      writeLines(rtemp[k,1],con=fasta_temp)
      writeLines(rtemp[k,2],con=fasta_temp)
   }
   close(fasta_temp)
   system("ABYSS -k40 temp_to_assemble.fa -o temp_assembled.fa")
   system("wc -l temp_assembled.fa > temp_numerello.txt")
   
   numzi=read.table("temp_numerello.txt",header=F)
   
   if (numzi[1,1]>=2) {
      assembled=read.table("temp_assembled.fa",header=F,comment.char="",colClasses="character",sep="\t")
      if (nrow(assembled)==2) {
         fasta_fs_filt_def=rbind(fasta_fs_filt_def,unique(fs_tab[fs_tab[,3]==lev[i],1]))
         fasta_fs_filt_def=rbind(fasta_fs_filt_def,as.character(unique(fs_tab[fs_tab[,3]==lev[i],2])))
         fasta_out_def=rbind(fasta_out_def,paste(">assembled_es_",lev[i],sep=""))
         fasta_out_def=rbind(fasta_out_def,as.character(assembled[2,1]))
      }
      if (nrow(assembled)>2) {
         assembled_int=assembled[seq(1,nrow(assembled)-1,by=2),1]
         assembled_seq=assembled[seq(2,nrow(assembled),by=2),1]
         assembled_tab=data.frame(assembled_int,assembled_seq,stringsAsFactors=FALSE)
         assembled_tab_sort=assembled_tab[order(nchar(assembled_tab[,2]),decreasing=T),]
         fasta_fs_filt_part=rbind(fasta_fs_filt_part,unique(fs_tab[fs_tab[,3]==lev[i],1]))
         fasta_fs_filt_part=rbind(fasta_fs_filt_part,as.character(unique(fs_tab[fs_tab[,3]==lev[i],2])))
         fasta_out_part=rbind(fasta_out_part,paste(">assembled_es_1_",lev[i],sep=""))
         fasta_out_part=rbind(fasta_out_part,as.character(assembled_tab_sort[1,2]))
         fasta_out_part=rbind(fasta_out_part,paste(">assembled_es_2_",lev[i],sep=""))
         fasta_out_part=rbind(fasta_out_part,as.character(assembled_tab_sort[2,2]))
      }
   }
}

system("rm temp*")
write.table(fasta_fs_filt_def,paste("modern_specific_filled_sites_for_archaic_empty_sites_",sense,".fa",sep=""),col.names=F,row.names=F,quote=F)
write.table(fasta_out_def,paste("archaic_empty_sites_for_modern_specific_filled_sites_",sense,".fa",sep=""),col.names=F,row.names=F,quote=F)
write.table(fasta_fs_filt_part,paste("possible_modern_specific_filled_sites_for_archaic_empty_sites_",sense,".fa",sep=""),col.names=F,row.names=F,quote=F)
write.table(fasta_out_part,paste("possible_archaic_empty_sites_for_modern_specific_filled_sites_",sense,".fa",sep=""),col.names=F,row.names=F,quote=F)
print(paste(Sys.time()," -> finished all operations on file ",args[3]," and its associates!!! YEAHYEAHYEAH!!!",sep=""))
   
