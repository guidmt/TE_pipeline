# Rscript generate_longer_reads_from_blasted_db_archaic_definitive_through_denovo.R args[1] args[2] args[3] args[4] args[5] args[6] args[7] args[8] args[9]
# args[1] = blastn 5prime output table (with no header and default table columns + "qlen" and "slen")
# args[2] = blastn 3prime output table (with no header and default table columns + "qlen" and "slen")
# args[3] = filtered 5prime reads database (in fasta format)
# args[4] = filtered 3prime reads database (in fasta format)
# args[5] = empty sites database (in fasta format)
# args[6] = minimum match to consider in the blastn output table
# args[7] = minimum mismatch to consider in the blastn output table
# args[8] = choose "sense" or "antisense"
# args[9] = unique prefix of the temp files to produce



args <- commandArgs(trailingOnly = TRUE)
prefisso=as.character(args[9])
print(paste(Sys.time()," -> running script generate_longer_reads_from_blasted_db_archaic_definitive_through_denovo.R in batch ",prefisso,sep=""))
library(data.table)
# inputfile=file(as.character(args[1]),"r")
# nomeinput=as.character(args[1])
# numerello=0
# while (length(blast_tab<-scan(inputfile,what="character",nlines=1,quiet=T))>0) {
#    numerello=numerello+1
#    if ((numerello%%100000)==0) {
#       print(paste(Sys.time()," -> scanned ",numerello," lines of file ",as.character(args[1]),sep=""))
#    }
# }
# close(inputfile)
# tab_5p_filt=read.table(args[1],header=F,comment.char="",colClasses=c("character","character","numeric","integer","integer","integer","integer","integer","integer","integer","numeric","numeric","integer","integer"))
tab_5p_filt=fread(args[1],header=F,sep="\t")
tab_5p_filt=as.data.frame(tab_5p_filt,comment.char="",colClasses=c("character","character","numeric","integer","integer","integer","integer","integer","integer","integer","numeric","numeric","integer","integer"))
# livelli=as.factor(tab_5p[,2])
# tab_5p_filt=NULL
# print(paste(Sys.time()," -> there are ",length(livelli)," sites to filter...",sep=""))
# for (i in 1:length(livelli)) {
#    if ((i%%100)==0) print(paste(Sys.time()," -> scanned and filtered ",i," sites...",sep=""))
#    temp_tab=tab_5p[tab_5p[,2]==lev[i],]
#    if nrow(temp_tab<=100) {
#       tab_5p_filt=rbind(tab_5p_filt,temp_tab)
#    } else {
#       temp_tab=temp_tab[order(as.numeric(temp_tab[,4]),decreasing=T),]
#       tab_5p_filt=rbind(tab_5p_filt,temp_tab[,1:50])
#       temp_tab=temp_tab[order(as.numeric(temp_tab[,13])-as.numeric(temp_tab[,4]),decreasing=T),]
#       tab_5p_filt=rbind(tab_5p_filt,temp_tab[,1:50])
#    }
# }
# print(paste(Sys.time()," -> finished scanning and filtering the 5prime tab!",sep=""))
# rm(livelli)
# rm(tab_5p)
# rm(temp_tab)
mmatch=as.numeric(args[6])
mmele=as.numeric(args[7])
senso=as.character(args[8])
tab_3p=read.table(args[2],header=F,comment.char="",colClasses=c("character","character","numeric","integer","integer","integer","integer","integer","integer","integer","numeric","numeric","integer","integer"))
# tab_3p=fread(args[2],header=F,sep=" ")
# save.image("archaic_assembly_part1.RData")
print(paste(Sys.time()," -> finished scanning the 5prime and 3prime tabs in batch ",prefisso,sep=""))
reads_5p=read.table(args[3],header=F,comment.char="",colClasses="character")
reads_3p=read.table(args[4],header=F,comment.char="",colClasses="character")
reads_es=read.table(args[5],header=F,comment.char="",colClasses="character")
lunghezza_3p=nrow(reads_3p)
reads_3p_int=reads_3p[seq(1,lunghezza_3p-1,by=2),1]
reads_3p_seq=reads_3p[seq(2,lunghezza_3p,by=2),1]
reads_3p_tab=data.frame(reads_3p_int,reads_3p_seq,stringsAsFactors=FALSE)
rm(reads_3p)
lunghezza_5p=nrow(reads_5p)
reads_5p_int=reads_5p[seq(1,lunghezza_5p-1,by=2),1]
reads_5p_seq=reads_5p[seq(2,lunghezza_5p,by=2),1]
reads_5p_tab=data.frame(reads_5p_int,reads_5p_seq,stringsAsFactors=FALSE)
rm(reads_5p)
lunghezza_es=nrow(reads_es)
reads_es_int=reads_es[seq(1,lunghezza_es-1,by=2),1]
reads_es_seq=reads_es[seq(2,lunghezza_es,by=2),1]
reads_es_tab=data.frame(reads_es_int,reads_es_seq,stringsAsFactors=FALSE)
rm(reads_es)
# save.image("archaic_assembly_part2.RData")
print(paste(Sys.time()," -> finished reading and tabularizing reads and empty sites in batch ",prefisso,sep=""))
split_es=unlist(strsplit((reads_es_tab[,1]),split="_"))
nomi_es=split_es[seq(2,length(split_es),by=2)]
reads_es_tab=cbind(reads_es_tab,nomi_es)
rm(split_es)
rm(nomi_es)
# split_3p=unlist(strsplit(reads_3p_tab[,1],split="_"))
# nomi_3p=split_3p[seq(3,length(split_3p),by=3)]
# reads_3p_tab=cbind(reads_3p_tab,nomi_3p)
# rm(split_3p)
# rm(nomi_3p)
# split_5p=unlist(strsplit(reads_5p_tab[,1],split="_"))
# nomi_5p=split_5p[seq(3,length(split_5p),by=3)]
# reads_5p_tab=cbind(reads_5p_tab,nomi_5p)
# rm(split_5p)
# rm(nomi_5p)
tab_5p_filt=tab_5p_filt[tab_5p_filt[,4]>=mmatch,]
tab_5p_filt=tab_5p_filt[(tab_5p_filt[,13]-tab_5p_filt[,4])>=mmele,]
tab_3p=tab_3p[tab_3p[,4]>=mmatch,]
tab_3p=tab_3p[(tab_3p[,13]-tab_3p[,4])>=mmele,]
split_3p_tab=unlist(strsplit(tab_3p[,2],split="_"))
nomi_3p_tab=split_3p_tab[seq(2,length(split_3p_tab),by=2)]
tab_3p=cbind(tab_3p,nomi_3p_tab)
rm(split_3p_tab)
rm(nomi_3p_tab)
split_5p_tab=unlist(strsplit(tab_5p_filt[,2],split="_"))
nomi_5p_tab=split_5p_tab[seq(2,length(split_5p_tab),by=2)]
tab_5p_filt=cbind(tab_5p_filt,nomi_5p_tab)
rm(split_5p_tab)
rm(nomi_5p_tab)
tab_3p[,1]=paste(">",tab_3p[,1],sep="")
tab_3p[,2]=paste(">",tab_3p[,2],sep="")
tab_5p_filt[,1]=paste(">",tab_5p_filt[,1],sep="")
tab_5p_filt[,2]=paste(">",tab_5p_filt[,2],sep="")
# blast_tab=rbind(tab_5p_filt,tab_3p)
# blast_tab[,1]=paste(">",blast_tab[,1],sep="")
# blast_tab[,2]=paste(">",blast_tab[,2],sep="")
# rm(tab_5p_filt)
# rm(tab_3p)
gc()
save.image(paste(prefisso,"_archaic_assembly_part1.RData",sep=""))
print(paste(Sys.time()," -> finished importing, tabularizing and preparing all the arguments, now starting the assembly in batch ",prefisso,sep=""))
# load(as.character(args[9]))
# tab_3p=fread(as.character(args[2]),sep=" ")
# tab_5p_filt=fread(as.character(args[1]),sep=" ")
# tab_3p=as.data.frame(tab_3p)
# tab_5p_filt=as.data.frame(tab_5p_filt)
reads_es_tab=reads_es_tab[((is.element(reads_es_tab[,3],tab_3p[,15]))&(is.element(reads_es_tab[,3],tab_5p_filt[,15]))),]
tab_3p=tab_3p[is.element(tab_3p[,15],reads_es_tab[,3]),]
tab_5p_filt=tab_5p_filt[is.element(tab_5p_filt[,15],reads_es_tab[,3]),]
# print(paste(Sys.time()," -> finished loading tabularized and prepared arguments for batch ",prefisso,"and now starting to assemble!",sep=""))
# lev_5p=levels(as.factor(tab_5p_filt[,15]))
# lev_3p=levels(as.factor(tab_3p[,15]))
lev=unique(reads_es_tab[,3])
# lev_3p=unique(tab_3p[,15])
# save.image(paste(prefisso,"_archaic_assembly_part4.RData",sep=""))
print(paste(Sys.time()," -> there are ",length(lev)," sites to assemble in batch ",prefisso,"...",sep=""))
assembled_3p_out=NULL
assembled_5p_out=NULL
empty_sites_out=NULL
assembled_3p_boh=NULL
assembled_5p_boh=NULL
empty_sites_boh=NULL

for (i in 1:length(lev)) {
   
   temp_tab_3p=tab_3p[tab_3p[,15]==as.character(lev[i]),]
   
   temp_3p=reads_3p_tab[is.element(reads_3p_tab[,1],temp_tab_3p[,1]),]
   
   temp_tab_5p=tab_5p_filt[tab_5p_filt[,15]==as.character(lev[i]),]
   
   temp_5p=reads_5p_tab[is.element(reads_5p_tab[,1],temp_tab_5p[,1]),]
   
   fasta_temp_3p=file("temp_3p_to_assemble.fa","w")
   
   fasta_temp_5p=file("temp_5p_to_assemble.fa","w")
   
   fasta_temp_both=file("temp_both_to_assemble.fa","w")

   
   for (k in 1:nrow(temp_3p)) {
      
      writeLines(temp_3p[k,1],con=fasta_temp_3p)
      
      writeLines(temp_3p[k,2],con=fasta_temp_3p)
      
   }
   
   close(fasta_temp_3p)
   
   for (j in 1:nrow(temp_5p)) {
      
      writeLines(temp_5p[j,1],con=fasta_temp_5p)
      
      writeLines(temp_5p[j,2],con=fasta_temp_5p)
      
   }
   
   close(fasta_temp_5p)

   
   system("cap3 temp_3p_to_assemble.fa")
   system("cap3 temp_5p_to_assemble.fa")

   system("seqkit temp_3p_to_assemble.fa.cap.contigs > temp_3p_to_assemble.fa.cap.contigs.tab.fa")
   system("seqkit temp_5p_to_assemble.fa.cap.contigs > temp_5p_to_assemble.fa.cap.contigs.tab.fa")
   
   system("wc -l temp_3p_to_assemble.fa.cap.contigs.tab.fa > temp_3p_numerello.txt")
   num3p=read.table("temp_3p_numerello.txt",header=F,)
   

   system("wc -l temp_5p_to_assemble.fa.cap.contigs.tab.fa > temp_5p_numerello.txt")
   num5p=read.table("temp_5p_numerello.txt",header=F)
   
   #
   # Criteri
   #
   
   # se uno dei due siti e' vuoto andiamo avanti
   
   # 
   #  Controlla se c'e' almeno una riga nei due file
   # 
   if ((as.numeric(as.character(num3p[1,1]))>=1)&(as.numeric(as.character(num5p[1,1]))>=1)) {
   
      #
      # Se in entrambi i file le righe sono esattamente uguale a 1 -> omozigosi
      #
      if ((nrow(assembled_3p)==1)&(nrow(assembled_5p)==1)) {
         
         assembled_3p=read.table(paste(prefisso,"temp_3p_to_assemble.fa.cap.contigs.tab.fa",sep=""),header=F,comment.char="",colClasses="character",sep="\t")
         assembled_5p=read.table(paste(prefisso,"temp_5p_to_assemble.fa.cap.contigs.tab.fa",sep=""),header=F,comment.char="",colClasses="character",sep="\t")
         
         assembled_3p_out=rbind(assembled_3p_out,paste(">assembled_fs_3prime_",as.character(lev[i]),sep=""))
         assembled_3p_out=rbind(assembled_3p_out,as.character(assembled_3p[1,2]))
         
         assembled_5p_out=rbind(assembled_5p_out,paste(">assembled_fs_5prime_",as.character(lev[i]),sep=""))
         assembled_5p_out=rbind(assembled_5p_out,as.character(assembled_5p[1,2]))
         
      } else { # no omozigosi
         
      #
      # Se in entrambi i file le righe sono esattamente uguale a 2 -> eterozigosi
      #
         
         if ((nrow(assembled_3p)==2)&(nrow(assembled_5p)==2)) {
            
            assembled_3p=read.table(paste(prefisso,"temp_3p_to_assemble.fa.cap.contigs.tab.fa",sep=""),header=F,comment.char="",colClasses="character",sep="\t")
            assembled_5p=read.table(paste(prefisso,"temp_5p_to_assemble.fa.cap.contigs.tab.fa",sep=""),header=F,comment.char="",colClasses="character",sep="\t")
            
            temp_both<-rbind(temp_3p,temp_5p)
            
            for (k in 1:nrow(temp_both)) {
               
               writeLines(temp_both[k,1],con=fasta_temp_both)
               
               writeLines(temp_both[k,2],con=fasta_temp_both)
               
            }
            
            system("cap3 fasta_temp_both.fa")
            system("seqkit fasta_temp_both.fa.cap.contigs > fasta_temp_both.cap.contigs.tab.fa")
            
            system("wc -l fasta_temp_both.cap.contigs.tab.fa > temp_both_numerello.txt")
            num_both=read.table("temp_both_numerello.txt",header=F)
            
            if(nrow(num_both)==3){
               
               assembled_both=read.table(paste(prefisso,"fasta_temp_both.cap.contigs.tab.fa",sep=""),header=F,comment.char="",colClasses="character",sep="\t")
               
               assembled_both_out=rbind(assembled_both_out,paste(">assembled_seq_a_",as.character(lev[i]),sep=""))
               assembled_both_out=rbind(assembled_both_out,as.character(assembled_3p[1,2]))
               
               assembled_both_out=rbind(assembled_both_out,paste(">assembled_seq_b_",as.character(lev[i]),sep=""))
               assembled_both_out=rbind(assembled_both_out,as.character(assembled_3p[2,2]))
               
               assembled_both_out=rbind(assembled_both_out,paste(">assembled_seq_c_",as.character(lev[i]),sep=""))
               assembled_both_out=rbind(assembled_both_out,as.character(assembled_3p[3,2]))
            
            # else 3      
               
            } else {
               
               assembled_both=read.table(paste(prefisso,"fasta_temp_both.cap.contigs.tab.fa",sep=""),header=F,comment.char="",colClasses="character",sep="\t")
               
               for(ap0 in 1:nrow(assembled_both)){
                  
               assembled_both_boh=rbind(assembled_both_boh,paste(">assembled_fs_3prime_",ap0,"_",as.character(lev[i]),sep=""))
               assembled_both_boh=rbind(assembled_both_boh,as.character(assembled_both[ap0,2]))
               
               }
               
            } # fail eterozigosi
            
      # else 2      
      } else {
            
            assembled_3p=read.table(paste(prefisso,"temp_3p_to_assemble.fa.cap.contigs.tab.fa",sep=""),header=F,comment.char="",colClasses="character",sep="\t")
            assembled_5p=read.table(paste(prefisso,"temp_5p_to_assemble.fa.cap.contigs.tab.fa",sep=""),header=F,comment.char="",colClasses="character",sep="\t")
               
            assembled_3p_int=assembled_3p[seq(1,nrow(assembled_3p)-1,by=2),1]
            assembled_3p_seq=assembled_3p[seq(2,nrow(assembled_3p),by=2),1]
            assembled_5p_int=assembled_5p[seq(1,nrow(assembled_5p)-1,by=2),1]
            assembled_5p_seq=assembled_5p[seq(2,nrow(assembled_5p),by=2),1]
            
            assembled_3p_tab=data.frame(assembled_3p_int,assembled_3p_seq,stringsAsFactors=FALSE)
            assembled_3p_tab_sort=assembled_3p_tab[order(nchar(assembled_3p_tab[,2]),decreasing=T),]
            assembled_5p_tab=data.frame(assembled_5p_int,assembled_5p_seq,stringsAsFactors=FALSE)
            assembled_5p_tab_sort=assembled_5p_tab[order(nchar(assembled_5p_tab[,2]),decreasing=T),]

            for(ap in 1:nrow(assembled_3p)){
            
               assembled_3p_boh=rbind(assembled_3p_boh,paste(">assembled_fs_3prime_",ap,"_",as.character(lev[i]),sep=""))
               assembled_3p_boh=rbind(assembled_3p_boh,as.character(assembled_3p[ap,2]))
         
         }
         
         for(ap2 in 1:nrow(assembled_5p)){
            
            assembled_5p_boh=rbind(assembled_5p_boh,paste(">assembled_fs_5prime_",ap2,"_",as.character(lev[i]),sep=""))
            assembled_5p_boh=rbind(assembled_5p_boh,as.character(assembled_5p[ap2,2]))
         
         }
         
      }
   }
}

- TUTTI GLI OUT E BOH CREATI PRIMA DEVONO AVERE UNA CONNESSIONE DI USCITA 
- DEVONO ANCHE AVERE I RISPETTIVI EMPTY SITES ESTRATTI DAI SITI VUOTI DELLA REFERENCE, DOVREBBERO MANTENERE LO STESSO SUFFISSO
   
save.image(paste(prefisso,"_archaic_assembly_part2.RData",sep=""))
print(paste(Sys.time()," -> finished assemblying the 3primes and 5primes for batch ",prefisso,", now starting selection of empty sites and assembled filled sites portions...",sep=""))

ass_3p_out_int=assembled_3p_out[seq(1,(nrow(assembled_3p_out))-1,by=2),1]
ass_3p_out_seq=assembled_3p_out[seq(2,nrow(assembled_3p_out),by=2),1]
ass_3p_out_tab=data.frame(ass_3p_out_int,ass_3p_out_seq,stringsAsFactors=F)
split_ass_3p_out=unlist(strsplit(ass_3p_out_tab[,1],split="_"))
nomi_ass_3p_out=split_ass_3p_out[seq(4,length(split_ass_3p_out),by=4)]
ass_3p_out_tab=cbind(ass_3p_out_tab,nomi_ass_3p_out)
ass_5p_out_int=assembled_5p_out[seq(1,(nrow(assembled_5p_out))-1,by=2),1]
ass_5p_out_seq=assembled_5p_out[seq(2,nrow(assembled_5p_out),by=2),1]
ass_5p_out_tab=data.frame(ass_5p_out_int,ass_5p_out_seq,stringsAsFactors=F)
split_ass_5p_out=unlist(strsplit(ass_5p_out_tab[,1],split="_"))
nomi_ass_5p_out=split_ass_5p_out[seq(4,length(split_ass_5p_out),by=4)]
ass_5p_out_tab=cbind(ass_5p_out_tab,nomi_ass_5p_out)
ass_3p_boh_int=assembled_3p_boh[seq(1,(nrow(assembled_3p_boh))-1,by=2),1]
ass_3p_boh_seq=assembled_3p_boh[seq(2,nrow(assembled_3p_boh),by=2),1]
ass_3p_boh_tab=data.frame(ass_3p_boh_int,ass_3p_boh_seq,stringsAsFactors=F)
split_ass_3p_boh=unlist(strsplit(ass_3p_boh_tab[,1],split="_"))
nomi_ass_3p_boh=split_ass_3p_boh[seq(4,length(split_ass_3p_boh),by=4)]
ass_3p_boh_tab=cbind(ass_3p_boh_tab,nomi_ass_3p_boh)
ass_5p_boh_int=assembled_5p_boh[seq(1,(nrow(assembled_5p_boh))-1,by=2),1]
ass_5p_boh_seq=assembled_5p_boh[seq(2,nrow(assembled_5p_boh),by=2),1]
ass_5p_boh_tab=data.frame(ass_5p_boh_int,ass_5p_boh_seq,stringsAsFactors=F)
split_ass_5p_boh=unlist(strsplit(ass_5p_boh_tab[,1],split="_"))
nomi_ass_5p_boh=split_ass_5p_boh[seq(4,length(split_ass_5p_boh),by=4)]

ass_5p_boh_tab=cbind(ass_5p_boh_tab,nomi_ass_5p_boh)
outesdef=file(paste("modern_empty_sites_for_archaic_specific_filled_sites_",senso,"_",prefisso,".fa",sep=""),"w")

reads_es_tab_filt=reads_es_tab[((is.element(reads_es_tab[,3],ass_3p_out_tab[,3]))&(is.element(reads_es_tab[,3],ass_5p_out_tab[,3]))),]
ass_3p_out_tab=ass_3p_out_tab[is.element(ass_3p_out_tab[,3],reads_es_tab_filt[,3]),]
ass_5p_out_tab=ass_5p_out_tab[is.element(ass_5p_out_tab[,3],reads_es_tab_filt[,3]),]
reads_es_tab_boh=reads_es_tab[((is.element(reads_es_tab[,3],ass_3p_boh_tab[,3]))&(is.element(reads_es_tab[,3],ass_5p_boh_tab[,3]))),]
ass_3p_boh_tab=ass_3p_boh_tab[is.element(ass_3p_boh_tab[,3],reads_es_tab_boh[,3]),]
ass_5p_boh_tab=ass_5p_boh_tab[is.element(ass_5p_boh_tab[,3],reads_es_tab_boh[,3]),]
system(paste("rm ",prefisso,"_temp_*",sep=""))

- fare un codice simile a reads_es_tab_filt=reads_es_tab[((is.element(reads_es_tab[,3],ass_3p_out_tab[,3]))&(is.element(reads_es_tab[,3],ass_5p_out_tab[,3]))),]
per both 

save.image(paste(prefisso,"_archaic_assembly_part3.RData",sep=""))
print(paste(Sys.time()," -> selected the empty and filled sites for batch",prefisso,", now producing fasta outputs...",sep=""))
outesdef=file(paste("modern_empty_sites_for_archaic_specific_filled_sites_",senso,"_",prefisso,".fa",sep=""),"w")
out3pdef=file(paste("archaic_specific_3prime_filled_sites_for_modern_empty_sites_",senso,"_",prefisso,".fa",sep=""),"w")
out5pdef=file(paste("archaic_specific_5prime_filled_sites_for_modern_empty_sites_",senso,"_",prefisso,".fa",sep=""),"w")

- AGGIUNGERE  OUTPUT PER GLI ETEROZIGOTI (archaic_specific_ete_both_filled_sites)
- equivalenete di outesdef per i both in eterozigosi

for (m in 1:nrow(reads_es_tab_filt)) {
   writeLines(as.character(reads_es_tab_filt[m,1]),con=outesdef)
   writeLines(as.character(reads_es_tab_filt[m,2]),con=outesdef)
   writeLines(as.character(ass_3p_out_tab[m,1]),con=out3pdef)
   writeLines(as.character(ass_3p_out_tab[m,2]),con=out3pdef)
   writeLines(as.character(ass_5p_out_tab[m,1]),con=out5pdef)
   writeLines(as.character(ass_5p_out_tab[m,2]),con=out5pdef)
   - AGGIUNGERE 
}

close(outesdef)
close(out3pdef)
close(out5pdef)

outesboh=file(paste("possible_modern_empty_sites_for_archaic_specific_filled_sites_",senso,"_",prefisso,".fa",sep=""),"w")
out3pboh=file(paste("possible_archaic_specific_3prime_filled_sites_for_modern_empty_sites_",senso,"_",prefisso,".fa",sep=""),"w")
out5pboh=file(paste("possible_archaic_specific_5prime_filled_sites_for_modern_empty_sites_",senso,"_",prefisso,".fa",sep=""),"w")

-aggiungere output both boh

for (n in 1:nrow(reads_es_tab_boh)) {
   writeLines(as.character(reads_es_tab_boh[n,1]),con=outesboh)
   writeLines(as.character(reads_es_tab_boh[n,2]),con=outesboh)
   writeLines(as.character(ass_3p_boh_tab[n,1]),con=out3pboh)
   writeLines(as.character(ass_3p_boh_tab[n,2]),con=out3pboh)
   writeLines(as.character(ass_5p_boh_tab[n,1]),con=out5pboh)
   writeLines(as.character(ass_5p_boh_tab[n,2]),con=out5pboh)
}

close(outesboh)
close(out3pboh)
close(out5pboh)
# system("rm archaic_assembly_part*")
print(paste(Sys.time()," -> batch ",prefisso, " IS DONE!     Buahahahahahaha!",sep=""))
