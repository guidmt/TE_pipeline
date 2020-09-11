library(data.table)

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
# args[10] = cap3 (clipping 0 No 1 yes)

args <- commandArgs(trailingOnly = TRUE)
prefix=as.character(args[9])
print(paste(Sys.time()," -> running script generate_longer_reads_from_blasted_db_archaic_definitive_through_denovo.R in batch ",prefix,sep=""))

tab_5p_filt=fread(args[1],header=F,sep="\t",data.table=T)

# get parameters minimum match and mismatch
mmatch=as.numeric(args[6])
mmele=as.numeric(args[7])
sense=as.character(args[8])
reads_5p_tab=fread(args[2],header=F,sep="\t",data.table=T)

# get the file filtered 5prime, filtered 3prime
reads_5p=fread(args[3],header=F,data.table=T)
reads_3p_tab=fread(args[4],header=F,data.table=T)

# get the file with empty site
reads_es_tab=fread(args[5],header=F,data.table=T)

tab_5p_filt=tab_5p_filt[tab_5p_filt[,4]>=mmatch,]
tab_5p_filt=tab_5p_filt[(tab_5p_filt[,13]-tab_5p_filt[,4])>=mmele,]

tab_3p=tab_3p[tab_3p[,4]>=mmatch,]
tab_3p=tab_3p[(tab_3p[,13]-tab_3p[,4])>=mmele,]

split_3p_tab=unlist(strsplit(tab_3p[,2],split="_"))
name_3p_tab=split_3p_tab[seq(2,length(split_3p_tab),by=2)]

tab_3p=cbind(tab_3p,name_3p_tab)
rm(split_3p_tab)
rm(name_3p_tab)

split_5p_tab=unlist(strsplit(tab_5p_filt[,2],split="_"))
name_5p_tab=split_5p_tab[seq(2,length(split_5p_tab),by=2)]

tab_5p_filt=cbind(tab_5p_filt,name_5p_tab)

rm(split_5p_tab)
rm(name_5p_tab)

tab_3p[,1]=paste(">",tab_3p[,1],sep="")
tab_3p[,2]=paste(">",tab_3p[,2],sep="")

tab_5p_filt[,1]=paste(">",tab_5p_filt[,1],sep="")
tab_5p_filt[,2]=paste(">",tab_5p_filt[,2],sep="")

reads_es_tab=reads_es_tab[((is.element(reads_es_tab[,3],tab_3p[,15]))&(is.element(reads_es_tab[,3],tab_5p_filt[,15]))),]
tab_3p=tab_3p[is.element(tab_3p[,15],reads_es_tab[,3]),]
tab_5p_filt=tab_5p_filt[is.element(tab_5p_filt[,15],reads_es_tab[,3]),]

lev=unique(reads_es_tab[,3])

print(paste(Sys.time()," -> there are ",length(lev)," sites to assemble in batch ",prefix,"...",sep=""))

assembled_3p_out=NULL
assembled_5p_out=NULL
empty_sites_out=NULL

assembled_3p_out2=NULL
assembled_5p_out2=NULL
empty_sites_out2=NULL

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

   system("seqkit fx2tab temp_3p_to_assemble.fa.cap.contigs > temp_3p_to_assemble.fa.cap.contigs.tab.fa")
   system("seqkit fx2tab temp_5p_to_assemble.fa.cap.contigs > temp_5p_to_assemble.fa.cap.contigs.tab.fa")
   
   system("wc -l temp_3p_to_assemble.fa.cap.contigs.tab.fa > temp_3p_number_var.txt")
   num3p=fread("temp_3p_number_var.txt",header=F,data.table=T)
   

   system("wc -l temp_5p_to_assemble.fa.cap.contigs.tab.fa > temp_5p_number_var.txt")
   num5p=fread("temp_5p_number_var.txt",header=F,data.table=T)
   
   #
   # Criteri
   #
   

   # 
   #  Controlla se c'e' almeno una riga nei due file
   # 
   if ((as.numeric(nrow(num3p))>=1)&(as.numeric(nrow(num5p))>=1)) {
   
      #
      # Se in entrambi i file le righe sono esattamente uguale a 1 -> omozigosi
      #
      if ((nrow(assembled_3p)==1)&(nrow(assembled_5p)==1)) {
         
         ass_3p_out_tab=fread(paste(prefix,"temp_3p_to_assemble.fa.cap.contigs.tab.fa",sep=""),header=F,data.table=F)
         ass_5p_out_tab=fread(paste(prefix,"temp_5p_to_assemble.fa.cap.contigs.tab.fa",sep=""),header=F,data.table=F)
         
         ass_3p_out_tab[1,1]<-paste(">assembled_fs_3prime_",as.character(lev[i]),sep="")
         ass_5p_out_tab[1,1]<-paste(">assembled_fs_5prime_",as.character(lev[i]),sep="")
         
      } else { # no omozigosi
         
      #
      # Se in entrambi i file le righe sono esattamente uguale a 2 -> eterozigosi
      #
         
         if ((nrow(assembled_3p)==2)&(nrow(assembled_5p)==2)) {
            
            assembled_3p=fread(paste(prefix,"temp_3p_to_assemble.fa.cap.contigs.tab.fa",sep=""),header=F,data.table=T)
            assembled_5p=fread(paste(prefix,"temp_5p_to_assemble.fa.cap.contigs.tab.fa",sep=""),header=F,data.table=T)
            
            temp_both<-rbind(temp_3p,temp_5p)
      
            #    Qui vedere output CAP
            
            for (k in 1:nrow(temp_both)) {
               
               writeLines(temp_both[k,1],con=fasta_temp_both)
               
               writeLines(temp_both[k,2],con=fasta_temp_both)
               
            }
            
            #  dove fasta_temp_both.fa e' genereato
            system(paste("cap3 -k", args[10],"fasta_temp_both.fa"))


            system("seqkit fx2tab fasta_temp_both.fa.cap.contigs > fasta_temp_both.cap.contigs.tab.fa")
            
            system("wc -l fasta_temp_both.cap.contigs.tab.fa > temp_both_number_var.txt")
            
            num_both=fread("temp_both_number_var.txt",data.table=T,header=T)
            
            if(nrow(num_both)==3){
               
               assembled_both=fread(paste(prefix,"fasta_temp_both.cap.contigs.tab.fa",sep=""),data.table=T,header=F)
               
               assembled_both_out=rbind(assembled_both_out,paste(">assembled_seq_a_",as.character(lev[i]),sep=""))
               assembled_both_out=rbind(assembled_both_out,as.character(assembled_3p[1,2]))
               
               assembled_both_out=rbind(assembled_both_out,paste(">assembled_seq_b_",as.character(lev[i]),sep=""))
               assembled_both_out=rbind(assembled_both_out,as.character(assembled_3p[2,2]))
               
               assembled_both_out=rbind(assembled_both_out,paste(">assembled_seq_c_",as.character(lev[i]),sep=""))
               assembled_both_out=rbind(assembled_both_out,as.character(assembled_3p[3,2]))
            
            # else 3      
               
            } else {
               
               assembled_both=fread(paste(prefix,"fasta_temp_both.cap.contigs.tab.fa",sep=""),header=F,comment.char="",colClasses="character",sep="\t")
               
               for(ap0 in 1:nrow(assembled_both)){
                  
               assembled_both_out2=rbind(assembled_both_out2,paste(">assembled_fs_3prime_",ap0,"_",as.character(lev[i]),sep=""))
               assembled_both_out2=rbind(assembled_both_out2,as.character(assembled_both[ap0,2]))
               
               }
               
            } # fail eterozigosi
            
      # else 2      
      } else {
            
         assembled_3p_tab=fread(paste(prefix,"temp_3p_to_assemble.fa.cap.contigs.tab.fa",sep=""),header=F,data.table=T)
         assembled_5p_tab=fread(paste(prefix,"temp_5p_to_assemble.fa.cap.contigs.tab.fa",sep=""),header=F,data.table=T)
               
         assembled_3p_tab_sort=assembled_3p_tab[order(nchar(assembled_3p_tab[,2]),decreasing=T),]
         assembled_5p_tab_sort=assembled_5p_tab[order(nchar(assembled_5p_tab[,2]),decreasing=T),]

            for(ap in 1:nrow(assembled_3p)){
            
               assembled_3p_out2=rbind(assembled_3p_out2,paste(">assembled_fs_3prime_",ap,"_",as.character(lev[i]),sep=""))
               assembled_3p_out2=rbind(assembled_3p_out2,as.character(assembled_3p[ap,2]))
         
         }
         
         for(ap2 in 1:nrow(assembled_5p)){
            
            assembled_5p_out2=rbind(assembled_5p_out2,paste(">assembled_fs_5prime_",ap2,"_",as.character(lev[i]),sep=""))
            assembled_5p_out2=rbind(assembled_5p_out2,as.character(assembled_5p[ap2,2]))
         
         }
         
      }
   }
   }
}

print(paste(Sys.time()," -> finished assemblying the 3primes and 5primes for batch ",prefix,", now starting selection of empty sites and assembled filled sites portions...",sep=""))

# Define a function to save the results
get_tab_from_seq<-function(x){
   
   x_int=x[seq(1,(nrow(x))-1,by=2),1]
   x_seq=x[seq(2,nrow(x),by=2),1]
   x_tab=data.frame(x_int,x_seq,stringsAsFactors=F)
   x_out2=unlist(strsplit(x_tab[,1],split="_"))
   x_out2=x_out2[seq(4,length(x_out2),by=4)]
   x_tab=cbind(x_tab,x_out2)
   
   return(x_tab)
      
}

# Define a function to perform subselection of two data.frames
save_res_3p_5p<-function(es_out,ass_3p_out,ass_5_out,outesdef,out3pdef,out5pdef){
   
   es_out_filt=x[((is.element(x[,3],ass_3p_out[,3]))&(is.element(x[,3],ass_5_out[,3]))),]
   ass_3p_out=ass_3p_out[is.element(ass_3p_out[,3],es_out[,3]),]
   ass_5_out=ass_5_out[is.element(ass_5_out[,3],es_out_filt[,3]),]
   
   for(m in 1:nrow(ass_3p_out)){
      
      writeLines(as.character(es_out_filt[m,1]),con=outesdef)
      writeLines(as.character(es_out_filt[m,2]),con=outesdef)
      writeLines(as.character(ass_3p_out[m,1]),con=out3pdef)
      writeLines(as.character(ass_3p_out[m,2]),con=out3pdef)
      writeLines(as.character(ass_5_out[m,1]),con=out5pdef)
      writeLines(as.character(ass_5_out[m,2]),con=out5pdef)
      
   }
   
   close(outesdef)
   close(out3pdef)
   close(out5pdef)
   
}


# output 3p_out_tab
split_ass_3p_out=unlist(strsplit(ass_3p_out_tab[,1],split="_"))
name_ass_3p_out=split_ass_3p_out[seq(4,length(split_ass_3p_out),by=4)]
ass_3p_out_tab=cbind(ass_3p_out_tab,name_ass_3p_out)

# output 5p_out_tab
split_ass_5p_out=unlist(strsplit(ass_5p_out_tab[,1],split="_"))
name_ass_5p_out=split_ass_5p_out[seq(4,length(split_ass_5p_out),by=4)]
ass_5p_out_tab=cbind(ass_5p_out_tab,name_ass_5p_out)

# output 3p_out2
ass_3p_out2_tab<-get_tab_from_seq(assembled_3p_out2)

# output 5p_out2
ass_5p_out2_tab<-get_tab_from_seq(assembled_5p_out2)


outesdef=file(paste("modern_empty_sites_for_archaic_specific_filled_sites_",sense,"_",prefix,".fa",sep=""),"w")

#
# Processing empty sites out 1 
#

outesdef=file(paste("modern_empty_sites_for_archaic_specific_filled_sites_",sense,"_",prefix,".fa",sep=""),"w")
out3pdef=file(paste("archaic_specific_3prime_filled_sites_for_modern_empty_sites_",sense,"_",prefix,".fa",sep=""),"w")
out5pdef=file(paste("archaic_specific_5prime_filled_sites_for_modern_empty_sites_",sense,"_",prefix,".fa",sep=""),"w")


save_res_3p_5p(es_out=reads_es_tab,
               ass_3p_out=ass_3p_out_tab,
               ass_3p_out=ass_5p_out_tab,
               outesdef=outesdef,
               out3pdef=out3pdef,
               out5pdef=out5pdef)

#
# Processing empty sites out 2
#

outesout2=file(paste("possible_modern_empty_sites_for_archaic_specific_filled_sites_",sense,"_",prefix,".fa",sep=""),"w")
out3pout2=file(paste("possible_archaic_specific_3prime_filled_sites_for_modern_empty_sites_",sense,"_",prefix,".fa",sep=""),"w")
out5pout2=file(paste("possible_archaic_specific_5prime_filled_sites_for_modern_empty_sites_",sense,"_",prefix,".fa",sep=""),"w")
   
save_res_3p_5p(es=reads_es_tab,
              out2_3p=ass_3p_out2_tab,
              out2_5p=ass_5p_out2_tab,
              outesout2=outesout2,
              out3pout2=out3pout2,
              out5pout2=out5pout2)

#
# Processing empty sites out 3
#

outesout2=file(paste("both_modern_empty_sites_for_archaic_specific_filled_sites_",sense,"_",prefix,".fa",sep=""),"w")
out3pout2=file(paste("both_archaic_specific_3prime_filled_sites_for_modern_empty_sites_",sense,"_",prefix,".fa",sep=""),"w")
out5pout2=file(paste("both_archaic_specific_5prime_filled_sites_for_modern_empty_sites_",sense,"_",prefix,".fa",sep=""),"w")

save_res_3p_5p(es=reads_es_tab,
               out2_3p=assembled_3p_out2,
               out2_5p=assembled_3p_out2,
               outesout2=outesout2,
               out3pout2=out3pout2,
               out5pout2=out5pout2)

print(paste(Sys.time()," -> batch ",prefix, " IS DONE!     Buahahahahahaha!",sep=""))

#- TUTTI GLI OUT E out2 CREATI PRIMA DEVONO AVERE UNA CONNESSIONE DI USCITA 
#- DEVONO ANCHE AVERE I RISPETTIVI EMPTY SITES ESTRATTI DAI SITI VUOTI DELLA REFERENCE, DOVREBBERO MANTENERE LO STESSO SUFFISSO

