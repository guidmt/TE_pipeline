# Rscript deduplicate_modern_flankings_and_filled_sites.R args[1] args[2] args[3] args[4]
# args[1] = modern 3prime flanking fasta file prefix (before "_sense/antisense.fa")
# args[2] = modern 5prime flanking fasta file prefix (before "_sense/antisense.fa")
# args[3] = modern filled sites fasta file prefix (before "_sense/antisense.fa")
# args[4] = choose "sense" or "antisense"



args <- commandArgs(trailingOnly = TRUE)
print(paste(Sys.time()," -> running script deduplicate_modern_flankings_and_filled_sites.R",sep=""))
senso=as.character(args[4])

filled_sites=read.table(paste(args[3],"_",senso,".fa",sep=""),header=F,comment.char="",colClasses="character")
sites_5p=read.table(paste(args[2],"_",senso,".fa",sep=""),header=F,comment.char="",colClasses="character")
sites_3p=read.table(paste(args[1],"_",senso,".fa",sep=""),header=F,comment.char="",colClasses="character")

print(paste(Sys.time()," -> finished importing all the arguments",sep=""))
lun_fs=nrow(filled_sites)
fs_int=filled_sites[seq(1,lun_fs-1,by=2),1]
fs_seq=filled_sites[seq(2,lun_fs,by=2),1]
fs_tab=data.frame(fs_int,fs_seq,stringsAsFactors=FALSE)

# fs_tab$fs_int=sub(">","",fs_tab$fs_int)
lun_5p=nrow(sites_5p)
int_5p=sites_5p[seq(1,lun_5p-1,by=2),1]
seq_5p=sites_5p[seq(2,lun_5p,by=2),1]
tab_5p=data.frame(int_5p,seq_5p,stringsAsFactors=FALSE)

# tab_5p$int_5p=sub(">","",tab_5p$int_5p)
lun_3p=nrow(sites_3p)
int_3p=sites_3p[seq(1,lun_3p-1,by=2),1]
seq_3p=sites_3p[seq(2,lun_3p,by=2),1]
tab_3p=data.frame(int_3p,seq_3p,stringsAsFactors=FALSE)

print(paste(Sys.time()," -> finished tabularizing the fasta files",sep=""))
vettore=!duplicated(tab_5p[,1])
tab_3p_sep_filt=tab_3p[vettore,]
tab_5p_sep_filt=tab_5p[vettore,]
fs_tab_sep_filt=fs_tab[vettore,]

print(paste(Sys.time()," -> finished de-duplicating the fasta files, now producing output files...",sep=""))
fasta_fs_filt=file(paste(args[3],"_",senso,".fa",sep=""),"w")
fasta_5p_filt=file(paste(args[2],"_",senso,".fa",sep=""),"w")
fasta_3p_filt=file(paste(args[1],"_",senso,".fa",sep=""),"w")

apply(fs_tab_sep_filt,1,FUN=function(X){
     writeLines(X[1],con=fasta_fs_filt)
     writeLines(X[2],con=fasta_fs_filt)
}
)

apply(tab_5p_sep_filt,1,FUN=function(X){
        writeLines(X[1],con=fasta_5p_filt)
        writeLines(X[2],con=fasta_5p_filt)
}
)

apply(tab_3p_sep_filt,1,FUN=function(X){
        writeLines(X[1],con=fasta_3p_filt)
        writeLines(X[2],con=fasta_3p_filt)
}
)

close(fasta_fs_filt)
close(fasta_5p_filt)
close(fasta_3p_filt)
print(paste(Sys.time()," -> DONE!",sep=""))
