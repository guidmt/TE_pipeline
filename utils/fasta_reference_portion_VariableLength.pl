#!/usr/bin/perl -w

use strict;

# ARGV[0]=poslist (tab separated chr\tpos\tlength\n)
# ARGV[1]=output name.fasta
# ARGV[2]=casual number, not important
# ARGV[3]= up or down
my $i;
my$j;
my$k;
#my $chr=$ARGV[1];
my @chrlines;
#my $chrlines;
my $string="";
#my $out=$chr."recAfterAncestral.txt";
#my $anc=$chr."ancestralAllele.txt";
#open (ANC, ">$anc");
open (OUT, ">$ARGV[1]");
my $chr="-999";
my $pos;
open (SNPS, "<$ARGV[0]");
my @snplines= <SNPS>;
chomp @snplines;

for $i (0 .. $#snplines){
    #print STDERR "hey! $i\n";
    my @splitsnp= split (/\s+/, $snplines[$i]);
    if ($splitsnp[0] ne $chr){
	my $file="MyNewChr".$splitsnp[0].".fasta";
	$chr=$splitsnp[0];
	$string="";
	open (CHR, "<$file");
	my @chrlines= <CHR>;
	chomp@chrlines;
		for $i (1 .. $#chrlines){
    		$string= $string.$chrlines[$i];
		}
	close (CHR);
	#next;
    }
    #if ($splitsnp[2] != $chr){
    #close (CHR);
    #@chrlines= ();
	#$chr= $splitsnp[2];
    $ARGV[2]=$splitsnp[2];
    #print STDERR "For $splitsnp[0] and $splitsnp[1] my Distance is $ARGV[2]!\n";
    if ($ARGV[3] eq "up"){
    	 $pos= ($splitsnp[1] - 1 - $ARGV[2]);
	}

 if ($ARGV[3] eq "down"){
         $pos= ($splitsnp[1] - 1 );
        }
#print STDERR "your pos is $splitsnp[1] and my pos is $pos\n";
my $end= $pos + $ARGV[2];
    my $base= uc(substr $string, $pos, $ARGV[2]);
    #if ($base eq "N" or $base eq "-" or $base eq "."){
# next;
    #}
my $h=">".$splitsnp[0].":".($pos+1)."-".$end;
print OUT "$h\n$base\n";
    #print ANC "$splitsnp[0]\t$base\n";
    #print OUT "$splitsnp[1]\t$base\n";
}
close (SNPS);
close (OUT);
#close (ANC);
print STDERR "FINISHED!\n";
    
    

