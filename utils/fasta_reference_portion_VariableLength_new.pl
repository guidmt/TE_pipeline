#!/usr/bin/perl -w

use strict;

# ARGV[0]=poslist (tab separated chr\tpos\tlength\n)
# ARGV[1]=output name.fasta
# ARGV[2]= up or down

my $i;
my $j;
my $k;

my @chrlines;
my $string="";
my $length_var="";

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
    $length_var=$splitsnp[2];
    #print STDERR "For $splitsnp[0] and $splitsnp[1] my Distance is $ARGV[2]!\n";
    if ($ARGV[2] eq "up"){
    	 $pos= ($splitsnp[1] - 1 - $length_var);
	}

 if ($ARGV[2] eq "down"){
         $pos= ($splitsnp[1] - 1 );
        }
#print STDERR "your pos is $splitsnp[1] and my pos is $pos\n";
my $end= $pos + $length_var;
    my $base= uc(substr $string, $pos, $length_var);
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
    
    

