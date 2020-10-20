#! /usr/bin/perl

# temp line, hard code input file since test is on chr21 only
$file = "all_variants_chr21.tsv";
open READ, $file or die "Can't open $file\n";
open BED, ">all_variants_chr21.bed" or die "Can't open all_variants_chr21.bed\n";
while (<READ>)
{
    chomp $_;
    if($_ =~ /^chr/){next;}
    @data =  split /\t/,$_;
    #removing BND that was accidently add to the all_variants_chr21.tsv
    if ($data[3] eq "BND"){next;}
    print BED "$data[0]\t$data[1]\t$data[2]\ttype=$data[3];dbVarID=$data[4]\n";
}
close READ;
