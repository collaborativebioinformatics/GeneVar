#! /usr/bin/perl
$file = "gencode.v35.annotation.gff3";
# split gencode annoation by type, so bedtools intersect will run faster
open CDS, ">cds_annotation.txt";
open EXON, ">exon_annotation.txt";
open GENE, ">gene_annotation.txt";
open START, ">start_codon_annotation.txt";
open STOP, ">stop_codon_annotation.txt";
open TRANSCRIPT, ">transcript_annotation.txt";
open UTR5, ">5utr_annotation.txt";
open UTR3, ">3utr_annotation.txt";
open STOP2, ">stop_codon_selenocysteine_annotation.txt";


open READ, $file or die print "Can't open $file\n";
while (<READ>)
{
	chomp $_;
	#skip header	 
	if($_ =~ /^#/){next;}

	@data = split/\t/,$_;
	if ($_ =~ /^chr/){$_ =~ s/^chr//;}
	if($data[2] eq "CDS")
	{
		print CDS "$_\n";
	}
	elsif($data[2] eq "exon")
	{
		print EXON "$_\n";
	}
	elsif($data[2] eq "gene")
	{
		print GENE "$_\n";
	}
	elsif($data[2] eq "stop_codon_redefined_as_selenocysteine")
	{
		print STOP2 "$_\n";
	}
	elsif($data[2] eq "start_codon")
	{
		print START "$_\n";
	}
	elsif($data[2] eq "stop_codon")
	{
		print STOP "$_\n";
	}
	elsif($data[2] eq "transcript")
	{
		print TRANSCRIPT "$_\n";
	}
	elsif($data[2] eq "five_prime_UTR")
	{
		print UTR5 "$_\n";
	}
	elsif($data[2] eq "three_prime_UTR")
	{
		print UTR3 "$_\n";
	}	
	else{print "Error - type unknown $_\n";}
}
close READ;

