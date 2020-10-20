#! /usr/bin/perl

#temp line with hard coded input file since this is for testing chr21
$file = "all_variants_chr21_clinvar.txt";
open READ, $file or die "Can't open $file\n";
open RESULTS, ">final_clinvar_dbvar_results.txt" or die "Can't open final_clinvar_dbvar_results.txt\n";
open TF, ">final_clinvar_dbvar_results_TF.txt" or die "Can't open final_clinvar_dbvar_results_TF.txt\n";
while (<READ>)
{
    chomp $_;
    @data =  split /\t/,$_;
    $data[3] =~ s/\;dbVarID\=(.*)//; $dbvar = $1; 
    $data[11] =~ s/CLNSIG=(.*)\;CLNVC//; $clinvar = $1; $clinvar =~ s/\;CLNVC.*$//; $clinvar =~ s/\;CLNSIGCONF.*$//; 
    if ($clinvar =~ /^ns/)
    {
	# these records have no CLNSIG, instead CLNSIGINCL
	$data[11] =~ s/.*CLNSIGINCL=(.*)//; $cli = $1; $cli =~ s/^.*://; $clinvar = $cli;
    }
    # print dbvar id and the clinvar significance  
    print RESULTS "$dbvar\t$clinvar\n";
    if ( ($clinvar eq "Pathogenic") || ($clinvar eq "Likely_pathogenic") || ($clinvar eq "Pathogenic/Likely_pathogenic") )
    {
	$truth = "TRUE";
    }
    else{$truth = "FALSE";}
    # print dbvar id and simple True or false if it is P/LP
    print TF "$dbvar\t$truth\n";
}
close READ;
