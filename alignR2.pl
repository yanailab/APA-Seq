#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Statistics::Descriptive;
use IO::Uncompress::Gunzip qw($GunzipError);

my ($sampleFile) = @ARGV; 
my $dir = '/illumina1/YanaiLab/eitan/UTRome/wholeEmb/samples/';
my $genedir = '/illumina1/YanaiLab/eitan/UTRome/genes/';

my %sample = read_samples($sampleFile);
my @samples = sort keys %sample;

for (my $s = 0; $s < @samples; $s++){
	
	my $idName = $samples[$s];
	my $currDir = $dir.$idName.'/';
	print STDERR "current directory: $currDir\n";

	my $htseqFile = $currDir.$idName.'.htseq_count.txt';
	my $samFile = $currDir.$idName.'.sam';
	my $R2_plus_dir = $currDir.'R2_plus';
	unless (-d $R2_plus_dir) {system "mkdir $R2_plus_dir";}
	
	open (IN, $htseqFile) || die print STDERR "cannot open $htseqFile for reading\n"; # append in case there are multiple lanes 
	while (<IN>){
		chomp;
		my ($gene,$count) = split;
		my $r2_geneFile = $R2_plus_dir.'/'.$gene;
		
		if ($count > 0) {	
			system "grep $gene $samFile \| awk \'{print \$10}\' \>  $r2_geneFile";
			my $geneFastaFile = $genedir.$gene.'.fasta';
			my $btFilesPathPrefix = $genedir.$gene;
			if (-e $geneFastaFile){
				my $r2_samFile = $R2_plus_dir.'/'.$gene.'.sam';
				my $r2_samstatFile = $r2_samFile.'.stat';
				system "bowtie2 -x $btFilesPathPrefix -r $r2_geneFile > $r2_samFile 2>$r2_samstatFile";
			}
		}
	} 
	close (IN);
}
	
sub read_samples {
	my ($infile) = @_;
	my %hash;
	open(IN, $infile) || die print STDERR "cannot open $infile for reading\n";
	my $line = <IN>;
	my @header = split ("\t", $line); 
	
	while (<IN>){
		chomp;	
		my @m = split ("\t", $_);  my $sampleName = $m[3]; my $eName = $m[0]; my $idName = $m[2];
		for (my $i=4;$i <@header;$i++)
		{$hash{$idName}->{$header[$i]} = $m[$i];}
	}
	close (IN);
	
	return %hash;
}
