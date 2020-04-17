#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
#use Statistics::Descriptive;

my ($sampleFile, $outfile) = '/illumina1/YanaiLab/eitan/UTRome/wholeEmb/data/CE_Samples.tab' # @ARGV;
my $dir = '/home/eitan/UTRome/samples/';
my $barcodeFile = '/illumina1/YanaiLab/eitan/UTRome/barcodes_cel-seq.tab';

my $bc = read_barcode($barcodeFile);

## Reads Sample information and generates a) samfile and b) htseq file for each
my %sample;
open(IN, $sampleFile) || die print STDERR "cannot open $sampleFile for reading\n";
open (OUT, ">$outfile") || die print STDERR "cannot open $outfile for writing\n";
my $line = <IN>;
my @header = split ("\t", $line); 

while (<IN>)
{
	chomp;
	my @m = split ("\t", $_); my $eName = $m[0]; my $idName = $m[2]; my $sampleName = $m[3];
	for (my $i=4;$i <@header;$i++)
	{$sample{$idName}->{$header[$i]} = $m[$i];}
	 
	my $barcode = $sample{$idName}->{'barcode'};
	my $bcSeq = $bc->{$barcode};
	print STDERR "id=$idName,barcode = $bcSeq\n";	
		
	my $currDir = $dir.$idName.'/';
	print STDERR "current directory: $currDir\n";
	
	# read number of good barcodes
	my $r1_r2File = $currDir.'R1_R2.data';
	my $goodCounts = `grep -c SBS123 $r1_r2File`;
	
	#read barcode error file
	my $bcErrorFile = $currDir.'R1_barcode_error.txt';
	open (INNER, $bcErrorFile) || die print STDERR "cannot open $bcErrorFile for reading\n";
	my %bcError;
		while (<INNER>){
			chomp;
			my $header = $_;
			my ($id,@rest) = split(' ', $header);
			my $seqLine = <INNER>; 
			if ($seqLine =~ /^(\w{8}).*/){ 
				$bcError{$1}++;
			}
			my $lines = <INNER>;
			$lines .= <INNER>;
		}
		close (INNER); 
		
		my $sum = 0;
		print STDERR "printing error barcodes to file\n";
		print OUT ">$bcSeq\n";
		while( my( $key, $value ) = each %bcError ){
    		print OUT "$key\t$value\n";
    		$sum += $value;
		}
		my $perc = sprintf('%.2f', ($sum/($goodCounts+$sum)) * 100);
		print OUT "sum=$sum reads, $perc % of $idName sample \n";
}
close (IN);	
close (OUT);


sub read_barcode {
	
	my ($infile) = @_;
	my %hash;
	
	open(IN, $infile) || die print STDERR "cannot open $infile for reading\n";
	my $header = <IN>;
	
	while (<IN>){
		chomp;
		my ($id,$sequence) = split;
		$hash{$id} = $sequence;
	}
	close(IN);
	
	return \%hash;
}
