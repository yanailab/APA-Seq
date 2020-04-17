#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
#use Statistics::Descriptive;

my ($sampleFile) = @ARGV;
my $dir = '/illumina1/YanaiLab/eitan/UTRome/wholeEmb/samples/';
my $annotFile = '/illumina1/YanaiLab/new_fs/refs/annotations/CE/WS230/c_elegans.WS230_spikein.annotations_trimmed.spikes_and_lincs.gff3';
my $fastqDir = '/illumina1/YanaiLab/new_fs/projects/metazome/raw_files/CE/CEL-Seq/';
my $barcodeFile = '/illumina1/yanailab/eitan/UTRome/barcodes_cel-seq.tab';

my $bc = read_barcode($barcodeFile);

## Reads Sample information and generates a) samfile and b) htseq file for each
my %sample;
open(IN, $sampleFile) || die print STDERR "cannot open $sampleFile for reading\n";
my $line = <IN>;
my @header = split ("\t", $line); 

while (<IN>)
{
	chomp;
	my @m = split ("\t", $_); my $eName = $m[0]; my $idName = $m[2]; my $sampleName = $m[3];
	for (my $i=4;$i <@header;$i++)
	{$sample{$idName}->{$header[$i]} = $m[$i];}
	 
	my $currDir = $dir.$idName.'/';
	print STDERR "current directory: $currDir\n";
	unless (-d $currDir) {system "mkdir $currDir";}
	 
	my $initSam = $currDir.$idName.'.init.sam';
	my $osam = $currDir.$idName.'.sam';	
	
	unless ((-e $initSam) || (-e  $osam)){
		print STDERR "generating initSam $initSam\n"; 
		system "samtools view -h $sample{$idName}->{'BWA_mapping_file'} > $initSam ";
	}

	unless (-e $osam){
		print STDERR "$osam doesn't exist\n";
		print STDERR "generating sam $osam and htseq-count file for $sampleName\n";
		my $htseqFile = $currDir.$idName.'.htseq_count.txt';
		system "htseq-count -o $osam $initSam $annotFile > $htseqFile";
		system "rm $initSam";
	}
	
	my $r2file = $currDir.'R2.data';
	unless (-e $r2file){
		system "cat $osam \| awk \'{print \$1, \$NF}\' \> $r2file";
	}	
	
	
	#read fatsq file instead of R1 file (check that barcode is as it should be, control for multiplexing)
	my $combFile = $currDir.'R1_R2.data';
	my $bcErrorFile = $currDir.'R1_barcode_error.txt';

	unless (-e $combFile){
	open (BCERR, ">$bcErrorFile") || die print STDERR "cannot open $bcErrorFile for writing\n";
	open (OUTER, ">$combFile") || die print STDERR "cannot open $combFile for writing\n";		
		print STDERR "generating $combFile\n";
		
		my @lanes = ('L007', 'L008');
		my $index = $sample{$idName}->{'index_id'}; 
		if ($index ==9){$index = '09_GATCAG_';} elsif ($index == 10){$index = '10_TAGCTT_';} 
			else {die print STDERR "no index received for $idName\n";}
		my $barcode = $sample{$idName}->{'barcode'};
		my $bcSeq = $bc->{$barcode};
		print STDERR "id=$idName,barcode = $bcSeq\n";
		
		my %read2xf;
		open (INNER, $r2file) || die print STDERR "cannot open $r2file for reading\n";
		while (<INNER>){
			my ($read,$xf) = split;
			$read2xf{$read} = $xf;
		}
		close (INNER);
	 
		my $fastqFile1 = $fastqDir.'index'.$index.$lanes[0].'_R1_001.fastq';
		open (INNER, $fastqFile1) || die print STDERR "cannot open $fastqFile1 for reading\n";
		print STDERR "reading $fastqFile1 for $idName\n";
		
		while (<INNER>){
			chomp;
			my $header = $_;
			my ($id,@rest) = split(' ', $header);
			$id =~ s/\@//;
			my $lines = '';
			my $seqLine = <INNER>; $lines .= $seqLine;
			$lines .= <INNER>;
			$lines .= <INNER>;
		
        	if (defined $read2xf{$id}) {
	        	if ($seqLine =~ /^$bcSeq/) { 
	        	print OUTER "$header\n$lines $read2xf{$id}\n\n";
        		}else{
	       		 	print BCERR "$header\n$lines";
	       		 } 	
			}
		}
		close (INNER); 
		
		my $fastqFile2 = $fastqDir.'index'.$index.$lanes[1].'_R1_001.fastq';
		open (INNER, $fastqFile2) || die print STDERR "cannot open $fastqFile2 for reading\n";
		print STDERR "reading $fastqFile2 for $idName\n";
		
		while (<INNER>){
			chomp;
			my $header = $_;
			my ($id,@rest) = split(' ', $header);
			$id =~ s/\@//;
			my $lines = '';
			my $seqLine = <INNER>; $lines .= $seqLine;
			$lines .= <INNER>;
			$lines .= <INNER>;
		
        	if (defined $read2xf{$id}) {
	        	if ($seqLine =~ /^$bcSeq/) { 
	        	print OUTER "$header\n$lines $read2xf{$id}\n\n";
        		}else{
	       		 	print BCERR "$header\n$lines";
	       		 } 	
			}
		}
		close (INNER); close (OUTER); close (BCERR);
		
	} # unless combfile
	
	# dividing reads into gene files
	my $plusDir = $currDir.'plus';
	my %geneSeqHash;
	unless (-d $plusDir){
		system ("mkdir $plusDir");
		open (INNER, $combFile) || die print STDERR "cannot open $combFile for reading\n";
		print STDERR "dividing 	$combFile into gene files\n";
		my $line = <INNER>;
		my $fileName;
		my $index = 0;
		while (<INNER>){
			chomp;
			my $seq_line = $_;
			$seq_line =~ s/^.*?T{8,}//;
			$index++;
			$line = <INNER>; # +
			$line = <INNER>; # quality
			$line = <INNER>; # XF
			chomp $line;
			my @xf = split (':', $line);
			my $xf_name = $xf[-1];
			#print STDERR "seq=$seq_line,xf_name=$xf_name,xxx\n";
			$line = <INNER>; # empty line
			$line = <INNER>; # header
			unless ($xf_name =~ /no/){
				push @{$geneSeqHash{$xf_name}},$seq_line;
			}
			print STDERR ":" if $index % 100000 == 0;	
		}
		close (INNER); 
		
		print STDERR "\nprinting gene files\n";
		my @keys = keys %geneSeqHash;
		for (my $ki = 0; $ki < @keys; $ki++){
			my $gname = $keys[$ki];
			my @seqs = @{$geneSeqHash{$gname}};
			$fileName = $plusDir.'/'.$gname;
			open (XF, ">$fileName") || die print STDERR "cannot open $fileName for writing\n";
				for (my $vi = 0; $vi < @seqs; $vi++){
					print XF "$seqs[$vi]\n";
				}
			close (XF);
		}
		
	} # plus dir 
	print STDERR "finished processing for $idName\n";
}
close (IN);


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
