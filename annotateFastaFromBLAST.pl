#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my ($fastafile, $outfile) ;

GetOptions ("out=s" => \$outfile , 
	    "fasta=s" => \$fastafile );

my (%blasthit, %annot) ;

foreach my $blastfile (@ARGV) {
    ( my $blastdb = $blastfile ) =~ s/.+\.(\w+)/$1/ ;

    print "Parsing blast results from $blastfile, looking up hits in blastdb $blastdb\n";

    ## FOR EACH BLAST ALIGNMENT, EXTRACT ACCESSION NUMBER OF TOP (ie, FIRST) HIT, SEND TO FILE...
    open BLAST, "$blastfile";
    open TOPHITS, "> $outfile\_tophitaccessions.txt";

    while (<BLAST>) {
	chomp;
	my @line = split "\t";
	if (exists $blasthit{$line[0]}) {next}
	my @idstring = split /\|/, $line[1] ;
	$blasthit{$line[0]}{accession} = $idstring[3];
	$blasthit{$line[0]}{pident} = $line[2];
	$blasthit{$line[0]}{length} = $line[3];
	$blasthit{$line[0]}{ppos} = $line[12];
	$blasthit{$line[0]}{plus} = ( ($line[7]-$line[6]) * ($line[9]-$line[8]) ) > 0 ;
	print TOPHITS "$line[1]\n";
    }
    close TOPHITS;
    close BLAST;

    ## ... THEN MAKE ANNOTATION DUMP FILE AND READ IT IN...
    my $CMD = "blastdbcmd -db $blastdb -outfmt '%a %t (%l)' -entry_batch $outfile\_tophitaccessions.txt > $outfile\_annotations.txt";
    system $CMD;
    
    open ANNOT, "$outfile\_annotations.txt";
    while (<ANNOT>) {
	chomp;
	m/(\S+) (.+)/ ;
	$annot{$1} = $2;
    }
    close ANNOT;
}

## ... AND FINALLY ADD ANNOTATIONS TO TRINITY CONTIGS
open FASTA, "$fastafile" ;
open FASTA2, "> $outfile\_annotated.fasta" ;
while (<FASTA>) {
    chomp;
    if ( m/>(.+) len=(\d+) path=.+/ ) {
	if (exists $blasthit{$1}) { 
	    print FASTA2 ">$1 $2bp tophit: $blasthit{$1}{accession} alength=$blasthit{$1}{length} pident=$blasthit{$1}{pident} ppos=$blasthit{$1}{ppos}";

	    if (exists $annot{$blasthit{$1}{accession}} ) {
		print FASTA2 " $annot{$blasthit{$1}{accession}}" ;
	    }
	    if (! $blasthit{$1}{plus} ) { print FASTA2 " RC" }
	    print FASTA2 "\n";
	}
	else { print FASTA2 ">$1 $2bp no blast hits\n" }
    }
    else { print FASTA2 "$_\n" }
}
close FASTA;
close FASTA2;

