#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my ($fastafile, $outfile, $verbose) ;

GetOptions ("out=s" => \$outfile , 
	    "fasta=s" => \$fastafile ,
	    "verbose" => \$verbose
);

sub rev_complement_IUPAC {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
        return $revcomp;
}

my (%blasthit, %annot) ;

foreach my $blastfile (@ARGV) {
    ( my $blastdb = $blastfile ) =~ s/.+\.(\w+)/$1/ ;

    if ($verbose) { print "Parsing blast results from $blastfile, looking up hits in blastdb $blastdb\n" }

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
my $header;
my $length;
my $sequence;
while (<FASTA>) {
    chomp;
    if ( m/>(.+) len=(\d+) path=.+/ ) {
	# "if $header" checks to see if this is NOT the very first sequence in the file.
	if ($header) {
	    if (exists $blasthit{$header}) { 
		print FASTA2 ">$header " . $length ."bp tophit: $blasthit{$header}{accession} alength=$blasthit{$header}{length} pident=$blasthit{$header}{pident} ppos=$blasthit{$header}{ppos}";

		if (exists $annot{$blasthit{$header}{accession}} ) { 
		    print FASTA2 " $annot{$blasthit{$header}{accession}}" ;
		}

		if (! $blasthit{$header}{plus} ) {
		    print FASTA2 " RC";
		    $sequence = rev_complement_IUPAC $sequence;
		}
		print FASTA2 "\n$sequence\n" ;
	    }
	    else { print FASTA2 ">$header ".$length."bp no blast hits\n$sequence\n" }
	    $header = '';
	    $length = '';
	    $sequence = '';
	}
	else {
	    $header = $1;
	    $length = $2;
	}
    }
    else {
	$sequence .= $_ ;
    }
}

if (exists $blasthit{$header}) { 
    print FASTA2 ">$header ".$length."bp tophit: $blasthit{$header}{accession} alength=$blasthit{$header}{length} pident=$blasthit{$header}{pident} ppos=$blasthit{$header}{ppos}";
    
    if (exists $annot{$blasthit{$header}{accession}} ) { 
	print FASTA2 " $annot{$blasthit{$header}{accession}}" ;
    }
    
    if (! $blasthit{$header}{plus} ) { 
	$sequence = rev_complement_IUPAC $sequence;
    }
    print FASTA2 "\n$sequence\n" ;
}
else { print FASTA2 ">$header ".$length."bp no blast hits\n$sequence\n" }


close FASTA;
close FASTA2;

