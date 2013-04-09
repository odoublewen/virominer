#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use File::Basename;
use Getopt::Long;
use Benchmark;
use POSIX qw(strftime);

my ( @s1, @q1, @s2, @q2, $file1, $file2, $basename, @basenames, $paired, $CMD, $counter,$rounded, $REDO, $outfile, @METAS, @FILTERS, $ADAPTERS, @BLAST) ;
my %filterhash;
my $OUTDIR = "virominer_out";
my $THREADS = 24 ;
my $BIN_PATH = dirname(__FILE__);
my $CWD = getcwd;
my $INDEX_PATH = "/home/owen/indices/" ;
my $BLAST_PATH = "" ;
my $filters = 'culex_genbank20121129,silva111_rrna' ;
my $metas = 'refseq_Eupath,refseq_Firmicutes,refseq_Fungi,refseq_OtherProks,refseq_ProteobacteriaABGP,refseq_viruses20130220' ;
my $metatrinity ;
my $blast = 'blastn,blastx';
my $BLASTNDB = 'nt';
my $BLASTXDB = 'nr';
my $call = "$0 @ARGV";
my ($TimeSampleStart, $TimeStageStart, $TimeStageEnd);

sub seqs2table {
    my $out = shift;
    open OUT, "$out";
    my $i;
    if ($paired) { for ($i = 0; $i < scalar @s1; $i++) { unless (exists $filterhash{$i}) {print OUT "$i\t$s1[$i]\t$q1[$i]\t$s2[$i]\t$q2[$i]\n";} } }
    else         { for ($i = 0; $i < scalar @s1; $i++) { unless (exists $filterhash{$i}) {print OUT "$i\t$s1[$i]\t$q1[$i]\n";} } }
    close OUT;
}

sub seqs2fasta {
    my $out = shift;
    my $counter = 0;
    open OUT, "$out";
    my $i;
    if ($paired) { for ($i = 0; $i < scalar @s1; $i++) { unless (exists $filterhash{$i}) {print OUT ">$i/1\n$s1[$i]\n>$i/2\n$s2[$i]\n"; $counter ++ ; } } }
    else         { for ($i = 0; $i < scalar @s1; $i++) { unless (exists $filterhash{$i}) {print OUT ">$i\n$s1[$i]\n"; $counter ++ ; } } }
    close OUT;
    return $counter;
}

GetOptions (
	    "redo" => \$REDO,
	    "out=s" => \$OUTDIR,
	    "threads=i" => \$THREADS,
	    "adapters=s" => \$ADAPTERS,
	    "blast=s" => \$blast,
	    "blastndb=s" => \$BLASTNDB,
	    "blastxdb=s" => \$BLASTXDB,
	    "indexpath=s" => \$INDEX_PATH,
	    "blastpath=s" => \$BLAST_PATH,
	    "filters=s" => \$filters,
	    "metagenomics=s" => \$metas,
	    "metatrinity=s" => \$metatrinity
    ) ;

@FILTERS = split ',', $filters ;
@METAS = split ',', $metas ;
@BLAST = split ',', $blast ;

my %TRINITY ;
if ( $metatrinity ) { %TRINITY = map { $_ => 1 } split ',', $metatrinity }
else                { $TRINITY{$METAS[-1]} = 1 }

## SANITY CHECKS GO HERE
if (not -d $OUTDIR) { mkdir $OUTDIR }
elsif (not $REDO)   { die "$OUTDIR exists.  Specify --redo to reuse previous output and regenerate reports, or choose a different output name with the --out option\n"}

foreach (@BLAST) { if ( $_ ne "blastn" and $_ ne "blastx" ) {die "blast option $_ is not blastn or blastx\n"}}
foreach ((@FILTERS, @METAS)) { if ( ! -s "$INDEX_PATH$_.1.bt2" ) { die "Cannot find bowtie2 index files at $INDEX_PATH$_.1.bt2\n" } }

open TIME, "> $OUTDIR\_timing.txt";
print TIME "$call\n";

$counter = 0;
for $file1 ( @ARGV ) {
    my ( %metahits, %counts );
    %filterhash = ();
    $counter ++ ;
    $TimeSampleStart = new Benchmark;

    print "# # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n";
    print "BEGINNING WORK ON FILE $file1 ($counter of " . scalar @ARGV . ")\n"; 
    if ( -s $file1 ) { print "Found file $file1"; }
    else             { print "File $file1 not found\n" ; next ; }
    $file2 = $file1 ;
    if ($file2 =~ s/_1\.fq/_2.fq/ and -s $file2 ) { $paired = 1; print " - Found file $file2 - running in paired end mode\n"; }
    else                                          { $paired = 0; print " - running in single end mode\n"; $file2 = ''; }

    $basename = $file1 ;
    if ($paired) {$basename =~ s/(.*\/)*(.*?)_1.fq(\.gz)*/$2/ }
    else         {$basename =~ s/(.*\/)*(.*?).fq(\.gz)*/$2/ }
    push @basenames, $basename;

    #### SEQUENCE CLEANUP AND FORMATTING
    print " STAGE 1: Sequence QC and formatting\n";
    $TimeStageStart = new Benchmark;
    (@s1, @q1, @s2, @q2) = () ;
    my $cutadaptstring;
    my $cutadaptoptions = "-q 25 -e 0 -O 8 -N";
    if ($ADAPTERS) {
	open PRIMERS, "$ADAPTERS";
	while (my $primer = <PRIMERS>) {
	    chomp $primer ;
	    $primer =~ s/.*\t(\w+)/$1/ ;
	    $primer = uc $primer;
	    $cutadaptstring .= "-b $primer ";
	    $primer =~ tr/ATGC/TACG/ ;
	    $primer = reverse $primer ;
	    $cutadaptstring .= "-b $primer ";
	}
	close PRIMERS;
	open IN, "cutadapt $cutadaptstring $cutadaptoptions --info-file $basename\_cutadapt_1.txt $file1 |" ;
    }
    elsif ($file1 =~ m/\.gz$/) {open IN, '-|', 'zcat', $file1 }
    else                       {open IN, $file1 }

    while ( <IN> ) {
	if    ($. % 4 == 0) { chomp ; push @q1, $_ ; $filterhash{scalar(@q1)-1} = 1 if (length($_) < 25) ;}
	elsif ($. % 2 == 0) { chomp ; push @s1, $_ ;}
    }
    close IN;

    if ($paired) {
	if ($ADAPTERS) { open IN, "cutadapt $cutadaptstring $cutadaptoptions --info-file $basename\_cutadapt_2.txt $file2 |" }
	elsif ($file2 =~ m/\.gz$/) {open IN, '-|', 'zcat', $file2 }
	else                    {open IN, $file2 }
	while ( <IN> ) {
	    if    ($. % 4 == 0) { chomp ; push @q2, $_ ; $filterhash{scalar(@q2)-1} = 1 if (length($_) < 25) ;}
	    elsif ($. % 2 == 0) { chomp ; push @s2, $_ ;}
	}
	close IN;
    }

    $counts{initial} = scalar @s1;

    $TimeStageEnd = new Benchmark;
    print TIME "Sample $basename\tInitialSeqs $counts{initial}\tStage 1 (INPUT)\t", strftime("%a, %d %b %Y %H:%M:%S %z", localtime(time())), "\n";
    print TIME "Sample $basename\tRejectedSeqs ", scalar keys %filterhash, "\tStage 1 (CUTADAPT)\t", timestr(timediff($TimeStageEnd, $TimeStageStart), 'all'), "\n";

    #### FILTERING
    if (scalar @FILTERS == 0) { print " STAGE 2: Skipping host filter -- all reads will be passed to metagenomics filters\n" }
    else { 
	print " STAGE 2: Host filter\n";
	$TimeStageStart = new Benchmark;
	for my $index ( @FILTERS ) {
	    $outfile = "$OUTDIR/$basename\_filter\_$index.txt";
	    $CMD = "| bowtie2 --sensitive-local -p $THREADS --12 - -x $INDEX_PATH$index | samtools view -S -F 4 - | cut -f 1 > $outfile" ;
	    if ($REDO and -s "$outfile") { print "  Skipping bowtie2 filter because file $outfile exists\n"}
	    else {
		seqs2table($CMD);
	    }
	    open HITS, "$OUTDIR/$basename\_filter\_$index.txt";
	    $counts{$index} = 0;
	  HITS: while (<HITS>) { 
	      chomp; 
	      if (exists $filterhash{$_}) {next HITS}
	      $filterhash{$_} = 1 ;
	      $counts{$index} ++;
	  }
	    close HITS;
	    $rounded = sprintf "%.2f", (100*($counts{$index}/$counts{initial}));
	    print " There were $counts{$index} hits to $index.  $counts{$index}/$counts{initial} = $rounded% of initial reads.\n";
	}
	$TimeStageEnd = new Benchmark;
	print TIME "Sample $basename\tFilteredSeqs ", scalar keys %filterhash, "\tStage 2 (HOSTFILT)\t", timestr(timediff($TimeStageEnd, $TimeStageStart), 'all'), "\n";
    }
    
    #### METAGENOMICS
    if (scalar @METAS == 0) { print " STAGE 3:  Skipping metagenomics filter -- all reads will be passed to Trinity\n" }
    else {
	print " STAGE 3: Metagenomics filter\n";
	$TimeStageStart = new Benchmark;
	for my $index ( @METAS ) {
	    $outfile = "$OUTDIR/$basename\_metagenomics_$index.txt";
	    $CMD = "| bowtie2 -p $THREADS --12 - -x $INDEX_PATH$index | samtools view -S -F 4 - | cut -f 1,3 > $outfile";
	    if ($REDO and -s "$outfile") { print "  Skipping bowtie2 metegenomics because file $outfile exists\n"}
	    else {
		seqs2table($CMD);
	    }
	    if (exists $TRINITY{$index}) {
		open HITS, "$OUTDIR/$basename\_metagenomics\_$index.txt";
		$counts{$index} = 0;
		while (<HITS>) { 
		    chomp; 
		    my @line = split "\t"; 
		    $counts{$index} ++ ;
		}
		close HITS;
	    }
	    else {
		open HITS, "$OUTDIR/$basename\_metagenomics\_$index.txt";
		$counts{$index} = 0;
	      HITS: while (<HITS>) { 
		    chomp; 
		    my @line = split "\t"; 
		    if (exists $filterhash{$line[0]}) { next HITS }
		    $filterhash{$line[0]} = 1 ;
		    $counts{$index} ++ ;
		}
		close HITS;
	    }
	    $rounded = sprintf "%.2f", (100*($counts{$index}/$counts{initial}));
	    print " There were $counts{$index} hits to $index.  $counts{$index}/$counts{initial} = $rounded% of initial reads.\n";
	    
	}
    
	## (THE SORT IS NECESSARY FOR KRONA)
	$CMD = "sort -n $OUTDIR/$basename\_metagenomics_*.txt > $OUTDIR/$basename\_metagenomics_sorted.txt";
	print "  CMD: $CMD\n";
	system $CMD;

	$TimeStageEnd = new Benchmark;
	print TIME "Sample $basename\tFilteredSeqs ", scalar keys %filterhash, "\tStage 3 (METAFILT)\t", timestr(timediff($TimeStageEnd, $TimeStageStart), 'all'), "\n";
	
	open COUNTS, ">$OUTDIR/$basename\_counts.txt";
	foreach (keys %counts) { print COUNTS "$_\t$counts{$_}\n"}
	close COUNTS;
    }

    ##### TRINITY ASSEMBLIES
    print " STAGE 4: Trinity Assembly\n";
    $TimeStageStart = new Benchmark;
    $counts{trinity} = seqs2fasta(">$OUTDIR/$basename\_trinityIN.fa");
    $outfile = "$OUTDIR/$basename\_trinity/Trinity.fasta";
    $CMD = "Trinity.pl --output $OUTDIR/$basename\_trinity --full_cleanup --min_contig_length 500 --seqType fa --JM 50G --single $OUTDIR/$basename\_trinityIN.fa --CPU $THREADS --inchworm_cpu $THREADS --bflyCalculateCPU";  ## --SS_lib_type RF
    if ($paired) { $CMD .= " --run_as_paired" }
    print "  CMD: $CMD\n";
    if ($REDO and -s "$outfile") { print "  Skipping Trinity assembly because file $outfile exists\n"}
    else {system $CMD}
    $TimeStageEnd = new Benchmark;
    print TIME "Sample $basename\tTrinitySeqs $counts{trinity}\tStage 4 (TRINITY)\t", timestr(timediff($TimeStageEnd, $TimeStageStart), 'all'), "\n";

    #### BLAST TRINITY ASSEMBLIES
    my $blastcounter = 1;
    my ( %blasthit, %annot ) ;
    foreach my $blaststage ( @BLAST ) {
	print " STAGE 5: BLAST of Trinity Contigs (blast program: $blaststage )\n";
	$TimeStageStart = new Benchmark;
	my ($blastOptions, $blastdb);
	if    ($blaststage eq "blastn") { $blastdb = $BLAST_PATH . $BLASTNDB; $blastOptions = "-task blastn -evalue 0.00001" }
	else                            { $blastdb = $BLAST_PATH . $BLASTXDB; $blastOptions = "-evalue 10" }
	my $seqEmitter;
	if ($blastcounter == 1) { $seqEmitter = "cat" }
	else                    { $seqEmitter = "$BIN_PATH/fastaFilter.pl $outfile" }

	$outfile = "$OUTDIR/$basename\_darkmatter_$blaststage.txt";

	$CMD = "$seqEmitter $OUTDIR/$basename\_trinity/Trinity.fasta | parallel --block 2k --recstart '>' --pipe -j $THREADS \"$blaststage -db $blastdb $blastOptions -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos' -max_target_seqs 5 -query -\" > $outfile";
	print "  CMD: $CMD\n";
	if ($REDO and -s "$outfile") { print "  Skipping BLAST search because file $outfile exists\n"}
	else { system $CMD } 

	## FOR EACH BLAST ALIGNMENT, EXTRACT ACCESSION NUMBER OF TOP (ie, FIRST) HIT, SEND TO FILE...
	open BLAST, "$outfile";
	open TOPHITS, "> $OUTDIR/$basename\_darkmatter_$blaststage\_tophitaccessions.txt";
	while (<BLAST>) {
	    chomp;
	    my @line = split "\t";
	    if (exists $blasthit{$line[0]}) {next}
	    my @idstring = split /\|/, $line[1] ;
	    $blasthit{$line[0]}{accession} = $idstring[3];
	    $blasthit{$line[0]}{pident} = $line[2];
	    $blasthit{$line[0]}{length} = $line[3];
	    $blasthit{$line[0]}{ppos} = $line[12];
	    print TOPHITS "$line[1]\n";
	}
	close TOPHITS;
	close BLAST;

	#### ... THEN MAKE ANNOTATION DUMP FILE AND READ IT IN...
	$CMD = "blastdbcmd -db $blastdb -outfmt '%a %t (%l)' -entry_batch $OUTDIR/$basename\_darkmatter_$blaststage\_tophitaccessions.txt > $OUTDIR/$basename\_darkmatter_$blaststage\_annotations.txt";
	system $CMD;

	open ANNOT, "$OUTDIR/$basename\_darkmatter_$blaststage\_annotations.txt";
	while (<ANNOT>) {
	    chomp;
	    m/(\S+) (.+)/ ;
	    $annot{$1} = $2;
	}
	close ANNOT;

	$TimeStageEnd = new Benchmark;
	print TIME "Sample $basename\t\tStage 5 BLAST ($blaststage)\t", timestr(timediff($TimeStageEnd, $TimeStageStart), 'all'), "\n";

	$blastcounter ++ ;
    }

    ### ... AND FINALLY ADD ANNOTATIONS TO TRINITY CONTIGS
    $TimeStageStart = new Benchmark;
    
    open FASTA, "$OUTDIR/$basename\_trinity/Trinity.fasta" ;
    open FASTA2, "> $OUTDIR/$basename\_Trinity_annotated.fasta" ;
    while (<FASTA>) {
	chomp;
	if ( m/>(.+) len=(\d+) path=.+/ ) {
	    if (exists $blasthit{$1}) { print FASTA2 ">$1 $2bp tophit: $blasthit{$1}{accession} alength=$blasthit{$1}{length} pident=$blasthit{$1}{pident} ppos=$blasthit{$1}{ppos} $annot{$blasthit{$1}{accession}}\n" }
	    else { print FASTA2 ">$1 $2bp no blast hits\n" }
	}
	else { print FASTA2 "$_\n" }
    }
    close FASTA;
    close FASTA2;

    $CMD = "fastaLengths.pl $OUTDIR/$basename\_Trinity_annotated.fasta > $OUTDIR/$basename\_Trinity_lengths.txt";
    print "  CMD: $CMD\n";
    system $CMD;

    $CMD = "cat $OUTDIR/$basename\_darkmatter_blast[nx].txt > $OUTDIR/$basename\_darkmatter_blast.txt";
    print "  CMD: $CMD\n";
    system $CMD;

    $TimeStageEnd = new Benchmark;
    print TIME "Sample $basename\t\tStage 6 FASTA ANNOTATION\t", timestr(timediff($TimeStageEnd, $TimeStageStart), 'all'), "\n";

}

print TIME "Sample $basename\t\tStages 1-6\t", timestr(timediff($TimeStageEnd, $TimeSampleStart), 'all'), "\n";


#### BUILDING KRONA PLOTS
print "# # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n";
my (@metabasenames, $metabasenamestring);
unless (scalar @METAS == 0) {
    print "BUILDING KRONA SUMMARY OF METAGENOMIC HITS\n";
    $TimeStageStart = new Benchmark;
    @metabasenames = map { "$OUTDIR/$_"."_metagenomics_sorted.txt"} @basenames ;
    $metabasenamestring = join ' ', @metabasenames ;
    $CMD = "ktImportBLAST -o $OUTDIR\_krona_metagenomic.html $metabasenamestring";
    print "CMD: $CMD\n";
    system $CMD;
}
print "BUILDING KRONA SUMMARY OF BLAST HITS\n";
@metabasenames = map { "$OUTDIR/$_\_darkmatter_blast.txt:$OUTDIR/$_\_Trinity_lengths.txt"} @basenames ;
$metabasenamestring = join ' ', @metabasenames ;
my @urls = map { "$CWD/$OUTDIR/$_\_Trinity_annotated.fasta"} @basenames ;
my $urlstring = "/cgi-bin/kronaQuery.pl?data=" . join ',', @urls ;
$CMD = "ktImportBLAST -o $OUTDIR\_krona_darkmatterblast.html -x 75 -y 256 -qp $urlstring $metabasenamestring";
print "CMD: $CMD\n";
system $CMD;

$TimeStageEnd = new Benchmark;
print TIME "Sample ALL\t\tStage 7 (KRONA)\t", timestr(timediff($TimeStageEnd, $TimeStageStart), 'all'), "\n";

close TIME;








=header comment

code scraps go here

=cut

