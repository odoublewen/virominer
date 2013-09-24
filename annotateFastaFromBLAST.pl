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
	    $blasthit{$line[0]}{plus} = ( ($line[7]-$line[6]) * ($line[9]-$line[8]) ) > 0 ;
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
    
    open FASTA, "$OUTDIR/$basename\_trinity.Trinity.fasta" ;
    open FASTA2, "> $OUTDIR/$basename\_Trinity_annotated.fasta" ;
    while (<FASTA>) {
	chomp;
	if ( m/>(.+) len=(\d+) path=.+/ ) {
	    if (exists $blasthit{$1}) { 
		print FASTA2 ">$1 $2bp tophit: $blasthit{$1}{accession} alength=$blasthit{$1}{length} pident=$blasthit{$1}{pident} ppos=$blasthit{$1}{ppos} $annot{$blasthit{$1}{accession}}" ;
		if (! $blasthit{$1}{plus} ) { print FASTA2 " RC" }
		print FASTA2 "\n";
	    }
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
