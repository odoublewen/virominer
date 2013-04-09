#! /usr/bin/perl -w
use strict ;
use Getopt::Long;
use CGI qw(:standard);

my $cgi = new CGI; 
my @parameters = $cgi->param;

my $servername = $ENV{'SERVER_NAME'};
my $docroot = $ENV{'DOCUMENT_ROOT'};
my $datadir = $ENV{'HTTP_REFERER'};

print "Content-type: text/plain\n\n";

#foreach (@parameters) {
#    print "$_ ";
#    print param("$_");
#    print "\n";
#}

my %headers = map { $_ => 1 } split(',',param("queries")) ;
my @files = split(',',url_param("data")) ;
my $datasetindex = int param("dataset");

my $header;
open FASTA, "$files[$datasetindex]";
while ( <FASTA> ) {
    if ( m/^>(\S+)/ ) {$header=$1}
    if (exists $headers{$header}) { print $_ }
}
