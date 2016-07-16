#!/share/apps/perl-5.14.2/bin/perl

use lib '/cluster/project8/vyp/vincent/libraries/perl/Text-CSV-1.21/blib/lib';
use Text::CSV;

use strict;

my $inputTable = $ARGV[0];
my $output = $ARGV[1];
my $chromosomeArg = $ARGV[2];

print "Output file: $output\n";

my $csv = Text::CSV->new({ sep_char => ',' });
my $outcsv = Text::CSV->new({ sep_char => ',' });

open (INP, " < $inputTable") or die "Cannot open $inputTable";

my $colnameFile = $output."_snpStats/colname_chr".$chromosomeArg.".tab";


### should not need the names below really
my $callFile = $output."_snpStats/calls.tab";
my $depthFile = $output."_snpStats/depth.tab";
my $rownameFile = $output."_snpStats/rowname.tab";
my $annotationFile = $output."_snpStats/annotations.csv";


open (COL, " > $colnameFile");

my ($Chr, $start, $ref, $alt) = (0,0,0,0);

my $currentChr = -1;

#### header line
$_ = <INP>;
if ($csv->parse($_)) {
    
    my @headers = $csv->fields();
    print $#headers."\n";
   

    my @newline;

    for (my $i = 0; $i <= 29; $i++) {
	$newline[ $i ] = $headers[ $i ];
	if ($headers[ $i ] eq "Chr") { $Chr = $i; };
	if ($headers[ $i ] eq "Start") { $start = $i;};
	if ($headers[ $i ] eq "Ref") { $ref = $i;};
	if ($headers[ $i ] eq "Obs") { $alt = $i;};
    }
    
    ##### Now print the annotation headers
    $outcsv->combine(@newline);
    my $headers = $outcsv->string."\n";
    
	######### Now looking at the actual data
    for (my $i = 30; $i <= $#headers; $i++) {print COL $headers[ $i ]."\n";}
    close (COL);

    my %hash_sig = ();

    my $row = 0;
    $_ = <INP>;
    do {
	if ($csv->parse($_)) {
	    $row++;
	    my @fields = $csv->fields();
	    my $chromosome = $fields[ $Chr ];
	    my $signature = $chromosome."_".$fields[ $start ]."_".$fields[ $ref ]."_".$fields[ $alt ];

	    if ( $chromosome != $currentChr) {
		print "Change of chromosome\n";
		$callFile = $output."_snpStats/calls_chr".$chromosome.".tab";
		$depthFile = $output."_snpStats/depth_chr".$chromosome.".tab";
		$rownameFile = $output."_snpStats/rowname_chr".$chromosome.".tab";
		$annotationFile = $output."_snpStats/annotations_chr".$chromosome.".csv";
		
		if ($currentChr != -1) {
		    close(OUT);
		    close(DEPTH);
		    close(ROW);
		    close(ANN);
		}
		
		print "Now writing the calls to $callFile\n";
		open (OUT, " > $callFile") or die $callFile;
		open (DEPTH, " > $depthFile") or die;
		open (ROW, " > $rownameFile") or die;
		open (ANN, " > $annotationFile") or die;
		
		print ANN $headers;
		$currentChr = $chromosome;
	    }

	    ####### Now dealing with the cleaner signature issue
	    my $cleanRef = $fields[ $ref ];
	    my $cleanAlt = $fields[ $alt ];	    
	    my $slen = 10;
	    my $cleanRef2 = substr($cleanRef, 0, $slen);
	    my $cleanAlt2 = substr($cleanAlt, 0, $slen);

	    my $clean_signature = $chromosome."_".$fields[ $start ]."_".$cleanRef2."_".$cleanAlt2;

	    my $extra = 1;
	    while (exists  $hash_sig{ $clean_signature }) {
		$extra++;
		$clean_signature = $chromosome."_".$fields[ $start ]."_".$cleanRef2."_".$cleanAlt2.".".$extra;
	    }
	    
	    if (exists  $hash_sig{ $clean_signature }) {
		print "Bug with signature\n";
		exit();
	    }

	    $hash_sig{ $clean_signature } = "OK";
	    #################

	    print ROW $clean_signature."\n";
	    print $row."\t".$signature."\t".$clean_signature."\n";
	    
	    ########## this below prints the annotation data
	    my @nnewline;	    
	    for (my $i = 0; $i <= 29; $i++)  {$nnewline[ $i ] = $fields[ $i ];}
	    $outcsv->combine(@nnewline);
	    print ANN $outcsv->string, "\n";

	    
	    my @spl = split(':', $fields[ 30 ]);
	    if ($spl[ 0 ] eq ".|.") {print OUT "NA";}
	    if ($spl[ 0 ] eq "0|0") {print OUT "0";}
	    if ($spl[ 0 ] eq "0|1") {print OUT "1";}
	    if ($spl[ 0 ] eq "1|1") {print OUT "2";}
	    
	    print DEPTH $spl[1];
	    
	    
	    for (my $i = 31; $i <= $#fields; $i++)  {
		my @spl = split(':', $fields[ $i ]);
		if ($spl[ 0 ] eq ".|.") {print OUT "\tNA";}
		if ($spl[ 0 ] eq "0|0") {print OUT "\t0";}
		if ($spl[ 0 ] eq "0|1") {print OUT "\t1";}
		if ($spl[ 0 ] eq "1|1") {print OUT "\t2";}

		print DEPTH "\t".$spl[1];
	    }
	    print OUT "\n";
	    print DEPTH "\n";
	}
	
	
    } while ( <INP> );

}


close(OUT);
close(DEPTH);
close(ROW);
close(ANN);
