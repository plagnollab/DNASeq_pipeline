#!/share/apps/perl-5.14.2/bin/perl

use lib '/cluster/project8/vyp/vincent/libraries/perl/Text-CSV-1.21/blib/lib';
use Text::CSV;

use strict;

print "Number of arguments $#ARGV\n";
if ($#ARGV != 6) {die "There should be 6 arguments as follows:\n .exec inputTable controlKeywords caseKeywords outputTable sampleList excludedControlList chooseExternalControlsRandom";}

##default option here
my $chooseExternalControlsRandom = "yes"; ### if set to "no" the external controls will be the non cases and non controls. if set to "yes" the external controls are randomly samples from the controls.

my $inputTable = $ARGV[0];
my $controlKeywords = $ARGV[1];
my $caseKeywords = $ARGV[2];
my $outputTable = $ARGV[3];
my $sampleList = $ARGV[4];
my $excludedControlList = $ARGV[5];
my $chooseExternalControlsRandom = $ARGV[6];


#if ( $excludedControlList eq "none" ) {$chooseExternalControlsNotRandom = 1;}

########
my %excludedControl = ();
if ($excludedControlList ne "none") {
    print "Using $excludedControlList \n";
    open (exclu, "< $excludedControlList") or die "Cannot open $excludedControlList";
    while (<exclu>) {
	chomp $_;
	print "Sample $_ cannot be used as control\n";
	$excludedControl { $_ } = 'yes';
    }
    close(exclu);
}

############## look at control keywords


open (DB, " < $controlKeywords");
my @keywords;
my $count =0;
while (<DB>) {
    chomp;
    if (length($_) > 0) {
	$keywords[ $count ] = $_;
	$count++;
    }
}
close (DB);

for (my $i = 0; $i <= $#keywords; $i++) {chomp $keywords[$i]; print "Control keyword: ".$keywords[$i]."\n";}


################
open (DB2, " < $caseKeywords");
my @cakeywords;
my $count2 =0;
while (<DB2>) {
    chomp;
    if (length($_) > 0) {
	$cakeywords[ $count2 ] = $_;
	$count2++;
    }
}
close (DB2);

for (my $i = 0; $i <= $#cakeywords; $i++) {chomp $cakeywords[$i]; print "Case keyword: ".$cakeywords[$i]."\n";}


############## open input file

my $csv = Text::CSV->new({ sep_char => ',' });


open (OUTS, " > $sampleList");
open (INP, " < $inputTable");

print OUTS "sample control externalControl case\n";

my $_ = <INP>;
my@spl = split(',', $_);

my %cafield = ();
my %cfield = ();
my %extfield = ();

#my @controlfields;
#my @casefields;

my $count = 0;
my $cacount = 0;

my $nsamples = 0;
my $ncontrols = 0;
my $ncases = 0;
my $nextcontrols = 0;
my @samples; 




##### First we parse the header line
if ($csv->parse($_)) {
    @samples = $csv->fields();
    print "Number of fields $#samples\n";

    for (my $i = 30; $i <= $#samples; $i++) {
	my $name = $samples[ $i ];
	if (exists $excludedControl { $name } ) {
	    print "This individual is excluded as a control: $name, perhaps because of relatedness\n";
	} else {
	    for (my $j = 0; $j <= $#keywords; $j++) {
		if ($name =~ /$keywords[$j]/) {
		    $count = $count + 1;
		    if ( ($count % 4 == 0) && ($chooseExternalControlsRandom eq "yes") ) {  ##there is a flag here to select the external controls randomly or not
			$extfield{ $i } = 'yes';
			$nextcontrols = $nextcontrols + 1;
		    } else {
			$cfield{ $i } = 'yes';
			$ncontrols = $ncontrols + 1;
			$j = $#keywords + 1;  ##to quit the loop
		    }
		} 
	    }
	}

	for (my $j = 0; $j <= $#cakeywords; $j++) {
	    if ($samples[ $i ] =~ /$cakeywords[$j]/) {
		
		if (exists $cfield{ $i }) {print $samples[ $i ]." is a case and a control together".$cakeywords[$j]."\n";exit;}
		$cacount = $cacount + 1;
		$cafield{ $i } = 'yes'; 
		$ncases = $ncases + 1;
		$j = $#cakeywords + 1;  ##to quit the loop
	    }
	}
	
	if ($chooseExternalControlsRandom eq "no" ) {
	    if (  (! exists ( $cafield{ $i } ) ) && (! exists ( $cfield{ $i } )) ) { 
		if (exists $excludedControl { $name } ) {
		    print "This individual is excluded as an external control: $name, perhaps because of relatedness\n";
		} else {
		    $extfield{ $i } = 'yes';
		    $nextcontrols = $nextcontrols + 1;
		    print "External control: ".$name."\n";
		}
	    }
	}
    


	if ($i > 29) {
        if (! exists $cfield{ $i }) {print OUTS $samples[ $i ]."\tno";} else {print OUTS $samples[ $i ]."\tyes";}  ##control flag
        if (! exists $extfield{ $i }) {print OUTS "\tno"; } else {print OUTS "\tyes";} 	    ###external control flag
        if (! exists $cafield{ $i }) {print OUTS "\tno\n";} else {print OUTS "\tyes\n";}  ###case flag
	}
    }

}
print "Number of controls identified: $ncontrols , external controls $nextcontrols and nb of cases $ncases\n";
close (OUTS);
print "Sample list in $sampleList\n";


my $nvariants = 0;
open (OUT, " > $outputTable");
do {

    if ($csv->parse($_)) {
	my @fields = $csv->fields();
	my @newline;
	my $newf = 0;
	my $totalCallControls = 0;
	my $totalCallExtControls = 0;
	my $nonRefCallControls = 0;
	my $nonRefCallExtControls = 0;
	my $MAFControls = 0;

	my $hetnames = "";
	my $homnames = "";
	my $nhets = 0;
	my $nhoms = 0;
	my $nExthets = 0;
	my $nExthoms = 0;

	####### the non sample columns
	for (my $i = 0; $i <= 29; $i++) {
	    $newline[ $newf ] = $fields[ $i ];
	    $newf++;
	}

	#### now the samples
	for (my $i = 30; $i <= $#fields; $i++) {
	    
	    if ( (exists $extfield{ $i } ) || (exists $cfield{ $i })) {  ##we have a control
		
		my @spl  = split(':', $fields[ $i ]);
		if ( ($spl[0] !~ /^\./ ) && ($spl[1] >= 5)) { ##Now we have a control call: depth greater than 5 and non missing
		    
		    if (exists $cfield{ $i }) {
			$totalCallControls++;
			if ($spl[0] eq "0|1" ) {$nonRefCallControls++; $hetnames = $hetnames.";".$samples[ $i ]; $nhets++;}
			if ($spl[0] eq "1|1" ) {$nonRefCallControls = $nonRefCallControls + 2; $homnames = $homnames.";".$samples[ $i ]; $nExthoms++;}
		    }
		    
		    
		    if (exists $extfield{ $i }) {
			$totalCallExtControls++;
			if ($spl[0] eq "0|1" ) {$nonRefCallExtControls++; $hetnames = $hetnames.";".$samples[ $i ]; $nExthets++;}
			if ($spl[0] eq "1|1" ) {$nonRefCallExtControls = $nonRefCallExtControls + 2; $homnames = $homnames.";".$samples[ $i ]; $nExthoms++;}
		    }
		}
	    } 

	    if (exists $cafield{ $i }) {  ##we have a case then
		#print $fields[ $i ]."\n";  #debug
		$newline[ $newf ] = $fields[ $i ];
		$newf++;
	    }
	}
	#exit; #debug
	my $nTotHets = $nhets + $nExthets;
	my $nTotHoms = $nhoms + $nExthoms;
	
	$hetnames = $nTotHets.$hetnames;
	$homnames = $nTotHoms.$homnames;
	
	if ($nvariants > 0) { ##non header line
	    $newline[ $newf ] = $totalCallControls; $newf++;
	    $newline[ $newf ] = $nonRefCallControls; $newf++;
	    if ($totalCallControls > 0) {$newline[ $newf ] = $nonRefCallControls/(2*$totalCallControls);} else {$newline[ $newf ] = 0;}; $newf++;
	    
	    if ($nonRefCallControls + $nonRefCallExtControls < 7) {
		$newline[ $newf ] = $hetnames; $newf++;
		$newline[ $newf ] = $homnames; $newf++;
	    } else {
		$newline[ $newf ] = $nhets; $newf++;
		$newline[ $newf ] = $nhoms; $newf++;
	    }
	    
	    $newline[ $newf ] = $totalCallExtControls; $newf++;
	    $newline[ $newf ] = $nonRefCallExtControls; $newf++;

	    if ($totalCallExtControls > 0) {$newline[ $newf ] = $nonRefCallExtControls/(2*$totalCallExtControls);} else {$newline[ $newf ] = 0;}; $newf++;

	} else { #header line
	    $newline[ $newf ] = "total.calls.controls"; $newf++;
	    $newline[ $newf ] = "non.ref.calls.controls"; $newf++;
	    $newline[ $newf ] = "freq.controls"; $newf++;
	    $newline[ $newf ] = "HetNames"; $newf++;
	    $newline[ $newf ] = "HomNames"; $newf++;
	    $newline[ $newf ] = "total.calls.external.controls"; $newf++;
	    $newline[ $newf ] = "non.ref.calls.external.controls"; $newf++;
	    $newline[ $newf ] = "freq.external.controls"; $newf++;
	}

	if ($nvariants % 1000 == 0) {print "Variant $nvariants Number of fields in squeezed array $newf\n";}
	
	

	if ($csv->combine(@newline)) {
	    my $string = $csv->string;
	    print OUT $string, "\n";
	    $nvariants++;
    } else {
	my $err = $csv->error_input;
	print "combine() failed on argument: ", $err, "\n";
    }
	
	
    } else {
	my $err = $csv->error_input;
	print "Failed to parse line: $err";
	exit;
    }
    
    #if ($nvariants == 100) {close(INP);}
} while (<INP> );

close(INP);
close (OUT);

