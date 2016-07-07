#!/share/apps/perl-5.14.2/bin/perl


use strict;

my $inputVCF = $ARGV[0];
my $outputVCF = $ARGV[1];
my $threshold = $ARGV[2];


my $variantExclusion = "/cluster/project8/vyp/exome_sequencing_multisamples/mainset/support/exclusion_lists/variant_exclusion.tab";
my %excluHash = ();

open (exclu, " < $variantExclusion");
while (<exclu>) {
    chomp;
    $excluHash { $_ } = 'yes';
    print "Finding variant $_ to be excluded\n";
}
close (exclu);


open (INP, " < $inputVCF") or die "Cannot open $inputVCF\n";
open (OUT, " > $outputVCF");






#GT:AD:DP:GQ:PL

while (<INP>) {

    if ($_ =~ /^\#/) {print OUT $_;} else {
	chomp $_;
	my @spl = split ('\t', $_);
	my $signature = $spl[0]."_".$spl[1]."_".$spl[3]."_".$spl[4];
	
	if (exists $excluHash { $signature }) {
	    print "Excluding variant $signature \n";
	    next;
	}
	

	print OUT $spl[0];
	for (my $i = 1; $i != 9; $i++) {print OUT "\t".$spl[ $i ];}
	
	my $ezSNP = 0;
	if ( (length($spl[ 3 ]) == 1) && (length($spl[ 4 ]) == 1) && ($spl[ 3 ] =~ /^[A-Z]/) && ($spl[ 4 ] =~ /^[A-Z]/) ) {
	    $ezSNP = 1;
	    #print $spl[ 3 ].$spl[4]."\n";
	}

	for (my $i = 9; $i <= $#spl; $i++) {
	    my @spl2 = split(':', $spl[ $i ]);

	    my $good = 1; ## by default we accept
	    my $locth = $threshold;

	    my $geno = $spl2[ 0 ];
	    my $AD = $spl2[ 1 ];
	    my $DP = $spl2[ 2 ];
	    my $GQ = $spl2[ 3 ];
	    my $PL = $spl2[ 4 ];
	    
	    #print $#spl2."   ".$spl[ $i ]."  ".$spl[8]."\n";
	    if ($#spl2 == 6) {
		$PL = $spl2[ 6 ]; ### some bug correction to deal with these odd PGT and PID fields
		#print "Warning ".$spl[8]."\t".$PL."\n";
	    }


	    if ($geno eq "0\/0") {
		if ($GQ <= $locth) {$good = 0;}
	    } else {
		my @PLs = split(',', $PL);
		my $PL0 = $PLs[ 0 ];
		if ($PL0 <= $locth) {$good = 0;}
	    }

	    #if ($geno eq "1\/1") {$locth = 9;}  ## for ALT/ALT calls, a 15 QUAL will do in any case
	    #if ($geno eq "2\/2") {$locth = 9;}  ## for ALT/ALT calls, a 15 QUAL will do in any case
	    #if ($geno eq "3\/3") {$locth = 9;}  ## for ALT/ALT calls, a 15 QUAL will do in any case

	    if ($ezSNP) {  ##if nice clean di-allelic SNP we add an extra filter
		if ($geno eq "0\/1") {
		    my @counts = split(',', $AD);
		    my $total = $counts[0] + $counts[1];
		    if ($counts[1] > $counts[0] + 1) {$locth = 9;}  ##if the alternate count is higher then we can be relax, because the alternative is probably 

		    if ($total >= 10) {
			my $ratio = $counts[1]/$total;
			my $limit = 0.18;
			if ($ratio <= $limit) {$good = 0;}
		    }
		}
	    }



	    if ($good) {
		#print OUT "\t".$spl[ $i ];
		print  OUT "\t".$spl2[0].":".$spl2[1].":".$spl2[ 2 ].":".$spl2[ 3 ].":".$PL;
	    } else { ## make the data missing
		#print  OUT "\t./.:".$spl2[1].":".$spl2[ 2 ].":".$spl2[ 3 ].":".$spl2[ 4 ];
		print  OUT "\t./.:".$spl2[1].":".$spl2[ 2 ].":".$spl2[ 3 ].":".$PL;
	    }
	    
	}

	
	print OUT "\n";
    }
    

}

close (INP);
close (OUT);

