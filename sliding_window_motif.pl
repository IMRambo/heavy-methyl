#!/usr/bin/perl
#use File::chdir;
#use Cwd;
use strict;

#Step 2 of pipeline
#Last updated: January 21, 2016
#Written by Kenneth Hoadley
#Modified by Ian Rambo
# - - - - - H E A D E R - - - - - - - - - - - - - - - - -

#This script takes sequences from a FASTA file and runs a running average for GC content, CG and GGG motifs for each nt within the sequence.
#This script also catalgoues the occurence of N values in each sequence for downstream QC. Output is placed into the 01-Data folder


#Script modified for bacterial methyltransferase motifs
# - - - - - U S E R    V A R I A B L E S - - - - - - - -
my $dir = "/net/biohen29/data/imrambo";
my $infolder = "$dir/00-ExtractData";
my $outfolder = "$dir/01-ClusterData";
my $outfile1  = "GCcontent.txt";
my $outfile2  = "CpGmotif.txt";
my $outfile3  = "GATCmotif.txt";
my $outfile4 = "GCmotif.txt";
my $outfile5 = "CpGnorm.txt";
my $outfile6 = "GpCnorm.txt";
my $outfile7 = "observed-expected.txt";
my $outfile8 = "ATcontent.txt";
my $outfile9 = "GATCexpected.txt";
#my $outfile10 = "CCWGGmotif.txt";
#my $outfile11 = "GGWCCmotif.txt";



# - - - - - G L O B A L  V A R I A B L E S  - - - - - -

my %Genes;
my @GCdata;
my @CGdata;
my @GCmdata;
my @GATCdata;
my @Ndata;
my @CpGnorm;
my @gpcnorm;
my @observedexpected;
my @ATdata;
my @GATC_expdata;
#my @CCWGGdata;
#my @GGWCCdata;
#my @CCWGG_expdata;
#my @GGWCC_expdata;
my $sequencecount = 0;
my $CpG = 0;
#my $GATC = 0;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

print "\n\n . . . Generating data . . .\n\n";

opendir(DATA, $infolder) or die "\n\nNADA $infolder, idiot . . . . \n\n";
my @filelist = readdir(DATA);
foreach my $file(@filelist)
{	if ($file =~ m/^\w/)
	{	print "processing $file\n";
		my $proccessedfolder = $file;
		$proccessedfolder =~ s/\.txt//;
		#print "$proccessedfolder\n";
		system("mkdir $outfolder/$proccessedfolder");
		open(OUT1, ">".$outfolder."/".$proccessedfolder."/".$outfile1);
		open(OUT2, ">".$outfolder."/".$proccessedfolder."/".$outfile2);
		open(OUT3, ">".$outfolder."/".$proccessedfolder."/".$outfile3);
		open(OUT4, ">".$outfolder."/".$proccessedfolder."/".$outfile4);
		open(OUT5, ">".$outfolder."/".$proccessedfolder."/".$outfile5);
		open(OUT6, ">".$outfolder."/".$proccessedfolder."/".$outfile6);
		open(OUT7, ">".$outfolder."/".$proccessedfolder."/".$outfile7);
		open(OUT8, ">".$outfolder."/".$proccessedfolder."/".$outfile8);
		open(OUT9, ">".$outfolder."/".$proccessedfolder."/".$outfile9);
		#open(OUT11, ">".$outfolder."/".$proccessedfolder."/".$outfile10);
		#open(OUT12, ">".$outfolder."/".$proccessedfolder."/".$outfile11);
		
		&FASTAread($infolder."/"."$file");
		foreach my $id (keys %Genes)
		{   if ($Genes{$id}{'SIZE'} == 550)
		    {   $sequencecount += 1;
			#start position
			my $position = 25;
			while ($position <= 524)
			{   my $window = substr($Genes{$id}{'ntseq'}, $position-25, 51); #sets window size for running average (20)
			    my $winsize = length($window);
			    my $screen = 'N';
			    if ($window =~ m/$screen/i){
				push (@GCdata, '-');
				push (@CGdata, '-');
				push (@GCmdata, '-');
				push (@GATCdata, '-');
				push (@CpGnorm, '-');
				push (@gpcnorm, '-');
				push (@observedexpected, '-');
				push (@ATdata, '-');
				push (@GATC_expdata, '-');
				#push (@CCWGGdata, '-');
				#push (@GGWCCdata, '-');
				$position += 1;
				}
			    else
			    {
				    # %GC calculator
				    my $count = 0;
				    $count = $window =~ tr/GC/GC/;
				    my $GC = $count/$winsize;
				    push (@GCdata, $GC);
				    
				    # CpG motif calculator
				    my $CGcount = 0;
				    my $char = 'CG';
				    my $pos = index($window, $char, 0);
				    while ($pos > 0)
				    {   $CGcount += (1/(1 + abs(($winsize/2)-$pos))); #allows for values to be weighed depending on distance from position within window
					$pos = index($window, $char, $pos+1);
				    }
				    push (@CGdata, &RoundOff(($CGcount/$winsize), 20));
				   
				    my $CpGdensitycount = 0;
				    my $CpGchart = 'CG';
				    my $CpGposit = index($window, $CpGchart, 0);
				    while ($CpGposit > 0)
				    {	$CpGdensitycount += 1;
					$CpGposit = index($window, $CpGchart, $CpGposit+1);
				    }
				    $CpG = $CpGdensitycount;
				    #------------------------------------------
				    #Nucleotide calculator
				    my $Gcount = 0;
				    my $Ccount = 0;
				    my $Tcount = 0;
				    my $Acount = 0;
				    $Gcount = $window =~ (tr/G/G/) + 0.001;
				    $Ccount = $window =~ (tr/C/C/) + 0.001;
				    $Tcount = $window =~ (tr/T/T/) + 0.001;
				    $Acount = $window =~ (tr/A/A/) + 0.001;
				    #------------------------------------------
				    #CpG O/E
				    my $OE = (($CpG/($Ccount*$Gcount))*$winsize);
				    #print "$CpG\t$C\t$G\t$OE\n";
				    push (@observedexpected, &RoundOff($OE, 4));
				    #------------------------------------------
				    #GC motif calculator
				    my $GCcount = 0;
				    my $chart = 'GC';
				    my $posit = index($window, $chart, 0);
				    while ($posit > 0)
				    {   $GCcount += (1/(1 + abs(($winsize/2)-$posit))); #allows for values to be weighed depending on distance from position within window
					$posit = index($window, $chart, $posit+1);
				    }
				    push (@GCmdata, &RoundOff(($GCcount/$winsize), 20));
				    #------------------------------------------
				    #normalized CpG - denominator - 0.00001 prevents you from dividing by zero
				    my $CpGnormalized = (&RoundOff(($CGcount/$winsize), 10)/(0.00001 + &RoundOff($GC, 10))); 
				    push (@CpGnorm, &RoundOff($CpGnormalized, 20));
				    #------------------------------------------
				    #gpcnorm
				    my $gpcnormalized = (&RoundOff(($GCcount/$winsize), 10)/(0.00001 + &RoundOff($GC, 10)));
				    push (@gpcnorm, &RoundOff($gpcnormalized, 20));
				    #------------------------------------------
				    #GATC motif calculator
				    my $GATCcount = 0;
				    my $GATCchar = 'GATC';
				    my $GATCpos = index($window, $GATCchar, 0);
				    while ($GATCpos > 0)
				    {   $GATCcount += (1/(1 + abs(($winsize/2)-$GATCpos))); #allows for values to be weighed depending on distance from position within window;
					$GATCpos = index($window, $GATCchar, $GATCpos+1);
				    }
				    push (@GATCdata, &RoundOff(($GATCcount/$winsize), 20)); #normalizes $GATCcount with respect to the size of $winsize
				    
				#    my $GATCdensitycount = 0;
				#    my $GATCchart = 'GATC';
				#    my $GATCposit = index($window, $GATCchart, 0);
				#    while ($GATCposit > 0)
				#    {	$GATCdensitycount += 1;
				#	$GATCposit = index($window, $GATCchart, $GATCposit+1);
				#    }
				#    $GATC = $GATCdensitycount;
				    #------------------------------------------
				    #GATC expected
				    #Calculate dinucleotide frequencies
				    my $GA = () = $window =~ /GA/g;
				    my $GAfreq = $GA/length($window);
				    my $AT = () = $window =~ /AT/g;
				    my $ATfreq = $AT/length($window);
				    my $TC = () = $window =~ /TC/g;
				    my $TCfreq = $TC/length($window);
				    
				    #Expected frequency of GATC
				    my $GATCexp = (($GAfreq * $ATfreq * $TCfreq)/($Gcount * $Ccount * $Tcount * $Acount)) * ($winsize - 4 + 1);
				    push (@GATC_expdata, &RoundOff($GATCexp, 20));
				    #------------------------------------------
				    #GATC observed/expected - this has not been tested
				    #Split the nucleotide window sequence into codons
				    #my @codons = unpack("(A3)*", $window);
				    ##Counts for codon GAT
				    #my $countGAT = grep (/GAT/, @codons);
				    ##Counts for codon ATC
				    #my $countATC = grep (/ATC/, @codons);
				    #my $GAT = $countGAT + 0.001;
				    #my $ATC = $countATC + 0.001;
				    #
				    #my $GATCoe = (($GATC/($GAT * $ATC)) * (scalar @codons));
				    #------------------------------------------
				    #CCWGG motif calculator
				#    my $CCWGGcount = 0;
				#    my $CCWGGchar = "CC[AT]GG";
				#    my $CCWGGpos = index($window, $CCWGGchar, 0);
				#    while ($CCWGGpos > 0)
				#    {   $CCWGGcount += (1/(1 + abs(($winsize/2)-$GATCpos))); #allows for values to be weighed depending on distance from position within window;
				#	$CCWGGpos = index($window, $GATCchar, $GATCpos+1);
				#    }
				#    push (@CCWGGdata, &RoundOff(($CCWGGcount/$winsize), 20));
				#    #------------------------------------------
				#    #CCWGG expected
				#    my $CC = () = $window =~ /CC/g;
				#    my $CCfreq = $CC/length($window);
				#    my $CW = () = $window =~ /C[AT]/g;
				#    my $CWfreq = $CW/length($window);
				#    my $WG = () = $window =~ /[AT]G/g;
				#    my $WGfreq = $WG/length($window);
				#    my $GG = () = $window =~ /GG/g;
				#    my $GGfreq = $GG/length($window);
				#    
				#    my $CCWGGexp = (($CCfreq * $CWfreq * $WGfreq * $GGfreq)/($Gcount * $Ccount * $Tcount * $Acount)) * ($winsize - 5 + 1);
				#    #------------------------------------------
				#    #GGWCC expected
				#    my $GW = () = $window =~ /G[AT]/g;
				#    my $WC = () = $window =~ /[AT]C/g;
				#    my $GWfreq = $GW/length($window);
				#    my $WCfreq = $WC/length($window);
				#    
				#    my $GGWCCexp = (($GGfreq * $GWfreq * $WCfreq * $CCfreq)/($Gcount * $Ccount * $Tcount * $Acount)) * ($winsize - 5 + 1);
				    #------------------------------------------
				    #AT calculator
				    my $ATcount = 0;
				    my $ATchar = 'AT';
				    my $ATpos = index($window, $ATchar, 0);
				    while ($ATpos > 0)
				    {   $ATcount += (1/(1 + abs(($winsize/2)-$ATpos))); #allows for values to be weighed depending on distance from position within window
					$ATpos = index($window, $char, $ATpos+1);
				    }
				    push (@ATdata, &RoundOff(($ATcount/$winsize), 20));
				    #-------------------------------------------
				    
				    
				    
				#    my $ATdensitycount = 0;
				#    my $ATchart = 'AT';
				#    my $ATposit = index($window, $ATchart, 0);
				#    while ($ATposit > 0)
				#    {	$ATdensitycount += 1;
				#	$ATposit = index($window, $ATchart, $ATposit+1);
				#    }
				#    $AT = $ATdensitycount;
				    
				    #Move to the next base
				    $position += 1;
				    
			    }
			    
			}
	    
			 print OUT1 ">$Genes{$id}{'HEAD'}\t",join("\t",@GCdata),"\n";
			 print OUT2 ">$Genes{$id}{'HEAD'}\t",join("\t",@CGdata),"\n";
			 print OUT3 ">$Genes{$id}{'HEAD'}\t",join("\t",@GATCdata),"\n";
			 print OUT4 ">$Genes{$id}{'HEAD'}\t",join("\t",@GCmdata),"\n";
			 print OUT5 ">$Genes{$id}{'HEAD'}\t",join("\t",@CpGnorm),"\n";
			 print OUT6 ">$Genes{$id}{'HEAD'}\t",join("\t",@gpcnorm),"\n";
			 print OUT7 ">$Genes{$id}{'HEAD'}\t",join("\t",@observedexpected),"\n";
			 print OUT8 ">$Genes{$id}{'HEAD'}\t",join("\t",@ATdata),"\n";
			 print OUT9 ">$Genes{$id}{'HEAD'}\t",join("\t",@GATC_expdata),"\n";
			 @GCdata = ();   #reset for next sequence
			 @CGdata = ();   #reset for next sequence
			 @GCmdata = ();  #reset for next sequence
			 @GATCdata = ();  #reset for next sequence
			 @CpGnorm = ();  #reset for next sequence
			 @gpcnorm = ();  #reset for next sequence
			 @observedexpected = ();  #reset for next sequence
			 @ATdata = ();
			 @GATC_expdata = ();
			 $CpG = 0;
			 
		    }
			 
		}
		close (OUT1);
		close (OUT2);
		close (OUT3);
		close (OUT4);
		close (OUT5);
		close (OUT6);
		close (OUT7);
		close (OUT8);
		close (OUT9);
		print "there are $sequencecount in the file $file\n";
		$sequencecount = 0;
		
		#N-value QC check
		

	undef %Genes;
	}
}
	 
	 
	 
print "\n\n\n. . . . . .DONE. . . . . . \n\n\n ";


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


sub FASTAread
{	# 1. Load FIlE . . . . . . . . . .
	$/=">";                                     # set input break string
	my $infile = $_[0];
	open(IN, "<$infile") or die "\n\nNADA $infile you FOOL!!!\n\n";
	my @DATA = <IN>; close(IN); shift(@DATA);	
	# 2. Parse sequence data . . . . . . . . . . . . .
	my $unid = 10000;                           # string to generate unique ids
	foreach my $entry (@DATA)
	{	my @data = split('\n', $entry);
		my $seq = '';
		foreach my $i (1..$#data)
		{	$seq .= $data[$i];  }
		$seq =~ s/>//;
		#unless ($seq =~ /[RWYSKMBDHVN\.\-]/i)          # filter for non ATGC base pairs
		{	$Genes{$unid}{'HEAD'}    = $data[0];       # store header
			$Genes{$unid}{'ntseq'}   = uc($seq);       # store sequence
			$Genes{$unid}{'SIZE'}    = length($seq);   # store length
			$unid += 1;
		}
	}
	$/="\n";
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

sub RoundOff
{	my $num = shift(@_);
	my $decimals = shift(@_);
	my $roundnum = int(($num * 10**$decimals) + 0.5)/(10**$decimals);
	return $roundnum;
}

# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - -
#
#
