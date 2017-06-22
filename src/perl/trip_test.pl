#!/usr/bin/perl

my %BCs;
my $fastq = 'SURFdrive/TRIP/trip0.5/normalization_reads.fastq';
my $pat1 = 'TACACAACTCGAG';
my $pat2 = 'TGATCCTGCATAC';
my $type = 'norm';
my $ind_len = 26;
my bc_len = 16;
my @stats_norm = bc_extract_exp1 ($fastq, "norm" , $ind_len , $bc_len , $pat1 , $pat2 , \%BCs, 0);



sub bc_extract_exp1 {
	my ($fastq,  $type , $ind_len , $bc_len , $pat1 , $pat2 , $BCs, $idx) = @_ ;
	
	#some initializations
	my $total = 0;				# total number of reads
	my $hits = 0;				# reads with successful match
	my $regex = make_regex ($ind_len, $bc_len, $pat1, $pat2, "exp");
	my $len1 = length $pat1;
	my $len2 = length $pat2;
	
	if ($debug) {
		print $error_fh "\nThe regex for $fastq is: $regex\n\n";
	}
	
	die "Cannot open the file  $fastq\n" unless open (FASTQ , "<" , $fastq);
	
	while (my $line = <FASTQ>) {
    
    	die "$fastq is not a proper fastq file\n" if (substr $line, 0, 1) ne "@";
   		$total++;
    	$ line = <FASTQ>;
    	chomp ($line);
    	my $BC_this = 0;
    
   		if ($line =~ $regex) { 	
    	$BC_this =  $1;
    	} elsif ( Hdist ($pat1, substr($line, $ind_len, $len1)) <= ($len1/7)) {
		    
		    if ( Hdist ($pat2, substr($line, ($ind_len + $len1 + $bc_len) , $len2)) < ($len2/7) ) {
		    	
		    		$BC_this = substr ($line , ($ind_len + $len1)  , $bc_len)
		    
		    } elsif ( Hdist ($pat2, substr($line, ($ind_len + $len1 + $bc_len -1) , $len2)) < ($len2/7) ) {
		    	
		    		$BC_this = substr ($line , ($ind_len + $len1)  , $bc_len - 1)
		    
		    } elsif ( Hdist ($pat2, substr($line, ($ind_len + $len1 + $bc_len + 1), $len2)) < ($len2/7)) {
		    	
		    		$BC_this = substr ($line , ($ind_len + $len1)  , $bc_len + 1)
		    }
		    	
   		}
    	if ($BC_this) { 
    		$hits++;
    		#$BCs -> {$BC_this} -> [$type eq "norm" ? 0 : 1] ++;
    		if ($type eq "norm")  {
				$BCs -> {$BC_this} -> [0] ++;
			} else {
				if (exists $BCs -> {$BC_this})  {  ($BCs -> {$BC_this} -> [$idx] ++) }
				}					
		}
    	<FASTQ>;
    	<FASTQ>;
	}
	close FASTQ;
	return ($total , $hits);
}


sub make_regex {
	my ($ind_len, $bc_len, $pat1, $pat2, $kind) = @_;
	
	my $len1 = length $pat1;	#get the length of pat1
	my $short_pat1 = $pat1;
	if  ( $len1  > 10 ) {  
		$short_pat1 = substr($pat1, ($len1 -10)) ;
	} 
	
	my $len2 = length $pat2;	#get the length of pat2
	my $short_pat2 = $pat2;
	if ( $len2  > 10 ) {  
		$short_pat2 = substr($pat2, 0, 10) ;
	} 
	
	# initialization of some important numbers for generating the final regexp
	my $regex;
	my $min1 = $ind_len + $len1 - (length $short_pat1) - 2;
	my $max1 = $ind_len + $len1 - (length $short_pat1) + 2;
	my $min2 = $bc_len - 1;
	my $max2 = $bc_len + 1;
	
	# dealing with special conditions where no index is used ± no first constant part
	if ($ind_len == 0 ) { 
		if ($short_pat1) {   
			$min1 = $len1 - length ($short_pat1);
			if ($min1 >= 2 ) { $min1 -= 2 }
			$max1 = $len1 - (length $short_pat1) + 2;
		} else { $min1 = $max1 = 0 }
	} 
	
	# now depending on the $kind it makes two different types of regex's
	if ($kind eq "exp" ) {
		$regex = "^.{" . $min1 . "," . $max1 . "}" . $short_pat1 . 
				 "(.{" . $min2 . "," . $max2 . "})" . $short_pat2 . ".+";
	} elsif ($kind eq "map" ) {
			$regex = "^.{" . $min1 . "," . $max1 . "}" . $short_pat1 . 
					 "(.{" . $min2 . "," . $max2 . "})" . $short_pat2 . "(.+)";	
		}
	return qr/$regex/;   # this qr makes a regex object which can be directly used
}

sub Hdist {
 	length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] )
}