#!usr/bin/perl;

use strict;
use warnings;
use Data::Dumper;
use List::Util qw(sum);
use Text::LevenshteinXS qw(distance);
use Getopt::Long;


############################	   INITIALIZATIONS 		##################################
# files
my @norm_files;        # files with normalization reads (or reads that you want to use for normalization)
my @exp_files;		   # files with expression counts
my $bc_file;		   # file containing the set of reference barcodes 
my $config_file; 		# configuration file
my $out_dir; 			# output dir

# logical arguments
my $verbose ='';
my $debug = '';

# arguments from the config file
my $ind_len;
my $bc_len;
my $exp_pat1;
my $exp_pat2;
my $hd;



parse_command_line (); # parse the command line

#print "@exp_files\n";
#print "@norm_files\n";


print "Parsing the configuration file ....\n\n" if ($verbose);
config_parse ($config_file); # parse the configuration file


my $debug_file = $out_dir . "debug_file.txt";
my $error_fh;
open $error_fh, ">", $debug_file or die "Cannot open debug_file.txt\n";

if ($debug) {
	print $error_fh my $time = localtime();
	print $error_fh "\nYour script has started working!!!!!!\n\n";
	# details of command line input
	print $error_fh "Here are the details of command line input:\n";
	if (@norm_files) { 
		print $error_fh "normFile:";
		print $error_fh "\t". $_ for @norm_files; print $error_fh "\n";
	}
	print $error_fh "expFile:";
	print $error_fh "\t". $_ for @exp_files; print $error_fh "\n";
	print $error_fh "config:\t\t" . $config_file . "\n";
	print $error_fh "out_dir:\t" . $out_dir . "\n";
	if ($verbose) {
		print $error_fh "The verbose option is on\n";
	}
	print $error_fh "\n";
	
	# details of config file
	print $error_fh "Here is what I got from your config file:\n";
	print $error_fh "index_length:\t" . $ind_len . "\n";
	print $error_fh "barcode_length:\t" . $bc_len . "\n";
	print $error_fh "pat1:\t" . $exp_pat1 . "\n";
	print $error_fh "pat2:\t" . $exp_pat2 . "\n";
	print $error_fh "hd:\t" . $hd . "\n\n";
}

######################## Reading the barcode file into a hash ############################
my %BCs;
print "Reading the barcode file $bc_file ....\n\n" if ($verbose);

%BCs = get_barcodes ($bc_file);

if ($debug) {
	print $error_fh "After initial parsing of $bc_file, the hash BCs looks like this:\n";
	my %h = map { $_ => $BCs {$_} } (sort keys %BCs)[0..5];
	print $error_fh Dumper \%h;
}


# filling up this hash with 0 values
foreach my $bc (keys %BCs) {
	foreach my $sample (0 .. ($#norm_files + $#exp_files + 1)) {
		$BCs {$bc}[$sample + 1] = 0;
	}
}

if ($debug) {
	print $error_fh "After filling the dummy values for counts, the hash BCs looks like this:\n";
	my %h = map { $_ => $BCs {$_} } (sort keys %BCs)[0..5];
	print $error_fh Dumper \%h;
}

############################## Reading normalization files ###############################

my @stats;
if (@norm_files) {
	foreach my $idx (0 .. $#norm_files) {
		
		print "Extracting barcodes from normalization reads file $norm_files[$idx] ....\n\n" if ($verbose);
		
		my @stats_this = bc_extract_exp2 ($norm_files [$idx],  ($idx + 1) , 
						  $ind_len , $bc_len , $exp_pat1 , $exp_pat2 , \%BCs);
		unshift (@stats_this , "norm_" . ($idx+1) );
		push @stats , [ @stats_this ];
	}
	if ($debug) {
		print $error_fh "After filling the normalization counts, the hash BCs looks like this:\n";
		my %h = map { $_ => $BCs {$_} } (sort keys %BCs)[0..5];
		print $error_fh Dumper \%h;
	}
}


############################## Reading expression files ##################################

foreach my $idx (0 .. $#exp_files) {
	
	print "Extracting barcodes from expression reads file $exp_files[$idx] ....\n\n" if ($verbose);
	
	my @stats_this = bc_extract_exp2 ($exp_files [$idx],  ($idx + scalar @norm_files + 1) , 
					  $ind_len , $bc_len , $exp_pat1 , $exp_pat2 , \%BCs);
	unshift (@stats_this , "exp_" . ($idx+1) );			
	push @stats , [ @stats_this ];
}

if ($debug) {
	print $error_fh "After filling the expression counts, the hash BCs looks like this:\n";
	my %h = map { $_ => $BCs {$_} } (sort keys %BCs)[0..5];
	print $error_fh Dumper \%h;
}

##############################   printing count table   ##################################

print "Writing results ....\n\n" if ($verbose);

my $out_put_file = $out_dir . "final_barcode_data_table.txt";

open (MAPOUT , ">" , "$out_put_file") or die "Cannot Open the outputfile $out_put_file: $!\n";
# printing the header
print MAPOUT "id" . "\t" . "barcode". "\t" . $stats[0][0] ;
print MAPOUT "\t" . $stats[$_][0] for (1 .. $#stats);
print MAPOUT "\n";

foreach my $bc (keys %BCs) {
	print MAPOUT $BCs{$bc}[0];
	print MAPOUT "\t$bc";
	print MAPOUT "\t".$_ for @{$BCs {$bc}}[1 .. (scalar @norm_files + scalar @exp_files )] ;
	print MAPOUT "\n";
}

close (MAPOUT);


##############################   printing stats table   ##################################

$out_put_file = $out_dir . "stats_table.txt";

open (STATOUT , ">" , "$out_put_file") or die "Cannot Open the outputfile $out_put_file: $!\n";

# printing the header
print STATOUT "File\tTotal_reads\tReads_with_BC\tExact_match\tPartial_match\n";

foreach my $array (@stats) {
	print STATOUT $_ . "\t" for  @$array [0 .. (scalar @$array - 2)];
	print STATOUT $$array[-1] . "\n";
}

close (STATOUT);


if ($debug) {
	print $error_fh "Your program completed successfully\n";
	print $error_fh my $time = localtime();
}
close $error_fh;  # close the error file handle	

unlink $debug_file unless ($debug);


##########################################################################################
#################************  SUBROUTINE bc_extract_exp2  ****************###############
##########################################################################################
# DESCRIPTION: 
# A subroutine to read the fastq file (typically the Single Read Illumina sequencing runs for
# normalization and expression counts). It reads the file and then extracts the barcode sequence 
# based on the given arguments. The barcodes are stored in a hash as keys. The value for each
# key (barcode) is an anonymous array. 
# If the reads are normalization, then the first element contains the normalization counts
# of that barcode. If the reads are expression than expression counts of that barcode are 
# added as the second element of this anonymous array. 
# 
# INPUT: 
# It takes seven arguments 
# 	a) $fastq -> The name of the fastq file to be read (with full path)
# 	b) $no    -> The number of element in array to put the counts in, when you have multiple files. 
# 	c) $ind_len -> The length of the index sequence. Zero if no index is part of the read. (It
# 				   is important to note that nothing is done with index sequence per se. This 
# 				   is just to get the coordinates of pattern match right.
# 	d) $bc_len -> The length of the barcode.
# 	e) $pat1 -> The sequence of the first constant part.
# 	f) $pat2 -> The sequence of the second constant part.
# 	g) $BCs -> The reference to the hash in which the barcodes and their counts are stored.
# 
# DEPENDS: 
# It depends on 
# 	Hdist (see below)
# 	make_regex (see below) 
# 
# WORKING: 
# It first determines the length of the $pat1 and the $pat2 (constant parts of the read). 
# It generates the regular expression for matching using subroutine make_regex. One example
# of such a regex is.
# 		^.{24,28}ACAACTCGAG(.{15,17})TGATCCTGCA.+
# One advantage of this regexp is that it allows for some (but not too much) flexibility in 
# exact position of the match of the first pattern. 
# Then it reads the fastq file line by line. The first line should start with @ otherwise it will die. 
# Then it reads in the next line which is the sequence read and matches the regexp pattern. 
# If the pattern is matched then the barcode is extracted using $1 variable (it store the recently
# matched pattern).  
# If the regexp pattern is not matched then, it finds the Hamming distance using Hdist subroutine between
# the $pat1 and the corresponding part of the read. If this is below (length $pat1)/7 then same 
# procedure is repeated for $pat2 and the correponding part of the read. Because the variance of
# one nucleotide in the barcode length is allowed (see above) so the $pat2 is matched at three
# possible places in the read. If the Hdist is below (length $pat2)/7 then the barcode sequence
# is extracted. 
# In the end the barcode sequence is fed as keys to the hash BCs. The value of each key (barcode) 
# is an anonymous array with first element is normalization count and the second one (if 
# the $type is "exp") is expression count.
# During the the procedure the subroutine also takes record of total reads and the reads
# which could be matched properly.
# 
# OUTPUT: 
# It returns the total reads and the reads matched.
# It populates the hash BCs with barcode sequences(as keys) and an anonymous array with 
# first element is normalization count and the second one (if the $type is "exp") is expression 
# count. One key of such an hash would look like this 
#	TAAAGTCCTGGCTTAA => [
#						345,
#						201
#						]
##########################################################################################

sub bc_extract_exp2 {
	my ($fastq,  $no , $ind_len , $bc_len , $pat1 , $pat2 , $BCs) = @_ ;
	
	#some initializations
	my $total = 0;				# total number of reads
	my $hits = 0;				# reads with successful match
	my $ex_match =0;			# barcodes with exact match
	my $part_match =0;			# barcodes which match only partially to only one of the reference barcodes
	my $regex = make_regex ($ind_len, $bc_len, $pat1, $pat2, "exp");
	my $len1 = length $pat1;
	my $len2 = length $pat2;
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
    		if (exists $BCs -> { $BC_this }) { 
    			$ex_match++;	
    			$BCs -> { $BC_this } -> [$no] += 1;	
    		} elsif ( $hd != 0 ) {
					my @mutants = grep { distance ($_ , $BC_this) <= $hd } keys %$BCs;
					#print "@mutants\n";
					if (scalar @mutants == 1 ) {
						$part_match++;
						$BC_this = shift @mutants;
						$BCs -> { $BC_this } -> [$no] += 1;
					}
				}
		}
    	<FASTQ>;
    	<FASTQ>;
	}
	close FASTQ;
	return ($total , $hits, $ex_match, $part_match);
}




##########################################################################################
#################***************   SUBROUTINE make_regex  ******************##############
##########################################################################################
# DESCRIPTION: 
# A subroutine to make a regular expression string from given arguments for proper matching
# in cases of different formats of read structure
# 
# INPUT: 
# It takes five arguments 
# 	a) $ind_len  -> The length of the index sequence (0 if no index is used)
# 	b) $bc_len   -> The length of the barcode
# 	c) $pat1	 -> The sequence of first constant part ('NA' if no constant part 1 is used)
# 	d) $pat1	 -> The sequence of second constant part (must be DNA string)
# 	e) $kind	 -> Either "exp" (for normalization/expression reads where only barcode is to be extracted)
# 					or "map" (for mapping reads where both the barcode and the neighboring genomic DNA
# 					is to be extracted) 
# 						
# DEPENDS: 
# no dependencies
# 
# WORKING: 
# From each of $pat1 and $pat2, it finds if they are longer than 10 bases. If so, it takes 
# only the last 10 bases of $pat2 and the first 10 bases of $pat2 are saved into two additional
# variables $short_pat1 and $short_pat2. 
# Next, it determines min and max values in order to put into the regex, as shown below
# ^{ $min1, $max1} $short_pat1 (.{(barcode length -1) , (barcode length + 1)} $short_pat2 .+
# Then it tries to make sure that the $min1 is not a negative number (as it can happen
# if no index is used). Also if it is above 2 (in cases of index absence where $pat1 is >1 longer
# than $short_pat1) it subtracts 2 from min1 to allow for some positional flexibility in 
# patten match.
# In the end it makes the regex from all the different components generated.
# If it is for normalization/expression reads then the last .+ does not get parenthesis. e.g.,
# ^.{24,28}ACAACTCGAG(.{15,17})TGATCCTGCA.+
# If it is for mapping reads then the last .+ gets parthesis around it. e.g.,
# ^.{24,28}ACAACTCGAG(.{15,17})TGATCCTGCA(.+)
# 
# OUTPUT: 
# The regular expression string
##########################################################################################

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
	return qr/$regex/;		# this qr makes a regex object which can be directly used
}



##########################################################################################
#################***************   SUBROUTINE Hdist  ******************###################
##########################################################################################
# DESCRIPTION: 
# A subroutine to find Hamming distance between two stings.
# 
# INPUT: 
# It takes two arguments 
# 	a) $_[0]  -> The first string.
# 	b) $_[1]  -> The second string.
# 
# DEPENDS: 
# no dependencies
# 
# WORKING: 
# It is not entirely clear how it works. It is taken from 'perlmonks.com'.
# Although it is claimed to work on strings of equal lengths, but it does not complain if the 
# lengths of the strings differ. However the Hamming distance is determined based on the length
# of first--that is HD = length $_[0] - difference between the two strings.
# If you have an addition/deletion on the end of second string and the rest of the two strings 
# are same then it will find the distance as the length of addition, but if the addition or 
# deletion is at the start of either string then hamming distance can be very large because 
# everything is shifted
# Not good for tranversions at all.
# OUTPUT: 
# The Hamming distance between the two strings. 
##########################################################################################

sub Hdist {
 	length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] )
}










sub get_barcodes {
	my ($filename) = @_;
	my %BCs;
	open ( BCFILE, "<" , $filename ) or die "Error: failed to open barcode file ($filename)\n";
	
	while (my $line = <BCFILE>) {
		next if ($line =~ m/^#/);
		chomp ($line);
		my ($ident, $barcode) = split (/\s+/ , $line);

		$barcode = uc($barcode);

		# Sanity checks on the barcodes
		die "Error: bad data at barcode file ($filename) line $.\n" unless defined $barcode;
		die "Error: bad barcode value ($barcode) at barcode file ($filename) line $.\n"
			unless $barcode =~ m/^[AGCTN]+$/;

		die "Error: bad identifier value ($ident) at barcode file ($filename) line $. (must be alphanumeric)\n" 
			unless $ident =~ m/^\w+$/;

		die "Error: badcode($ident => $barcode) is shorter or equal to maximum number of " .
		    "mismatches ($hd). This makes no sense. Specify fewer  mismatches.\n" 
		    	if length ($barcode) <= $hd;
		    	
	 	$BCs {$barcode } -> [0] = $ident;

	}
	close BCFILE;
	return %BCs;
}

























##||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||##
##**************************************************************************************##
######				  		Sub-routines for initializations  						######
##**************************************************************************************##
##||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||##


##########################################################################################
##############************   SUBROUTINE parse_command_line  ****************##############
##########################################################################################
# DESCRIPTION: 
# It parses the command line and does some sanity checks and prints out error messages in case
# something is wrong. 
# 
# INPUT: 
# It takes no arguments 
# 
# DEPENDS: 
# It depends on Getopt::Long
# 
# WORKING: 
# Parses the command line and puts the values into the variables declared before this 
# subroutine is called.
# 
# OUTPUT: 
# Changes the values of the declared variables according to the command line
##########################################################################################
sub parse_command_line {
	my $help;

	usage() if (scalar @ARGV==0);

	my $result = GetOptions ( "normFile|nf|n=s" => \@norm_files,
				  "expFile|ef|e=s"  => \@exp_files,
				  "bcFile|bc|b=s"  => \$bc_file,
				  "out_dir|od|o=s" => \$out_dir,
				  "config|c=s" => \$config_file,
				  "help|h" => \$help,
				  "verbose|v" => \$verbose,
				  "debug|d" => \$debug,
				  ) ;
	
	usage() if ($help);
	
	if (@norm_files) {
	
		@norm_files = split(/,/ , join(',', @norm_files));
		foreach my $file (@norm_files) {
			die "Error: cannot open $file\n" unless (-r $file);
		} 
	} else { print "\nNo normalization reads file/s specified .. \n\n"; }
	
	
	die "Error: expression reads file not specified (use '--expFile [FILENAME]')\n" unless @exp_files;
	
	@exp_files = split(/,/ , join(',', @exp_files));
	foreach my $file (@exp_files) {
		die "Error: cannot open $file\n" unless (-r $file);
	}
	
	if (scalar @norm_files == 1 && scalar @exp_files > 1) {
		print "\nWarning! there is only one normalization file and more than one expression files. 
All expression files will be normalized against this file\n\n";
	} elsif (scalar @norm_files != scalar @exp_files) {
			die "Error: Unequal number of expression and normalization files
Either specify only equal number of expression and normalization files or only one normalization file\n";
		}
	
	die "Error: the barcode file not specified (use '--bcFile [FILENAME]')\n" unless defined $bc_file;
	die "Error: cannot open $bc_file\n" unless (-r $bc_file);
	
	die "Error: configuration file not specified (use '--config [FILENAME]')\n" unless defined $config_file;
	die "Error: cannot open $config_file\n" unless (-r $config_file);
	
	$out_dir =~ s/\/*$/\//;
	die "Error: output dir not specified (use '--out_dir [DIRECTORYNAME]')\n" unless defined $out_dir;
	die "Error: the output dir $out_dir does not exist\n" unless  (-d $out_dir);

	#die "Error: the TRIP format not defined (use '--format [FORMAT]\n)" unless defined ;
	exit unless $result;
}


##########################################################################################
##############************   	 SUBROUTINE config_parse    ****************##############
##########################################################################################
# DESCRIPTION: 
# It parses the command line and does some sanity checks and prints out error messages in case
# something is wrong. 
# 
# INPUT: 
# It takes one argument
# 	a) $config -> The name of the configuration file
# 
# DEPENDS: 
# non dependencies
# 
# WORKING: 
# Reads in the file and removes the new lines. Removes comments and any white spaces.
# checks if something is left. If yes, it splits based on "=" to separate the argument name
# and argument value into a key value pair put into a hash.
# For this hash the individual arguments (the names of which are declared before calling 
# this subroutine) are extracted and sanity checked.
# OUTPUT: 
# Changes the values of the declared variables according to the configuration file
##########################################################################################

sub config_parse {
	my ($config) = @_;
	my %User_Preferences;

	open (CONFIG, "<", $config) or die "Cannot open configuration file: $config\n";

	while (<CONFIG>) {
		chomp;                  # no newline
		s/#.*//;                # no comments
		s/^\s+//;               # no leading white
		s/\s+$//;               # no trailing white
		next unless length;     # anything left?
		my ($var, $value) = split(/\s*=\s*/, $_, 2);
		$User_Preferences{$var} = $value;
	} 
	close CONFIG;
	
	# some sanity checks
	$ind_len = $User_Preferences { 'index_length' };
	die "Could not find any specification of index_length\n" unless defined $ind_len;
	die "The index_length given is not an integer\n" unless $ind_len =~ /^\d+$/;
	
	$bc_len = $User_Preferences { 'barcode_length' };
	die "Could not find any specification of barcode_length\n" unless defined $bc_len;
	die "The barcode_length given is not an integer\n" unless $bc_len =~ /^\d+$/;
	
	$exp_pat1 = $User_Preferences { 'pat1' };
	die "Could not find any specification of pat1\n" unless defined $exp_pat1;
	$exp_pat1 = uc $exp_pat1;
	if ($exp_pat1 =~ /^NA$/) { $exp_pat1 = ''} else {
		die "The pat1 contains non-DNA characters\n" if $exp_pat1 =~ /[^ACTG]/;
	}
	$exp_pat2 = $User_Preferences { 'pat2' };
	die "Could not find any specification of pat2\n" unless defined $exp_pat2;
	$exp_pat2 = uc $exp_pat2;
	die "The pat2 contains non-DNA characters\n" if $exp_pat2 =~ /[^ACTG]/;
	
	$hd = $User_Preferences { 'hd' };
	die "Could not find any specification of hd\n" unless defined $hd;
	die "The hd given is not an integer\n" unless $hd =~ /^\d+$/;
	
}




##########################################################################################
##############************   		 SUBROUTINE usage   	 ****************#############
##########################################################################################
# DESCRIPTION: 
# It prints out help info.
# 
# INPUT: 
# It takes no argument
# 
# DEPENDS: 
# non dependencies
# 
# WORKING: 
# prints out the statement between <<EOF and EOF and exits
# 
# OUTPUT: 
# Prints the help info on screen
##########################################################################################

sub usage()
{
print<<EOF;

An accessory Util for TRIP Aanalysis Software Kit (TASK), by Waseem Akhtar, Johann de Jong, Jelle  and Ludo Pagie  (w.akhtar\@nki.nl), 11sep2013
This program reads FASTQ files from any type of barcode dependent parallel reporter assay experiment and 
extracts the barcode sequence for each read based on read patterns described in the config file. These barcodes 
are then matched with the reference set of barcodes and their counts are determined.


usage: $0 --expFile [EXPRESSION_FILE] --normFile [NORMALIZATION_FILE] --bcFile [BARCODE_FILE] 
		  --config [CONFIGURATION_FILE] --out_dir [OUTPUT_DIRECTORY] [--help]

Arguments:

-e/-ef/--expFile	EXPRESSION_FILE		- file/s (multiple files should be comma separated) containing expression reads. Fastq (unzipped)

-n/-nf/--normFile	NORMALIZATION_FILE	- file/s (multiple files should be comma separated) containing normalization reads. Fastq (unzipped)

-b/-bcFile		BARCODE_FILE		- file containing barcode sequences

-c/-config		CONFIGURATION_FILE	- file containing extra arguments (for an example see below)

-o/-od/--out_dir	OUTPUT_DIRECTORY	- path to the directory where output should be saved

-h/--help					- Prints out usage info	
-v/--verbose					- Prints out running commentary during the execution of the script

Example of a Configuration File:

index_length	= 10  # length of the index (0 = if no index used) 
barcode_length	= 16  # length of the barcode
pat1		= GTCACAAGGGCCGGCCACAACTCGAG  # first constant part of normalization and expression reads
pat2		= TGATCCTGCAGTGTCACCTAAATCGTATGCGGCCGCGAATTCTTACTT  # second constant part of normalization and expression reads
hd		= 2  # the maximum distance from a given barcode for a query barcode to be included in the analysis
		
EOF

exit 1;
}


