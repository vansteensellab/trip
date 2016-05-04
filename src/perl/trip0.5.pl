#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use List::Util qw(sum);
use Text::LevenshteinXS qw(distance);
use Getopt::Long;
use File::Spec;
use File::Basename;

############################	   INITIALIZATIONS 		##################################
# expression files
my @exp_files;
my $norm_file;

# mapping files
my $map_file_for;
my $map_file_rev;

my $config_file; 		# configuration file
my $out_dir; 			# output dir
my $map_style = 'n';		# the mapping format, default means no mapping 

# arguments from the config file
my $ind_len;
my $bc_len;
my $exp_pat1;
my $exp_pat2;
my $hd;
my $map_pat1;
my $map_pat2;
my $map_pat_rev;
my $cores;
my $bowtie_base;
my $max_dist_for;
my $max_dist_rev;
my $min_counts;
my $verbose = '';
my $debug = '';

############################	 Getting vairables 		##################################

parse_command_line (); # parse the command line

print "\nParsing the configuration file ....\n\n" if ($verbose);
config_parse ($config_file); # parse the configuration file

# intializations for degug file
my $debug_file = $out_dir . "debug_file.txt";
my $error_fh;
open $error_fh, ">", $debug_file or die "Cannot open debug_file.txt\n";

if ($debug) {
	print $error_fh my $time = localtime();
	print $error_fh "\nYour script has started working!!!!!!\n\n";
	# details of command line input
	print $error_fh "Here are the details of command line input:\n";
	print $error_fh "normFile:\t" . $norm_file . "\n";
	print $error_fh "expFile:";
	print $error_fh "\t". $_ for @exp_files; print $error_fh "\n";
	print $error_fh "config:\t\t" . $config_file . "\n";
	print $error_fh "out_dir:\t" . $out_dir . "\n";
	print $error_fh "map:\t\t" . $map_style . "\n";
	if ($map_style ne 'n') {
		print $error_fh "mapFor:\t" . $map_file_for . "\n";
		if($map_style ne 'f') {
			print $error_fh "mapRev:\t" . $map_file_rev . "\n";
		}
	}
	if ($verbose) {
		print $error_fh "The verbose option is on\n";
	}
	print $error_fh "\n\n";
	# details of config file
	print $error_fh "Here is what I got from your config file:\n";
	print $error_fh "index_length:\t" . $ind_len . "\n";
	print $error_fh "barcode_length:\t" . $bc_len . "\n";
	print $error_fh "pat1:\t" . $exp_pat1 . "\n";
	print $error_fh "pat2:\t" . $exp_pat2 . "\n";
	print $error_fh "hd:\t" . $hd . "\n";
	if ($map_style ne 'n') {
		print $error_fh "map_pat1:\t" . $map_pat1 . "\n";
		print $error_fh "map_pat2:\t" . $map_pat2 . "\n";
		if ($map_style ne 'f') {
			print $error_fh "map_pat_rev:\t" . $map_pat_rev . "\n";
		}
		print $error_fh "cores:\t" . $cores . "\n";
		print $error_fh "bowtie_base:\t" . $bowtie_base . "\n";
		print $error_fh "max_dist_for:\t" . $max_dist_for . "\n";
		if ($map_style ne 'f') {
			print $error_fh "max_dist_rev:\t" . $max_dist_rev . "\n";
		}
	}
	print $error_fh "min_counts:\t" . $min_counts . "\n";
}

################# Handling the Normalization count file ##################################

my %BCs;

print "Extracting barcodes from the normalization file $norm_file ....\n\n" if ($verbose);

my @stats_norm = bc_extract_exp1 ($norm_file, "norm" , $ind_len , $bc_len , $exp_pat1 , $exp_pat2 , \%BCs, 0);

# if (exists $BCs {'AATAACCTGCCTTTCT'}) { print "Yes AATAACCTGCCTTTCT is $BCs{'AATAACCTGCCTTTCT'}->[0] \n"; }
# if (exists $BCs {'CAATATGCATTGTATC'}) { print "Yes CAATATGCATTGTATC is $BCs{'CAATATGCATTGTATC'}->[0] \n"; }

my $all_unique_bc = scalar (keys %BCs);
print "Total number of unique barcodes in $norm_file:\t" . $all_unique_bc . "\n\n" if ($verbose);

%BCs =  map { $_ => $BCs {$_} } grep { $BCs {$_} -> [0] >= $min_counts } keys %BCs; # only those BCs with count above 5

die "Error: there are no barcodes having reads more than your specified min_counts: $min_counts\nPlease lower this parameter\n"
	if (! keys %BCs);

my $bc_above_min_count = scalar (keys %BCs);
print "The barcodes above min_count in $norm_file:\t" . $bc_above_min_count . "\n\n" if ($verbose);

if ($debug) {
	print $error_fh "After initial parsing of $norm_file, the hash BCs looks like this:\n";
	my %h = map { $_ => $BCs {$_} } (sort keys %BCs)[0..5];
	print $error_fh Dumper \%h;
	print $error_fh "Total number of unique barcodes in $norm_file:\t" . $all_unique_bc . "\n\n";
	print $error_fh "The barcodes above min_count in $norm_file:\t" . $bc_above_min_count . "\n\n";
}

%BCs = bc_group (\%BCs, $hd);   # groups barcode based on hamming distance 2 or less
my $gen_bc = scalar (keys %BCs);
print "The total number of genuine barcodes:\t" . $gen_bc . "\n\n" if ($verbose);

if ($debug) {
	print $error_fh "After grouping the barcodes, the hash BCs looks like this:\n";
	my %h = map { $_ => $BCs {$_} } (sort keys %BCs)[0..5];
	print $error_fh Dumper \%h;
	print $error_fh "The total number of genuine barcodes:\t" . $gen_bc . "\n\n";
}

push (@stats_norm, sum map { $BCs{$_} -> [0] } keys %BCs) ;
unshift (@stats_norm, "Norm");

foreach my $idx (0 .. $#exp_files) {
	foreach my $key (keys %BCs) { $BCs {$key} -> [$idx + 1] = 0 } # add a dummy value 0 for the expression counts
}

if ($debug) {
	print $error_fh "After filling the dummy values for expression counts, the hash BCs looks like this:\n";
	my %h = map { $_ => $BCs {$_} } (sort keys %BCs)[0..5];
	print $error_fh Dumper \%h;
}

################# Handling the Expression count files ####################################

my %stats_exp;
foreach my $idx (0 .. $#exp_files) {
	
	print "Extracting barcodes from the expression file $exp_files[$idx] ....\n\n" if ($verbose);	
	my @stat_this = bc_extract_exp1 ($exp_files [$idx], "exp" , $ind_len , $bc_len , $exp_pat1 , $exp_pat2 , \%BCs, ($idx+1));	
	push (@stat_this, sum map { $BCs{$_} -> [$idx + 1] } keys %BCs) ;
	unshift (@stat_this, "Exp_". ($idx+1));
	$stats_exp {$idx} = [@stat_this];
	
	if ($debug) {
		print $error_fh "After analyzing expression file $exp_files[$idx]\nthe hash BCs looks like this:\n";
		my %h = map { $_ => $BCs {$_} } (sort keys %BCs)[0..5];
		print $error_fh Dumper \%h;
	}
}

################# Handling the Mapping reads file ########################################

# initialization of variables

my %IDs;
my @stats_map;
my %map_hash_for;
my %map_hash_rev;

					######### Forward Mapping ###############
if ($map_style ne 'n') {
	print "Extracting barcodes and genomic regions from the mapping file $map_file_for ....\n\n" if ($verbose);
	my $for_map_tbl =  $out_dir . "mapping_table_for.txt";	
	
	@stats_map = mapping_for1 ($map_file_for, $ind_len, $bc_len , $map_pat1, 
					$map_pat2, $for_map_tbl, \%IDs, \%BCs);
	
	my $sam_for = $out_dir . "samFor.sam";
	my $align_command = "bowtie2 -p " . $cores . " -t -f --very-sensitive -x " . 
						$bowtie_base . " -U " . $for_map_tbl ." -S " . $sam_for;
	print $align_command;
	print "TEST";				
	system ($align_command);
	
	if ($debug) {
		print $error_fh $align_command;
		print $error_fh "After analyzing the mapping file $map_file_for, the hash IDs looks like this:\n";
		my %h = map { $_ => $IDs {$_} } (sort keys %IDs)[0..5];
		print $error_fh Dumper \%h;
		print $error_fh "\nHere is the align_command:\t$align_command\n\n";
	}
	
	parse_sam ($sam_for , \%map_hash_for); 		# parse the sam file
	
	if ($debug) {
		print $error_fh "After parsing the sam file for forward reads, the hash map_hash_for looks like this:\n";
		my %h = map { $_ => $map_hash_for {$_} } (sort keys %map_hash_for)[0..5];
		print $error_fh Dumper \%h;
	}
	
	top_map (\%map_hash_for, $max_dist_for);		# refine and get the mapping info
	
	if ($debug) {
		print $error_fh "After applying sub top_map, the hash map_hash_for looks like this:\n";
		my %h = map { $_ => $map_hash_for {$_} } (sort keys %map_hash_for)[0..5];
		print $error_fh Dumper \%h;
		print $error_fh $for_map_tbl
	}
	
	# unlink $sam_for;							# delete the sam file

	# unlink $for_map_tbl;						# delete the table
	unshift (@stats_map, "Mapping");
	
					############ Reverse Mapping ###############
	if ($map_style eq 'b' | $map_style eq 'r') {
		
		my %re_map;  # a hash for collecting those cases which need re-alignment that is (>17)S\d{2}M
		
		print "Extracting genomic regions from the mapping file $map_file_rev ....\n\n" if ($verbose);
		my $rev_map_tbl =  $out_dir . "mapping_table_rev.txt";			
		fastq_subset ($map_file_rev, $map_pat_rev, $rev_map_tbl, \%IDs);
		
		my $sam_rev = $out_dir . "samRev.sam";
		$align_command = "bowtie2 -p " . $cores . " -t -f --very-sensitive-local -x " . 
						 $bowtie_base . " -U " . $rev_map_tbl ." -S " . $sam_rev;
						 
		system ($align_command);
		
		parse_sam ($sam_rev , \%map_hash_rev, \%re_map);  # parse the sam file
		
		if ($debug) {
			print $error_fh "\nHere is the align_command for mapping file $map_file_rev:\t$align_command\n\n";
			print $error_fh "After parsing the sam file for reverse reads, the hash map_hash_rev looks like this:\n";
			my %h = map { $_ => $map_hash_rev {$_} } (sort keys %map_hash_rev)[0..5];
			print $error_fh Dumper \%h;
		}
		
		if (keys %re_map) {
			my $sam_rev2 = $out_dir . "samRev2.sam";
			my $rev_map_tbl2 = $out_dir . "mapping_table_rev2.txt";
			
			print "Writing the fasta file for re-alignment\n" if ($verbose);
			open (REMAP, ">", $rev_map_tbl2) or die "Cannot open $rev_map_tbl2 for writing\n";
			foreach my $key (keys %re_map) {
				print REMAP ">". $key . "\n" . $re_map { $key } . "\n";
			}
			close (REMAP);
			
			if ($debug) {
				print $error_fh "\nHere is how re_map hash looks\n\n";
				my %h = map { $_ => $re_map {$_} } (sort keys %re_map)[0..5];
				print $error_fh Dumper \%h;
			}
			
			print "Starting re-alignment to the genome\n" if ($verbose);
			my $align_command2 = "bowtie2 -p " . $cores . " -t -f --very-sensitive-local -x " . 
						 		 $bowtie_base . " -U " . $rev_map_tbl2 ." -S " . $sam_rev2;
						 		
			system ($align_command2);
			
			parse_sam ($sam_rev2 , \%map_hash_rev, \%re_map);  # parse the sam file
			
			if ($debug) {
				print $error_fh "\nHere is the align_command for re-mapping of file $map_file_rev:\t$align_command2\n\n";
				print $error_fh "After parsing the sam file re-mapped reads, the hash map_hash_rev looks like this:\n";
				my %h = map { $_ => $map_hash_rev {$_} } (sort keys %map_hash_rev)[0..5];
				print $error_fh Dumper \%h;
			}
			# unlink $sam_rev2;							# delete the sam file
			# unlink $rev_map_tbl2;						# delete the table
			
		}
		
		
		top_map (\%map_hash_rev, $max_dist_rev);    # refine and get the mapping info
		
		if ($debug) {
			print $error_fh "After applying sub top_map, the hash map_hash_rev looks like this:\n";
			my %h = map { $_ => $map_hash_rev {$_} } (sort keys %map_hash_rev)[0..5];
			print $error_fh Dumper \%h;
		}
		
		# unlink $sam_rev;							# delete the sam file
		# unlink $rev_map_tbl;						# delete the table
	}
}
########################   Writing out the results #######################################



					##### Writing the TRIP data table  #######
print "Writing results ....\n" if ($verbose);
my $out_put_file = $out_dir . "final_TRIP_data_table.txt";

open (MAPOUT , ">" , "$out_put_file") or die "Cannot Open the outputfile $out_put_file: $!\n";

my @map_head = qw /chr ori pos reads mapq freq1 freq2/;
my $map_head_for = join ("\t" , map { $_ . "_f" } @map_head ); 
my $map_head_rev = join ("\t" , map { $_ . "_r" } @map_head ); 
my $exp_head = join ("\t" , map { "exp_" . $_ } 1.. scalar @exp_files);
my $unmapped_string = "NA\t" x 7;

# Printing the header 
print MAPOUT "barcode\t" . "norm\t" . $exp_head;

if ( $map_style eq 'f' ) {
	print MAPOUT "\t" . $map_head_for;
}
if ( $map_style eq 'r' ) {
	print MAPOUT "\t" . $map_head_rev;
}
if ( $map_style eq 'b' ) {
	print MAPOUT "\t" . $map_head_for . "\t" . $map_head_rev;
}

print MAPOUT "\n";

# Printing the data
foreach my $bc (sort { $BCs{$b} ->[0] <=> $BCs{$a} ->[0] } keys %BCs) {
	print MAPOUT $bc;
	print MAPOUT "\t" . $_ for  @{ $BCs {$bc} };

	if ($map_style eq 'b' | $map_style eq 'f') {
		if (exists $map_hash_for {$bc} ) {
			print MAPOUT "\t" . $_ for @{ $map_hash_for {$bc} };
		} else { print MAPOUT "\t" . $unmapped_string }
	}
	if ($map_style eq 'b' | $map_style eq 'r') {
		if (exists $map_hash_rev {$bc} ) {
			print MAPOUT "\t" . $_ for @{ $map_hash_rev {$bc} };
		} else { print MAPOUT "\t" . $unmapped_string }
	}
	print MAPOUT "\n";
}
close MAPOUT;


	
	
						##### Writing the stats  #######

$out_put_file = $out_dir . "stats.txt";

open (STATOUT , ">" , "$out_put_file") or die "Cannot Open the outputfile $out_put_file: $!\n";
# print in the header
print STATOUT "Type\t" . "Total\t" . "Matching\t" . "genBC_containing\n"; 
print STATOUT $_."\t" for @stats_norm;
print STATOUT "\n";

foreach my $idx (sort keys %stats_exp) {
	print STATOUT $_."\t" for @{ $stats_exp {$idx} };
	print STATOUT "\n";
}
if ($map_style) {
	print STATOUT $_."\t" for @stats_map;
	print STATOUT "\n";
}

close STATOUT;




######################### plotting using R script ########################################

my $path = dirname(File::Spec->rel2abs(__FILE__));

system ("Rscript $path/../R/trip_plot.R $out_dir");

if ($debug) {
	print $error_fh "Your program completed successfully\n";
	print $error_fh my $time = localtime();
}
close $error_fh;  # close the error file handle	

unlink $debug_file unless ($debug);

############################### End of the Script ########################################

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

	my $result = GetOptions ( "normFile|nf|n=s" => \$norm_file,
				  "expFile|ef|e=s"  => \@exp_files,
				  "mapFor|mf|f=s"  => \$map_file_for,
				  "mapRev|mr|r=s" => \$map_file_rev,
				  "out_dir|od|o=s" => \$out_dir,
				  "map|m=s" => \$map_style,
				  "config|c=s" => \$config_file,
				  "help|h" => \$help,
				  "verbose|v" => \$verbose,
				  "debug|d" => \$debug,
				  ) ;
	
	usage() if ($help);
	
	die "Error: normalization reads file not specified (use '--normFile [FILENAME]')\n" unless defined $norm_file;
	die "Error: cannot open $norm_file\n" unless (-r $norm_file);
	
	
	if (@exp_files) {
		@exp_files = split(/,/ , join(',', @exp_files));
		foreach my $file (@exp_files) {
			die "Error: cannot open $file\n" unless (-r $file);
		}
	} else { 
			die "Error: expression read file/s not specified (use '--expFile [FILE1, FILE2,...]')\n" 
		}
	
	die "Error: configuration file not specified (use '--config [FILENAME]')\n" unless defined $config_file;
	die "Error: cannot open $config_file\n" unless (-r $config_file);
	
	die "Error: output dir not specified (use '--out_dir [DIRECTORYNAME]')\n" unless defined $out_dir;
	die "Error: the output dir $out_dir does not exist\n" unless  (-d $out_dir);

	$out_dir =~ s/\/*$/\//;
	
	if ($map_style ne 'n' && $map_style ne 'f' && $map_style ne 'r' && $map_style ne 'b') {
		die "Error: invalid --map input	'$map_style'\nIt can be only:		n OR f OR r OR b\n";
	}
	
	if ($map_style eq 'b' | $map_style eq 'r') {
	
		die "Error: files for mapping reads not specified (use '--mapFor [FILENAME] --mapRev [FILENAME]')\n" 
			unless (defined $map_file_for and defined $map_file_rev);
		die "Error: cannot open $map_file_for\n" unless (-r $map_file_for);
		die "Error: cannot open $map_file_rev\n" unless (-r $map_file_rev);
		
	}
	
	if ($map_style eq 'f') {
		die "Error: mapping forward reads file not specified (use '--mapFor [FILENAME]')\n" unless defined $map_file_for;
		die "Error: cannot open $map_file_for\n" unless (-r $map_file_for);
	}

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
	
	$min_counts = $User_Preferences { 'min_counts' };
	die "Could not find any specification of min_counts\n" unless defined $min_counts;
	die "The min_counts given is not an integer\n" unless $min_counts =~ /^\d+$/;
	
	if ($map_style ne 'n') {
		
		$map_pat1 = $User_Preferences { 'map_pat1' };
		die "Could not find any specification of map_pat1\n" unless defined $map_pat1;
		$map_pat1 = uc $map_pat1;
		if ($map_pat1 =~ /^NA$/) { $map_pat1 = ''} else {
			die "The map_pat1 contains non-DNA characters\n" if $map_pat1 =~ /[^ACTG]/;
		}
	
		$map_pat2 = $User_Preferences { 'map_pat2' };
		die "Could not find any specification of map_pat2\n" unless defined $map_pat2;
		$map_pat2 = uc $map_pat2;
		die "The map_pat2 contains non-DNA characters\n" if $map_pat2 =~ /[^ACTG]/;
		
		$cores = $User_Preferences { 'cores' };
		die "Could not find any specification of cores\n" unless defined $cores;
		die "The cores given is not an integer\n" unless $cores =~ /^\d+$/;
	
		$bowtie_base = $User_Preferences { 'bowtie_base' };
		$bowtie_base =~ s/\/*$//;
		
		die "The bowtie_base for alignment $bowtie_base does not exist\n" unless (-e "$bowtie_base.1.bt2");
	
		$max_dist_for = $User_Preferences { 'max_dist_for' };
		die "Could not find any specification of max_dist_for\n" unless defined $max_dist_for;
		die "The max_dist_for given is not an integer\n" unless $max_dist_for =~ /^\d+$/;
	
		$max_dist_rev = $User_Preferences { 'max_dist_rev' };
		die "Could not find any specification of max_dist_rev\n" unless defined $max_dist_rev;
		die "The max_dist_rev given is not an integer\n" unless $max_dist_rev =~ /^\d+$/;
	}
	
	if ($map_style eq 'b' | $map_style eq 'r') {
	
		$map_pat_rev = $User_Preferences { 'map_pat_rev' };
		die "Could not find any specification of map_pat_rev\n" unless defined $map_pat_rev;
		$map_pat_rev = uc $map_pat_rev;
		die "The map_pat_rev contains non-DNA characters\n" if $map_pat_rev =~ /[^ACTG]/;
	}
	
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

trip version 0.3.5
A tool from TRIP Aanalysis Software Kit (TASK), by 
Waseem Akhtar, Johann de Jong, Jelle ten Hoeve and Ludo Pagie  (w.akhtar\@nki.nl), 11sep2013
This program reads FASTQ files from TRIP experiments, extracts out the barcodes, filters out mutant barcodes and makes
a table of barcode counts in normalization and expression reads. In addition it extracts the genomic neighborhood 
sequences from the mapping reads and finds the locations (and their frequency) associated with these barcodes.


usage: $0 --normFile [NORMALIZATION_FILE] --expFile [EXPRESSION_FILE] 
		  --mapFor [MAPPING_FORWARD] --mapRev [MAPPING_REVERSE]
		  --config [CONFIGURATION_FILE] --out_dir [OUTPUT_DIRECTORY] 
		  --map [n|f|r|b] [--help] [--vervose] [--debug]

Arguments:

-n/-nf/--normFile	NORMALIZATION FILE	- file containing normalization reads. Fastq (unzipped)

-e/-ef/--expFile	EXPRESSION_FILE		- file containing expression reads. Fastq (unzipped)

-f/-mf/--mapFor		MAPPING_FORWARD		- file containing forward mapping reads. Fastq (unzipped)

-r/-mr/--mapRev		MAPPING_REVERSE		- file containing reverse mapping reads. Fastq (unzipped)

-c/--config		CONFIGURATION_FILE	- file containing extra arguments (for an example see below)

-o/-od/--out_dir	OUTPUT_DIRECTORY	- path to the directory where output should be saved

-m/-map			[n OR b OR r OR f]	- should mapping be done from both forward and reverse 
						reads (b) or from forward reads only (f) or from reverse reads only (r). 
						If 'f' only --mapFor needs to be specified otherwise (b or r) both 
						--mapFor and --mapRev needs to be specified.
						The default value is undefined which means no mapping.
										  
-h/--help					- Prints out usage info	
-v/--verbose					- Prints out running commentary during the execution of the script
-d/--debug					- Prints out the progress into a file debug_file.txt in the output directory

Example of a Configuration File:

index_length	= 10  # length of the index (0 = if no index used) 
barcode_length	= 16  # length of the barcode
pat1		= GTCACAAGGGCCGGCCACAACTCGAG  # first constant part of normalization and expression reads
pat2		= TGATCCTGCAGTGTCACCTAAATCGTATGCGGCCGCGAATTCTTACTT  # second constant part of normalization and expression reads
hd		= 2  # hamming distance threshould to filter out barcode mutants
map_pat1	= GTCACAAGGGCCGGCCACAACTCGAG  # first constant part of forward mapping read
map_pat2	= TGATC  # second constant part of forward mapping read
map_pat_rev	= GTACGTCACAATATGATTATCTTTCTAGGGTTAA  # constant part of reverse mapping read
cores		= 2  # the number of processors to be used for Bowtie2 alignments 
bowtie_base	= /Users/wa/CoolShit/bowtie2-2.1.0/mm9  # the index of the genome to align against
max_dist_for	= 500  # max distance to cluster forward mapping positions
max_dist_rev	= 20  # max distance to cluster reverse mapping positions
min_counts	= 5  # minimum number of counts for a genuine barcode
		
EOF

exit 1;
}



##||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||##
##**************************************************************************************##
######				  	Sub-routines for barcode counts analysis 					######
##**************************************************************************************##
##||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||##


##########################################################################################
#################************  SUBROUTINE bc_extract_exp1  ****************###############
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
# It takes 8 arguments 
# 	a) $fastq -> The name of the fastq file to be read (with full path)
# 	b) $type -> Either "norm" (for normalization reads) or "exp" for expression reads. 
# 	c) $ind_len -> The length of the index sequence. Zero if no index is part of the read. (It
# 				   is important to note that nothing is done with index sequence per se. This 
# 				   is just to get the coordinates of pattern match right.
# 	d) $bc_len -> The length of the barcode.
# 	e) $pat1 -> The sequence of the first constant part.
# 	f) $pat2 -> The sequence of the second constant part.
# 	g) $BCs -> The reference to the empty hash (or hash with gDNA counts if the type is cDNA) 
# 				in which the barcodes and their counts are stored.
# 	h) $idx -> the position of the anonymous array in the hash where to add the count. 
# 				For $type eq "norm" it should be 0 and for others it should be 1 ..
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
# 	TAAAGTCCTGGCTTAA => [
# 						345,
# 						201
# 						]
##########################################################################################

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
	return qr/$regex/;   # this qr makes a regex object which can be directly used
}



##########################################################################################
#################***************   SUBROUTINE bc_group  ******************################
##########################################################################################
# DESCRIPTION: 
# A subroutine to group barcodes in a hash table where keys are barcodes and values are their
# counts.
# 
# INPUT: 
# It takes two arguments 
# 	a) $BCs  -> A reference to the hash that contains barcodes as keys and their counts as values
#				What is really IMPORTANT here is that read counts are the first element of an
#				anonymous array. this is like 
# 					%hash = {
#           				'CTAGACCATAGATTCA' => [
#                  		   				              1
#                  					               ],
#                			 }
# 
# 	b) $hd   -> Hamming distance threshould
# 	
# DEPENDS: 
# It depends on 
# Text::LevenshteinXS qw(distance)
# 
# WORKING: 
# Starting from the most frequent barcodes it finds and removes all the barcodes with 
# a hamming distance lower than hd (the second argument to subroutine) + 1
# OUTPUT: 
# The hash with only genuine barcodes and their counts as the first 
##########################################################################################

sub bc_group {
	my ($BCs, $hd) = @_ ;
	my @bc_array = sort { $BCs -> {$b}->[0] <=> $BCs -> {$a}->[0] } keys %$BCs ; #keys sorted based on values
	my @new_array;   # an empty array

	do {
		my $bc = shift @bc_array; # remove from old array
		push (@new_array , $bc) ; # add to the new array
		@bc_array = grep { ! (distance ($bc, $_) < ($hd + 1)) } @bc_array;

	} while (@bc_array) ;
	my (%BCs) = map { $_ => $BCs ->{$_} } @new_array;

	return %BCs;
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


##########################################################################################
##############************   	 SUBROUTINE time_suffix    ****************##############
##########################################################################################
# DESCRIPTION: 
# It creates a date time stamp, which can be added to filenmaes. 
# 
# INPUT: 
# It takes no arguments
# 
# DEPENDS: 
# It depends on DateTime;
# 
# WORKING: 
# Date time string is created for the current time using DateTime module function now
# This string is further modified to give a date time suffix
# OUTPUT: 
# Returns a string with date and time like this: 130817_165114
##########################################################################################

sub time_suffix {
	require DateTime;
	my $ts = DateTime->now(time_zone => 'local');
	$ts =~ s/[-:]//g;
	$ts =~ s/T/_/;
	return substr ($ts, 2, 13);
}


##||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||##
##**************************************************************************************##
######				 	 	Sub-routines for mapping analysis 						######
##**************************************************************************************##
##||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||##


##########################################################################################
#################*************** SUBROUTINE mapping_for1  ******************##############
##########################################################################################
# DESCRIPTION: 
# A subroutine to read the fastq file (typically the forward reads file from a paired End Illumina
# run. It reads the file and then extracts the genomic DNA and barcode sequence based on the given
# arguments and writes a fasta file where id is the variable part of read id joined with the 
# barcode sequence with "_" whereas the sequence is the genomic DNA component of the read.
# 
# INPUT: 
# It takes eight arguments 
# 	a) $fastq -> The name of the fastq file to be read (with full path)
# 	b) $ind_len -> The length of the index sequence. Zero if no index is part of the read. (It
# 				   is important to note that nothing is done with index sequence per se. This 
# 				   is just to get the coordinates of pattern match right.
# 	c) $bc_len -> The length of the barcode.
# 	d) $pat1 -> The sequence of the first constant part.
# 	e) $pat2 -> The sequence of the second constant part.
# 	f) $out_table -> The name of the file to write fasta (should be with full path)
# 	g) $IDs -> The reference to the empty hash in which the variable parts of the sequence ids 
# 				(as keys) and barcode sequences (as values) are deposited.
#	h) $BCs -> The reference to hash BCs coming from norm/exp data. It is used to check if 
#			   the read contains a genuine barcode or not. If yes the genuine barcode counter
# 			   gets one added to it.
# 
# DEPENDS: 
# It depends on 
# 	Hdist (see below)
# 	make_regex (see above)
# 
# WORKING: 
# It first determines the length of the $pat1 and the $pat2 (constant parts of the read). 
# It generates the regular expression for matching using subroutine make_regex. 
# One example of such an regexp with the $bc_length being 16 would look like this: 	
# 	^.{24,28}ACAACTCGAG(.{15,17})TGATC(.+)
# Then it reads the fastq file line by line. The first line should start with @ otherwise it will die. 
# The variable part of the Illumina id is extracted. Then it reads in the next line which is 
# the sequence read and matches the regexp pattern. 
# If the pattern is matched then the barcode and the genomic sequence is extracted using $1 
# and $1 variables (they store the recently matched patterns).  
# If the regexp pattern is not matched then, it finds the Hamming distance using Hdist subroutine between
# the $pat1 and the corresponding part of the read. If this is below (length $pat1)/7 then same 
# procedure is repeated for $pat2 and the correponding part of the read. Because the variance of
# one nucleotide in the barcode length is allowed (see above) so the $pat2 is matched at three
# possible places in the read. If the Hdist is below (length $pat2)/7 then the barcode and the
# genomic sequence are extracted. 
# In the end the id and the barcode sequence is fed to the hash IDs. And the id_barcode as the 
# name and the genomic DNA as the sequence is printed onto the fasta file $out_table. 
# During the the procedure the subroutine also takes record of total reads and the reads
# which could be matched properly. 
# 
# OUTPUT: 
# It returns the total reads and the reads matched.
# It populates the hash IDs with unique parts of sequence ids of matched reads (as keys) and 
# the barcode sequences (as values).
# Writes the fasta file from reads with successful patten matches. The name of each fasta entry 
# is the unique part of id and the barcode sequence separated by "_". The sequence is the 
# genomic DNA. One such entry is shown below
# >8:1101:3220:1977_CATATTCCGGCCTGAT
# CAGAATCCCAGGTAGTCCAGGCTGACCTCAGACTTACTTACGT
##########################################################################################

sub mapping_for1 {
	
	my ($fastq, $ind_len , $bc_len , $pat1 , $pat2, $out_table , $IDs, $BCs) = @_;
	
	#some initializations
	my $len1 = length $pat1;	# get the length of pat1
	my $len2 = length $pat2;	# get the length of pat2
	my $regex = make_regex ($ind_len, $bc_len, $pat1, $pat2, "map");
	my $total = 0; 		# total number of reads
	my $hits = 0;		# reads with a matched pattern
	my $genuine = 0;	# reads with a genuine barcode
	
	
	if ($debug) {
		print $error_fh "\nThe regex for $fastq is: $regex\n\n";
	}
	
	die "Cannot open the file $fastq\n" unless open (FASTQ , "<" , $fastq);
	open (FOUT , ">" , $out_table);

	while (my $line = <FASTQ>) {
	
		$total++;
		die "It is not a proper fastq file\n" if ($line !~ /^@.*:([^\s:]+:[^\s:]+:[^\s:]+:[^\s:]+):?(\s.*)?/);
		my $id_this = $1;
		
		$ line = <FASTQ>;
		chomp ($line);
		my $BC_this = 0;
		my $genomic_this;
		 
		if ($line =~ $regex) { 
		
			$BC_this =  $1 ;
			$genomic_this = $2 ; 
		
		
	
		} elsif ( Hdist ($pat1, substr($line, $ind_len, $len1)) <= ($len1/7) ) {
			
				if ( Hdist ($pat2, substr($line, ($ind_len + $len1 + $bc_len) , $len2)) < ($len2/7) ) {
				
					$BC_this = substr ($line , ($ind_len + $len1)  , $bc_len);
					$genomic_this = substr($line , ($ind_len + $len1 + $len2 + $bc_len) );
			
				} elsif ( Hdist ($pat2, substr($line, ($ind_len + $len1 + $bc_len -1) , $len2)) < ($len2/7) ) {
				
					$BC_this = substr ($line , ($ind_len + $len1)  , $bc_len - 1);
					$genomic_this = substr($line , ($ind_len + $len1 + $len2 + $bc_len -1) );
			
				} elsif ( Hdist ($pat2, substr($line, ($ind_len + $len1 + $bc_len + 1), $len2)) < ($len2/7) ) {
				
					$BC_this = substr ($line , ($ind_len + $len1)  , $bc_len + 1);
					$genomic_this = substr($line , ($ind_len + $len1 + $len2 + $bc_len + 1) );
				}
				
		}
		if ($BC_this) {
			$hits++ ;
			if (exists $BCs -> {$BC_this} ) { $genuine++ }
			
			$IDs -> {$id_this} = $BC_this ;
			print FOUT  ">" . $id_this . "_" . $BC_this . "\n" . $genomic_this . "\n";
		}
	
	
		<FASTQ>;
		<FASTQ>;
	}

	close FASTQ;
	close FOUT ;
	return ($total , $hits, $genuine);
}


##########################################################################################
#################*************** SUBROUTINE fastq_subset  ******************##############
##########################################################################################
# DESCRIPTION: 
# A subroutine to read the fastq file (typically the reverse reads file from a paired End Illumina
# run. It reads the file and then writes a fasta file for only those IDs which are present in the
# hash of IDs given as one of the arguments. It also removes the constant part from the read before
# writing the fasta file.  
# 
# INPUT: 
# It takes one argument 
# 	a) $fastq -> The name of the fastq file to be read (with full path)
# 	b) $pat  ->  The constant part (DNA string) to be matched and removed from the start of the reads
# 	c) $out_table -> The name of the file to write fasta (should be with full path)
# 	d) $IDs -> The reference to the hash of IDs (to subset the fastq file). This hash is the 
# 			   output of subroutine mapping_for1 (see above).
# 
# DEPENDS: 
# It depends on Hdist (see below) 
# 
# WORKING: 
# It first determines the length of the $pat (constant part of the read). If it is larger than 10
# the pattern is reduced to last 10 characters of $pat. Then it makes a regexp string based on $pat
# assuming that the $pat is at the start of the string. In the regexp, all the nucleotides after 
# the match are captured by using (.+). 
# Then it reads the fastq file line by line. The first line should start with @ otherwise it will die. 
# The variable part of the Illumina id is extracted. If this id exists in the hash IDs then it proceeds
# further otherwise reads in the next three lines and goes to the next round of while loop. In case 
# this id exists, it reads in the next line which is the sequence read and matches the regexp pattern. 
# If the pattern is matched then the extracted piece of read becomes the sequence for
# the fasta entry. The id for fasta is made by joining > , $id_this , "_" and the barcode sequence
# which is the value of $id_this in the hash IDs. 
# If the regexp pattern is not matched then, it finds the Hamming distance using Hdist subroutine between
# the $pat and the first length ($pat) characters of the read. If this is below (length $pat)/7 then the 
# rest of the read is taken out as the read entry for the fasta file.  
# 
# OUTPUT: 
# Writes the fasta file containing the genomic parts of reads of those Fastq entries whose ids are
# present in the hash IDs. The name of each fasta entry is the unique part of id and the barcode
# sequence separated by "_". One such entry is shown below
# >8:1101:3220:1977_CATATTCCGGCCTGAT
# CCCCAGCACTTGAGAAAGACAAGACAGCAGCATCAGGAGCTCAAGGAGGGCAGCTTTGAATACGTA
##########################################################################################

sub fastq_subset {
	my ($fastq, $pat, $out_table , $IDs) = @_;
	
	my $len = length $pat;	# get the length of pat
	my $short_pat = $pat;   # if the length of pat is larger than 10 then only last 10 bases are taken for matching
	if ( $len  > 10 ) {  
		$short_pat = substr($pat, ($len -10));
	} 
	my $regex = "^.{" . ($len - (length $short_pat) - 2) . "," . ($len - (length $short_pat) + 2) . "}" . 
				$short_pat . "(.+)" ;

	
	if ($debug) {
		print $error_fh "\nThe regex for $fastq is: $regex\n\n";
	}

	die "Cannot open the file $fastq\n" unless open (FASTQ , "<" , $fastq);
	open (FOUT , ">" , $out_table);
	
	while (my $line = <FASTQ>) {
		chomp ($line);
		
		die "It is not a proper fastq file\n" if ($line !~ /^@.*:([^\s:]+:[^\s:]+:[^\s:]+:[^\s:]+):?(\s.*)?/);
		my $id_this = $1;
		
		if (exists $IDs -> {$id_this}) {
			$line = <FASTQ>;
			chomp ($line);
			if ($line =~ /$regex/) {
				print FOUT ">" . $id_this . "_" . $IDs -> {$id_this} ."\n" . $1. "\n";
			}  elsif (Hdist ($pat, substr ($line, 0 , $len)) <  $len/7 ){
					print FOUT ">" . $id_this . "_" . $IDs -> {$id_this} . "\n" . substr ($line, $len) . "\n";	
 				}
		} else { <FASTQ> }
		
		<FASTQ>;
		<FASTQ>;
	}
	close FASTQ;
	close FOUT ;
}


##########################################################################################
#################**************** SUBROUTINE top_map ******************###################
##########################################################################################
# DESCRIPTION: 
# A subroutine to work on a hash of hashes (the output of subroutine parse_sam) it analyzes 
# the mapping data and then tells the top location, its frequency and other mapping statistics 
# which are put back into the same hash, thus destroying the old hash and creating a new one.
# 
# INPUT: 
# It takes two arguments 
# 	a) $map_hash -> a reference to a hash of hashes (the output of parse_sam subroutine). The 
# outer hash has the keys that are barcodes. The value of each key (barcode) is another hash
# whose keys are position strings (chr:ori:pos) and the values are anonymous arrays each
# consisting of two elements [read counts, total mapq]
# 	b) $max_dist -> an intiger value which is the max distance within which the mapped positions
# are merged together (the position with most reads now represents all the reads within 
# that distance)
# 
# DEPENDS: 
# It depends on 
# 	List::Util qw(sum), 
# 	top_kv (see below),
# 	map_refine (see below) 
# 
# WORKING: 
# It loops through each key (barcode in this case) of the hash
# 	For each barcode, if there is only one mapping information (that is the value is a hash 
# with only one key) then the key is split into chr, ori and pos. The value which is an 
# anonymous array, its first element becomes the total reads by dividing the second member 
# (which is sum of all the mapq scores) with total reads, one gets average mapq. The freq1 
# (the frequency of most frequent location) becomes one and the freq2 becomes zero.
# 	Now if there are more than one location mapped, then total reads are the sum of first 
# elements of all the values (anonymous arrays) in this hash total mapq is the sum of second
# elements of all the values in this hash.
# 	Now using the top_kv subroutine the top location, reads associated with the top location
# and the total mapq associated with the top location are determined. This subroutine also 
# removes the top location as a key from the hash, so hash is reduced in size.
# 	Then using map_refine subroutine the remaining keys in the hash are looked for locations 
# mapped closely (within $max_dist) to the top location. If this is the case then their reads and 
# mapq is added to the those of top location and those closely mapped keys are removed from the
# hash. At this stage the chr, ori and pos are extracted for the top location. Orientation is
# changed to "+" if it is 0, or to "-" if it is 16. Also average mapq (av_mapq) and freq1 are 
# calculated for top location. 
# 	Now if after all that, the hash is empty then freq2 becomes zero, otherwise second top 
# location is determined again using the top_kv subroutine and map_refine to finally get to 
# freq2. In the end the value of the current key of the hash is changed to a an anonymous array 
# consisting of all the mapping parameters.
# 
# OUTPUT: 
# Changes the input hash of hashes into a simple hash whose keys are barcodes and the value of 
# each key (barcode) is an anonymous array consisting of mapping information
# [chr, ori, pos, total reads, average mapq, freq1 and freq2]   
##########################################################################################

sub top_map {
 	my ($map_hash, $max_dist) = @_;
	foreach my $bc (keys %$map_hash) {
		my %this = %{ $map_hash -> {$bc} };
		my ($chr, $ori, $pos, $t_reads, $av_mapq, $freq1, $freq2);
		# if there is only one key in the hash, then it is simple, the only location is the top location	
		if (scalar keys %this == 1 ) {
			($chr, $ori, $pos) = map { split (":" , $_)  } keys %this;
			if ($ori == 0) { $ori = "+" } elsif ($ori == 16) { $ori = "-" }
			
			($t_reads, $av_mapq) =  map { @$_ } values %this;
			$av_mapq =  sprintf "%.2f" , $av_mapq/ $t_reads;
			$freq1 = 1.000;
			$freq2 = 0.000;
		} else {
				$t_reads = sum map { $_ ->[0] } values %this;
				my $t_mapq  = sum map { $_ ->[1] } values %this;
				my ($top1) = top_key (\%this); # finds top location and then refines it
				my ($reads1, $t_mapq1) = map_refine ($top1, \%this, $max_dist);
				($chr, $ori, $pos) = split (":" , $top1);
				
				if ($ori == 0) { $ori = "+" } elsif ($ori == 16) { $ori = "-" }
				$av_mapq = sprintf "%.2f" , $t_mapq1 / $reads1;
				$freq1   = sprintf "%.3f" , $reads1  / $t_reads;
				
				# If after map_refine, no more keys are left in hash then you do not have 2nd location
				my ($reads2, $top2, $t_mapq2);
				if (scalar keys %this < 1) {
					$freq2 = 0.000;
				} else {
						($top2) = top_key (\%this);
						($reads2, $t_mapq2) = map_refine ($top2, \%this, $max_dist);
						$freq2   = sprintf "%.3f" , $reads2  / $t_reads;
					}
		}
				
		$map_hash -> {$bc} = [$chr, $ori, $pos, $t_reads, $av_mapq, $freq1, $freq2];	
	}
 }



##########################################################################################
#################**************** SUBROUTINE map_refine ******************################
##########################################################################################
# DESCRIPTION: 
# A subroutine to work on a hash whose keys are genomic loation strings (chr:ori:pos) and values
# are anonymous arrays with two element in each [read counts, total mapq]. Using a reference key 
# given as an argument it finds other keys in the hash which represent a location within the 
# max_dist --another arugment-- of the given reference key. If such keys are found those keys are 
# removed from the hash and their values [read counts, total mapq] are added to the $read and $mapq
# arguments given (which actually are the values of the reference key). In the end it returns 
# the updated values of $read and $mapq with all the closely associated loations (keys) removed
# from the hash. 
# 
# INPUT: 
# It takes three arguments 
# 	a) $top -> It is a genomic loation string (chr:ori:pos) and is actually one of the keys in $hash
# 	b) $hash -> A reference to the hash on which this subroutine works.
# 	c) $max_dist -> An intiger value which is the max distance within which the mapped positions
# are merged together (the position with most reads now represents all the reads within 
# that distance)
# 
# DEPENDS: 
# no dependencies 
# 
# WORKING: 
# Using the reference key $top it gets the $reads and $mapq (the total reads and the map quality 
# associated with the top location). Then deletes the $top key from the hash.
# converts the $top key into an array of [chr, ori, pos] using split.
# Also converts the keys of the hash (now without $top key) into an array of arrays where each key
# is now ["chr:ori:pos"(the key) , chr, ori, pos]. 
# Then using grep it subsets the array (of arrays) to only those for which the location falls within
# $max_dist. 
# Finally it iterates through this subsetted array and using the first element of each sub-array 
# which is the key, it takes the read counts and map quality values from the hash and adds them to $reads
# and $mapq respectively. Also removes those keys from the hash.
#  
# OUTPUT: 
# Changes the input hash by deleting the $top key and possibly other keys if their positions are 
# within max_dist. Returns the updated values of $reads and $mapq. 
##########################################################################################

sub map_refine {
	my ($top, $hash, $max_dist) = @_;
	my ($reads, $mapq) = @{ $hash ->{$top} }; # using reference key get $reads and $mapq
	delete $hash ->{$top};	# remove the reference key from the hash
	my @top_arr = split (":" , $top);
	
	# make an array of arrays with each element being ["chr:ori:pos"(the key) , chr, ori, pos]
	my @positions = map { [ $_ , split (":" , $_) ] } keys %$hash;
	
	# filter out  those elements of this array where the location is not within max distance
	@positions = grep { $_ ->[1] eq $top_arr[0] and 
						$_ ->[2] == $top_arr[1] and
						abs ($_ ->[3] - $top_arr[2]) < $max_dist } @positions;
	
	# iterate through this array and using the first element of the sub-array (which is the key of the hash)
	# transfer the read counts and mapq to the $read and $mapq
	# also delete the corresponding keys form the hash
	foreach my $k (@positions) {
		$reads  += $hash ->{ $k ->[0] } ->[0];
		$mapq   += $hash ->{ $k ->[0] } ->[1];
		delete  $hash ->{ $k -> [0] };
	}
	return ($reads, $mapq);
}


##########################################################################################
#################**************** SUBROUTINE top_key ******************###################
##########################################################################################
# DESCRIPTION: 
# A subroutine that works on a hash and finds the key with highest value of first element of 
# the value (which is an anonymous array in this case)
# 
# INPUT: 
# It takes one argument 
# 	a) $hash -> A reference to the hash on which this subroutine works. Each key of this hash
# 				has a value which is an anonymous array to two elements. The number in the first
# 				element is used to find the key with the highest value.
# 
# DEPENDS: 
# no dependencies 
# 
# WORKING: 
# Iterates through the keys and values of the hash and updates the $top_key and $top_val variables
# if the value at [0] of the anonymous array is larger than the current $top_val. According to 
# Stack Overflow this is faster than just simply sorting the hash based on values.
# http://stackoverflow.com/questions/2886872/what-is-the-easiest-way-to-get-a-key-with-the-highest-value-from-a-hash-in-perl
# #
# OUTPUT: 
# Returns the the key with the highest value at [0] of the anonymous array 
##########################################################################################

sub top_key {
    my ($hash)  = @_;
    my $top_val = 0;
    my $top_key;
    while ((my $key, my @value) = each %$hash) {
    	#print Dumper \@value;
  		if ($value[0][0] > $top_val) {
    		$top_val = $value[0][0];
    		$top_key = $key;
  		}
  	}
  	return($top_key);
}



##########################################################################################
#################**************** SUBROUTINE parse_sam ******************#################
##########################################################################################
# DESCRIPTION: 
# A subroutine that reads in a sam file (the output of Bowtie alignment) line by line and 
# then stores the information into a hash of hashes (whose reference is given as the second 
# argument to the subroutine). 
# 
# INPUT: 
# It takes two arguments 
# 	a) $samFile -> the name of the same file to be read. The name should contain the whole path
# 	b) $map_hash -> A reference to the empty hash to which mapping information read from sam file
# 
# DEPENDS: 
# It depends on
# 	cigar_parse (see below)
# 
# WORKING: 
# It initializes a hash %cigar_hash for efficiently parsing the CIGAR string (sub cigar_parse)
# Everytime a new CIGAR string comes up it is parsed using cigar_parse and the result is stored 
# in a hash (CIGAR string => result). Next time when this CIGAR string comes up, then you do not
# need to parse it again, you just do a lookup of %cigar_hash and get the result from the value.
# 
# The sam file is read line by line skiping the header lines as they start with @.
# Extracts $id (ID) and other mapping features. 
# The $id is composed of read_ID and the barcode joined by "_". The barcode is extracted as $BC_this.
# $pos (the position of alignment to genome) is adjusted based on orientation and cigar. 
# Now the mapping location is joined into one string $mapping (chr:ori:pos).
# Then the $BC_this (the barcode) is looked up in the map_hash if it exists. If does then within
# the value of the key $BC_this, $mapping (chr:ori:pos) is looked up. If it exists then you add
# 1 to the first element of the anonymous array (this is read count) and add the $mapq (the alignment
# quality score) to the second element of this anonymous array. If the key does not exists then
# you initilize it. The resultant hash looks something like this.  
# 
# %map_hash = ('TAAAGTCCTGGCTTAA' => {
#                                   'chr16:0:92423745' => [
#                                                           8,
#                                                          '42'
#                                                         ]
#                                 },
#           'AATGCTAGCACTAAGC' => {
#                                   'chr19:16:32841919' => [
#                                                            100,
#                                                            '40'
#                                                          ],
#                                 
#                                  'chr19:16:32841999' => [
#                                                            1,
#                                                            '40'
#                                                          ],
#                                 'chr19:16:32841929' => [
#                                                            101,
#                                                            '38'
#                                                          ]
#                                 );
# 
# 
# 
# 
# OUTPUT: 
# The hash $map_hash (which was empty to start with) is now filled up with all the mapping 
# information as shown above.
##########################################################################################

sub parse_sam {
	my ($samFile, $map_hash, $re_map) = @_;
	my %cigar_hash; # a hash for efficiently parsing cigar strings
	my $mapping; # initialization of the string that harbors the mapping info (chr:ori:pos)
	my $BC_this; # the barcode for each read which would become the key of our hash %map_hash
	
	open (SAM , "<" , $samFile) or die "Cannot open the sam file $samFile: $!\n";
	
	while (<SAM>) {
		next if (substr $_, 0, 1) eq "@";  # skip all the header lines
		my ($id, $ori, $chr, $pos, $mapq, $cigar, undef, undef, undef, $query) = (split /\s/) [0..9];
		
		$BC_this = (split ("_" , $id))[-1];
		
		if ($cigar =~ /^\d{2,}M/ or $cigar =~ /^[12345]S/) {
			#if the read was aligned to -ve strand then you need to find the appropriate start pos of the read
			#for this you have to add the span (calculated from Cigar) - 1 to the -ve strand reads 
			if ($ori == 16) {			  # for -ve strand reads
				if (exists $cigar_hash {$cigar} ) {
				$pos += $cigar_hash {$cigar} - 1;
				} else {
						$cigar_hash {$cigar} = cigar_parse ($cigar, "MD" , "SI") ;
						$pos += $cigar_hash {$cigar} - 1;
					
					}
			}
			$mapping = join (":" , ($chr, $ori, $pos));
					
		} elsif ($cigar =~ /^\*$/) {
			$mapping = join (":" , ($chr, $ori, $pos));
			} elsif ($cigar =~ /^(\d+)S/ and $1 > 17) { # these are cases which need re-alignmnet
				$re_map -> {$id} = substr ($query, 0 , $1);
				}
		
		# put the information in the map_hash
		if (exists $map_hash -> {$BC_this} ) {
			if (exists $map_hash -> {$BC_this} -> {$mapping} ) {
				$map_hash -> {$BC_this} -> {$mapping} -> [0]++;
				$map_hash -> {$BC_this} -> {$mapping} -> [1] += $mapq;
			} else {
					$map_hash -> {$BC_this} -> {$mapping} = [1 , $mapq];
				} 
		} else {
				$map_hash -> {$BC_this} =  { $mapping => [1, $mapq] }; 
			}		
	}
	return;
}



##########################################################################################
#################*************** SUBROUTINE cigar_parse  ******************###############
##########################################################################################
# DESCRIPTION: 
# A subroutine to parse the CIGAR string from sam file output. The purpose is to get a value which
# is equal to the length of the genomic sequence span of the aligned region. 
# 	So you include M (matched) and D (deleted in the query) numbers and 
# exclude S (soft clipped) and I (insertion in the query)  numbers. Thus the sum of M and D 
# numbers is what is required. For instance from cigar 34M5S2I26D the result should be 60. 
# 
# INPUT: 
# It takes three arguments 
# 	a) $cigar -> The CIGAR string from sam file output
# 	b) $keep  -> The letters of CIGAR string whose numbers we are interested in. Typically M and D
# 				 for our mapping purpose.
# 	c) $discard -> The letters of CIGAR string that we are not interested in. Typically S and I.
# 
# DEPENDS: 
# It depends on 
# 	List::Util qw(sum)
# 
# WORKING: 
# First it checks the $cigar string if it is a valid entry: containing only [MSID]
# If no pattern is given for $keep, it complains and dies.
# If there is someting given for $discard, then if will remove those letters and associated 
# numbers from the CIGAR string. The remaining string will be splitted with the split being
# the letters in $keep. These leaves behind only an array of intigers which are summed together
# using List::Util qw(sum).
# OUTPUT: 
# Returns the sum of the numbers associated with the CIGAR letters given in $keep
##########################################################################################
sub cigar_parse {
	use List::Util qw(sum);
	my ($cigar, $keep , $discard) = @_;
	
	die "$cigar is not a valid CIGAR string\n" if ($cigar =~ /[^\dMSID]/);
	die "No pattern is given to look for\n" unless $keep;
		
	if ($discard) {
		$cigar =~ s/\d+[$discard]//g;
		return sum split ("[$keep]", $cigar);
	} else {
		return sum split ("[MSDI]" , $cigar);
		}
}







