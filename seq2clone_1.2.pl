#!/usr/bin/perl
#use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;

#seq2clone.pl - identify unique clones from sequence reads
#             - output position, unique clone name, #sequence reads, #aligned locations, #alignments at each position (for checking)
#version: 1.2: post-processing of capture regions alignment data

my $usage="SYNOPSIS: seq2clone.pl [-p] -f INPUTBEDFILE
OPTIONS: 
	-p : post processing (process the number of alignment in capReg)
	 -f : BED file\n";

my $postProcess;
my $inputBedFilename;

my %options=();
getopts("pf:", \%options);

if (defined($options{f}) && defined($options{p})) {
	$inputBedFilename = $options{f};
	$postProcess = 1;
} elsif (defined $options{f}) {
	$inputBedFilename = $options{f};
	$postProcess = 0;
} else {
	die("$usage");
}

# Define BED file attributes (only up to strand; other attributes not used):
my @BED_attributes = ('chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand');
# Define mod_BED file attributes + # alignment in capReg:
my @BED_capReg = ('chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', '', '', 'capReg');
# Define BEDgraph file attributes (only up to strand; other attributes not used):
my @BEDgraph_attributes = ('chrom', 'chromStart', 'chromEnd', 'score');

# Open & read input file:
my $input_file = &get_file_handle( $inputBedFilename );

my $A_headerLine = "";
my @input_data = &get_data ( $postProcess, $input_file, \$A_headerLine );
close( $input_file );

# hash of seqs pointing to clone_ids:
my %seqs;
# hash of locations pointing to clone_ids & number of seqs:
my %locs;
# array of clone_ids pointing to parent clone_id, number of seqs & number of locations:
my @clones;
# total number of clone_ids assigned so far:
my $nclone_ids = 0;

my $location;

for ($i=0;$i<=$#input_data;$i++) {
	$location = $input_data[$i]{'chrom'}."\t".$input_data[$i]{'chromStart'}."\t".$input_data[$i]{'chromEnd'}."\t".$input_data[$i]{'strand'};
    $locs{ $location }{'nseqs'}++;
	# print($input_data[$i]);
    my $clone_id;
    if ( defined( $seqs{ $input_data[$i]{'name'} } ) ) {
        $clone_id = $seqs{ $input_data[$i]{'name'} };
        if ( defined( $locs{ $location }{'clone_id'} ) ) {
            # check for colliding clone_id:
            if ( $locs{ $location }{'clone_id'} != $clone_id ) {
                # replace clone_id with most-upstream parent_clone_id (or leave it as current clone_id if not yet assigned):
                while ( defined( $clones[ $clone_id ]{'parent_clone_id'} ) ) { $clone_id = $clones[ $clone_id ]{'parent_clone_id'} }
                my $colliding_clone_id = $locs{ $location }{'clone_id'};
                # replace colliding_clone_id with most-upstream parent clone_id for colliding clone (or leave it as current colliding clone_id if not yet assigned):
                while ( defined( $clones[ $colliding_clone_id ]{'parent_clone_id'} ) ) { $colliding_clone_id = $clones[ $colliding_clone_id ]{'parent_clone_id'} }
                # assign parent_clone_id for higher-numbered upstream clone_id
                if ( $clone_id > $colliding_clone_id ) {
                    $clones[ $clone_id ]{'parent_clone_id'} = $colliding_clone_id;
                }
                elsif ( $colliding_clone_id > $clone_id ) {
                    $clones[ $colliding_clone_id ]{'parent_clone_id'} = $clone_id;
                }
            }
        }
        else {
            $locs{ $location }{'clone_id'} = $clone_id;
            $clones[ $clone_id ]{'nlocs'}++;
			if (defined ($input_data[$i]{'capReg'})) {
				$clones[ $clone_id ]{'capReg'}+=$input_data[$i]{'capReg'};
			}
		}
	}
    elsif ( defined( $locs{ $location }{'clone_id'} ) ) {
        $clone_id = $locs{ $location }{'clone_id'};
        $seqs{ $input_data[$i]{'name'} } = $clone_id;
        $clones[ $clone_id ]{'nseqs'}++;
    }
    else {
        $clone_id = ++$nclone_ids;
        $seqs{ $input_data[$i]{'name'} } = $clone_id;
        $clones[ $clone_id ]{'nseqs'}++;
        $locs{ $location }{'clone_id'} = $clone_id;
        $clones[ $clone_id ]{'nlocs'}++;
		if ( defined ($input_data[$i]{'capReg'}) ) { 
			$clones[ $clone_id ]{'capReg'}+=$input_data[$i]{'capReg'};
		}
	}
}

# now ensure parent_clone_id points to most-upstream parent, update parent clone data, & assign consective output_ids:
my $clone_index=0;
for ($i=0; $i<=$#clones; $i++ ) {
    if ( defined( $parent_clone_id = $clones[ $i ]{'parent_clone_id'} ) ) {
        # find most-upstream parent_clone_id
        while ( defined( $clones[ $parent_clone_id ]{'parent_clone_id'} ) ){ $parent_clone_id = $clones[ $parent_clone_id ]{'parent_clone_id'} };
        $clones[ $i ]{'parent_clone_id'} = $parent_clone_id;
        $clones[ $parent_clone_id ]{'nseqs'}+= $clones[ $i ]{'nseqs'};
        $clones[ $parent_clone_id ]{'nlocs'}+= $clones[ $i ]{'nlocs'};
		if ($postProcess==1) { $clones[ $parent_clone_id ]{'capReg'}+= $clones[ $i ]{'capReg'} };
    }
    else { $clones[ $i ]{'output_id'} = $clone_index++ }
}

# Write output data
foreach $loc ( sort keys %locs ) {
	my $capReg;
	my $clone_id = $locs{ $loc }{'clone_id'};
	@location = split( /\t/, $loc);
	my $chrom = $location[0];
	my $chromStart = $location[1];
	my $chromEnd = $location[2];
	my $strand = $location[3];
	if ($postProcess==1) {
		if ( defined( $clones[ $clone_id ]{'parent_clone_id'} ) ){ $clone_id = $clones[ $clone_id ]{'parent_clone_id'} };
		if ( defined( $clones[ $clone_id ]{'capReg'})) { $capReg = $clones[ $clone_id ]{'capReg'}};
		print "$chrom\t$chromStart\t$chromEnd\tclone_$clones[ $clone_id ]{'output_id'}\t1\t$strand\t$clones[$clone_id]{'nlocs'}\t$capReg\n";
	}
	else {
		if ( defined( $clones[ $clone_id ]{'parent_clone_id'} ) ){ $clone_id = $clones[ $clone_id ]{'parent_clone_id'} };
		print "$chrom\t$chromStart\t$chromEnd\tclone_$clones[ $clone_id ]{'output_id'}\t1\t$strand\t$clones[$clone_id]{'nseqs'}\t$clones[$clone_id]{'nlocs'}\t$locs{ $loc }{'nseqs'}\n";
	}
}

exit;


sub get_file_handle {
    
    my $filename = shift;
    my $file_handle;
    
    if ($filename eq 'stdin') { open ($file_handle, "<&STDIN") or die ("can't read from stdin\n") }
    else {
        if ($filename =~ /\.gz\/?\z/) {
            if (system("gzip -t $filename") == 0) {
                open ($file_handle, "gunzip -c $filename |") or die ("can't open $filename\n");
            }
            else { die ("gzip -t fail for $filename\n") }
        }
        else { open ($file_handle, "<$filename") or die ("can't open $filename\n") }
    }
    return $file_handle;
}


sub get_data {
    
    my @unsorted_data = (); #initialize as empty array (of hashes)

	my $postProcess = shift;

	my @attributes;
    
	if ($postProcess==1) {
		@attributes = @BED_capReg;
	}
	else {
		@attributes = @BED_attributes;
	}

    my $fileHandle = shift;
    my $headerLineRef = shift;

    my $line_number = 0;

    while ($currentLine = <$fileHandle>) {
        chomp ($currentLine);
        if ($currentLine =~ /^track/) { #if it's a track line determine the file type
            if ( defined( $headerLineRef )) { $$headerLineRef = $currentLine }; #save header line
            if ($currentLine =~ /type=bedgraph/i ) { @attributes = @BEDgraph_attributes }; # /i case-insensitive search
        }
        elsif ($currentLine =~ /^#/) { } #if it's a comment line (starting with '#') do nothing
        else { # data line
            my @values = split( /\s/, $currentLine); #split on whitespace \s
            my %param;
            my $index = 0;
            while ( ( $index <= $#values ) && ( $index <= $#attributes ) ) {
                $param{$attributes[$index++]} = $values[$index];
            }
            
            if ( defined( $param{'chrom'} ) && defined( $param{'chromStart'} ) && defined( $param{'chromEnd'} ) && defined( $param{'strand'} ) ) {
                %{ $unsorted_data[++$#unsorted_data] } = %param;
                #$unsorted_data[$#unsorted_data]{'line'} = $currentLine;
                #$unsorted_data[$#unsorted_data]{'line_number'} = $line_number++;
            }
        }
    }
    
    return @unsorted_data;
}


