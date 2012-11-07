#!/usr/bin/env perl
###############################################################################
#
#    parseCazyBlast
#    
#    Parse the output of a cazy based blast to determine best matching regions
#
#    Copyright (C) Michael Imelfort
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

#pragmas
use strict;
use warnings;

#core Perl modules
use Getopt::Long;
use Carp;
use POSIX;
use List::Util qw(min);

#CPAN modules
use Bio::SearchIO; 
use Data::Dumper;

#locally-written modules

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# edit here to log all external commands 
my $global_log_commands = 0;

# ext command failure levels
use constant {
    IGNORE_FAILURE => 0,
    WARN_ON_FAILURE => 1,
    DIE_ON_FAILURE => 2
};

# get input params and print copyright
printAtStart();
my $global_options = checkParams();

######################################################################
# CODE HERE
######################################################################
my $global_perc_match = overrideDefault(0.8, 'conserved_length');
# open the blast file
my $bformat = overrideDefault('blast', 'format');
my $in = new Bio::SearchIO(-format => $bformat, -file => $global_options->{'blast'} ) ;

# open the output file 
my $out_fh = openWrite($global_options->{'out'});

my %global_hit_hash = ();       # store only the best hits here. One for each region
my $region_id = 0;   

while( my $result = $in->next_result ) 
{
    my @hit_array = ();
    while( my $hit = $result->next_hit ) 
    {
        my $hsp_ok = 1; 
        my $some_hsps = 0;
        while( my $hsp = $hit->next_hsp ) 
        {
            $some_hsps = 1;
            if(defined $global_options->{'length'}) {
                if( $hsp->length('total') < $global_options->{'length'}) {
                    $hsp_ok = 0; next;
                }
            }
            if(defined $global_options->{'percid'}) {
                if ( $hsp->percent_identity < $global_options->{'percid'}) {
                    $hsp_ok = 0; next;
                }
            }
            if(defined $global_options->{'evalue'}) {
                if ( $hsp->significance > $global_options->{'evalue'}) {
                    $hsp_ok = 0; next;
                }
            }
            if(defined $global_options->{'bits'}) {   
                 if( $hsp->bits < $global_options->{'bits'}) {
                    $hsp_ok = 0; next;
                }
            }
            if(defined $global_options->{'mismatches'}) {
                my $mismatchcount = $hsp->length('total') - ($hsp->num_conserved + $hsp->gaps('total'));
                if ( $mismatchcount > $global_options->{'mismatches'}) {
                    $hsp_ok = 0; next;
                }
            }
            
            # we require that the hit be at least 80% as long as the reference
            if($hsp->length('hit')/$hsp->length('total') < $global_perc_match)
            {
                $hsp_ok = 0; next;
            }
            
            # print $result->query_name. " " .$hit->name."<--\n";
            # contig_1 lcl|ADY13293.1_Neopullulanase<--
        }
        if((1 == $hsp_ok) and (1 == $some_hsps))
        {
            # a keeper!
            push @hit_array, \$hit;
        }
    }

    # sort all the stored hits by starting pos on the contig
    # now the sort will work like we want it to!
    my @sorted_hit_array = sort {min(int(${$a}->{'_hsps'}[0]{'QUERY_START'}), int(${$a}->{'_hsps'}[0]{'QUERY_END'}))  <=> min(int(${$b}->{'_hsps'}[0]{'QUERY_START'}), int(${$b}->{'_hsps'}[0]{'QUERY_END'}))} @hit_array;

    my %current_best_hit_refs = (); # rid -> best hits for each region in this contig
    my %region_starts = ();    # rid -> region starts
    my %region_ends = ();      # rid -> region ends
    my $region_start = LONG_MAX;
    my $region_end = -1;
    my $region_best_eval = 100; 
    my $region_complete = 0;
    
    foreach my $hit_ref (@sorted_hit_array)
    {
        my $hit = ${$hit_ref};
        $hit->rewind;       # reset the internal counter or the folowing call to next_hsp will return NULL
        my $hsp = $hit->next_hsp;
        my $q_start = $hsp->start('query'); 
        my $q_end = $hsp->end('query');
        if($q_start > $q_end)
        {
            my $tmp = $q_start;
            $q_start = $q_end;
            $q_end = $tmp;
        }
        #print  $q_start."\t".$q_end."\n";
        if(($q_start > $region_end) and (-1 != $region_end))
        {
            # we have left the region. The best hit here is the best hit
            $region_start = LONG_MAX;
            $region_end = -1;
            $region_best_eval = 100; 
            $region_id++;
        }
        else
        {
            # update region boundaries
            if($q_start < $region_start)
            {
                $region_starts{$region_id} = $q_start; 
                $region_start = $q_start;
            }
            if($q_end > $region_end)
            {
                $region_ends{$region_id} = $q_end;
                $region_end = $q_end;
            }
            # test to see if we have a better match now
            if($hsp->significance < $region_best_eval) 
            {
                $current_best_hit_refs{$region_id} = \$hit;
                $region_best_eval = $hsp->significance;
            }
        }
    }
    
    # now we should have a list of best hits for this contig
    foreach my $rid (sort {$a <=> $b} (keys (%current_best_hit_refs)))
    {
        my $hit = ${$current_best_hit_refs{$rid}};
        $hit->rewind;       # reset the internal counter or the folowing call to next_hsp will return NULL
        my $hsp = $hit->next_hsp;
        print $out_fh join("\t", ($result->query_name, $hit->name, $hsp->start('query'), $hsp->end('query'), $hsp->strand('query'), $hsp->significance))."\n";
    }
    $region_id++;
}

close $out_fh;

######################################################################
# CUSTOM SUBS
######################################################################

######################################################################
# TEMPLATE SUBS

######################################################################
# PARAMETERS

sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "blast|b:s" ,"format|f:s", "length|l:i", "percid|p:i", "evalue|e:i", "bits|b:s", "mismatches|m:i", "conserved_length|c:i", "out|o:s" );
    my %options;

    # Add any other command line options, and the code to handle them
    # 
    GetOptions( \%options, @standard_options );

    # if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    # Compulsory items
    if(!exists $options{'blast'} ) { printParamError ("You need to supply a blast file"); }
    #if(!exists $options{''} ) { printParamError (""); }

    return \%options;
}

sub printParamError
{
    #-----
    # What to do if there's something wrong with a parameter
    #  
    my ($error) = @_;  
    print "**ERROR: $0 : $error\n"; exec("pod2usage $0");
}

sub overrideDefault
{
    #-----
    # Set and override default values for parameters
    #
    my ($default_value, $option_name) = @_;
    if(exists $global_options->{$option_name}) 
    {
        return $global_options->{$option_name};
    }
    return $default_value;
}

######################################################################
# FILE IO

sub openWrite
{
    #-----
    # Open a file for writing
    #
    my ($fn) = @_;
    open my $fh, ">", $fn or croak "**ERROR: could not open file: $fn for writing $!\n";
    return $fh;
}

sub openRead
{   
    #-----
    # Open a file for reading
    #
    my ($fn) = @_;
    open my $fh, "<", $fn or croak "**ERROR: could not open file: $fn for reading $!\n";
    return $fh;
}

######################################################################
# EXTERNAL COMMANDS
#
# checkAndRunCommand("ls", {
#                          -a => ""
#                          }, 
#                          WARN_ON_FAILURE);

sub checkFileExists {
    #-----
    # Does a file exists?
    #
    my ($file) = @_;
    unless(-e $file) {
        croak "**ERROR: $0 : Cannot find:\n$file\n";
    }
}

sub logExternalCommand
{
    #-----
    # Log a command line command to the command line!
    #
    if(1 == $global_log_commands) {
        print $_[0], "\n";
    }
}

sub isCommandInPath
{
    #-----
    # Is this command in the path?
    #
    my ($cmd, $failure_type) = @_;
    if (system("which $cmd |> /dev/null")) {
        handleCommandFailure($cmd, $failure_type);
    }
}

sub runExternalCommand
{
    #-----
    # Run a command line command on the command line!
    #
    my ($cmd) = @_;
    logExternalCommand($cmd);
    system($cmd);
}

sub checkAndRunCommand
{
    #-----
    # Run external commands more sanelier
    #
    my ($cmd, $params, $failure_type) = @_;
    
    isCommandInPath($cmd, $failure_type);
    
    # join the parameters to the command
    my $param_str = join " ", map {formatParams($_)} @{$params};
    
    my $cmd_str = $cmd . " " . $param_str;
    
    logExternalCommand($cmd_str);

    # make sure that all went well
    if (system($cmd_str)) {
         handleCommandFailure($cmd_str, $failure_type)
    }
}

sub formatParams {

    #---------
    # Handles and formats the different ways of passing parameters to 
    # checkAndRunCommand
    #
    my $ref = shift;
    
    if (ref($ref) eq "ARRAY") {
        return join(" ", @{$ref});
    } elsif (ref($ref) eq "HASH") {
        return join(" ", map { $_ . " " . $ref->{$_}} keys %{$ref});
    }
    croak 'The elements of the $params argument in checkAndRunCommand can ' .
        'only contain references to arrays or hashes\n';
}


sub handleCommandFailure {
    #-----
    # What to do when all goes bad!
    #
    my ($cmd, $failure_type) = @_;
    if (defined($failure_type)) {
        if ($failure_type == DIE_ON_FAILURE) {
            croak "**ERROR: $0 : " . $! . "\n";
        } elsif ($failure_type == WARN_ON_FAILURE) {
            carp "**WARNING: $0 : " . $! . "\n";
        }
    }
}


######################################################################
# MISC

sub printAtStart {
print<<"EOF";
---------------------------------------------------------------- 
 $0
 Copyright (C) Michael Imelfort
    
 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
---------------------------------------------------------------- 
EOF
}

__DATA__

=head1 NAME

    parseCazyBlast

=head1 COPYRIGHT

   copyright (C) Michael Imelfort

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 DESCRIPTION

   Parse the output of a cazy based blast to determine best matching regions

=head1 SYNOPSIS

    parseCazyBlast -blast|b BLAST_FILE

      -blast -b BLAST_FILE         Blast file to parse
      -out -o OUT_FILE             File to write results to
      [-conserved_length -c]       The min amount (%) of thre hit which must map to the contig (default 0.80)               
      [-length -l INTEGER]         The length of the match that the HSP must be higher than
      [-percid -p INTEGER]         The required percent identity that the HSP must be higher than
      [-evalue -e INTEGER]         The required e-value that the HSP must be lower than
      [-bits -b INTEGER]           The required bits score that the HSP must be higher than
      [-mismatches -m INTEGER]     The number of mismatches that the HSP must be lower than
      [-format -f FORMAT]          Format of the blast file (default: standard blast output)
      [-help -h]                   Displays basic usage information
         
=cut

