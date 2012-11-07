#!/usr/bin/env perl
###############################################################################
#
#    cmpparsed.pl
#    
#    Compare parsed CAZY blast results for shared components
#
#    Remove overlapped regions from multiple parsed cazy results 
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
use List::Util qw(max);

#CPAN modules

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
my @files = split ",", $global_options->{'files'};
if($#files < 1) { croak "**ERROR: You need AT LEAST 2 files\n"; }
# 
# We assume that the files are sorted by contigId and then by position
#
my %contig_2_pos  = ();
my %file_2_reduced = ();
my $file_counter = -1;
foreach my $file (@files)
{
    $file_counter ++;
    my @tmp_reduced = ();
    $file_2_reduced{$file_counter} = \@tmp_reduced;
    my $fh = openRead($file);
    while(<$fh>)
    {
        chomp $_;
        my @fields = split(/\t/,$_);
        if(!exists $contig_2_pos{$fields[0]}) { my @tmp = (); $contig_2_pos{$fields[0]} = \@tmp; }
        $_ .= "\t$file_counter";
        push @{$contig_2_pos{$fields[0]}}, $_;
    }
    close $fh;
}

foreach my $cid (keys %contig_2_pos)
{
    my @sorted = sort sortHash @{$contig_2_pos{$cid}};  
    my $current_end = -1;
    my $current_best_region = "unset";
    my $current_best_file = 0;
    foreach my $line (@sorted)
    {
        my @fields = split(/\t/, $line);
        if("unset" eq $current_best_region)
        {
            $current_best_region = $line;
            $current_end = $fields[3];
        }
        else
        {
            if($fields[2] < $current_end)
            {
                # pileup
                my @tmp_fields = split(/\t/, $current_best_region);
                if($fields[5] < $tmp_fields[5])
                {
                    $current_best_region = $line;
                    $current_best_file = $fields[-1];                   
                } 
                $current_end = max($fields[3], $tmp_fields[3]);
            }
            else
            {
                # ok
                # put away the previous region
                add2Best($current_best_file, $current_best_region);
                print ${$file_2_reduced{$current_best_file}}[-1]."\n";
        
                $current_best_region = $line;
                $current_best_file = $fields[-1];                
                $current_end = $fields[3];
            }
        }
    }
    # put the last one on
    if("unset" ne $current_best_region)
    {
        my @tmp_fields = split(/\t/, $current_best_region);
        add2Best($current_best_file, $current_best_region);
        print ${$file_2_reduced{$current_best_file}}[-1]."\n";
    }
}
# print the lot to file
foreach my $file_counter (keys %file_2_reduced)
{
    my $fn = $files[$file_counter].".reduced";
    my $fh = openWrite($fn);
    foreach my $line (@{$file_2_reduced{$file_counter}})
    {
        print $fh "$line\n";
    }
    close $fh;
}

######################################################################
# CUSTOM SUBS
######################################################################
sub add2Best
{
    my ($fc, $line) = @_;
    my @fields = split(/\t/, $line);
    pop @fields;
    $line = join("\t", @fields);
    push  @{$file_2_reduced{$fc}}, $line;    
}

sub sortHash
{
    my @a_fields = split("\t", $a);
    my @b_fields = split("\t", $b);
    return int($a_fields[2]) <=> int($b_fields[2]);
}
######################################################################
# TEMPLATE SUBS

######################################################################
# PARAMETERS

sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "files|f:s");
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
    if(!exists $options{'files'} ) { printParamError ("Need files to compare"); }
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

    cmpparsed.pl

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

   Insert detailed description here

=head1 SYNOPSIS

    cmpparsed.pl -files|f FILE1,FILE2[,FILE3,..]

      -files -f FILE1,FILE2[,FILE3,..]   File to parse and compare
      [-help -h]                         Displays basic usage information

      Parse N files and output 2 matricies:
      1 - (Nx1) Number unique to file N
      2 - (NxN) Pairwise shared 
         
=cut

