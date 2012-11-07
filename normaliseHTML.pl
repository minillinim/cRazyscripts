#!/usr/bin/perl
use strict;
use warnings;

use Bio::DB::GenBank;
use Bio::SeqIO::genbank;

my ($infn) = @ARGV;

my $domain = "UNSET";
my $p_name = "";
my $EC = "";
my $org = "";
my $in_line = 0;
my %gbs = ();

my $html_head = "<html><body>\n";
my $html_tail = "</body></html>\n";

open my $infh, "<", $infn or die $!;
my $meta_name = $infn;
$meta_name =~ s/_.*//;
`mkdir -p ../processed/META`;
open my $meta_fh, ">", "../processed/META/$meta_name.META.html" or die $!;
open my $chomped_fh, ">", "../processed/$meta_name.normalised.csv" or die $!;
print $chomped_fh join "\t", ("Domain", "Protein", "Organism", "EC", "GB");
print $chomped_fh "\n";
my $in_meta = 0;
my $in_seqs = 0;
while(<$infh>)
{
    if(0 == $in_meta)
    {
        next if not $_ =~ /<div><h1 style="color:#00688B;">/;
        $in_meta = 1;
        print $meta_fh $html_head;
        print $meta_fh $_;
    }
    elsif(1 == $in_seqs)
    {
        if($_ =~ /^   <tr class="royaume">/)
        {
            # <tr class="royaume"><td colspan="7" class="separateur1"><span class="tabulation">Archaea</span></td></tr><tr id="line_titre"><td>Protein Name</td> <td>EC#</td><td>Organism</td><td>GenBank</td><td>Uniprot</td><td>PDB/3D</td></tr>
            chomp $_;
            my @fields = split /span/, $_;
            # class="tabulation">Archaea</
            $fields[2] =~ s/.*>//;
            $fields[2] =~ s/<.*//;
            # Archaea
            $domain = $fields[2];
        }
        elsif($domain eq "UNSET")
        {
            # no point in jumping the gun
            next;
        }
        else
        {
            # we have a class set and we have a standard line
            chomp $_;
            if($_ =~ /^    <tr valign="top" onmouseover="this.bgColor='#F8FFD5';" /)
            {
                $in_line = 2;
                $_ =~ s/.*&nbsp;//;
                $_ =~ s/<\/t.*//;
                $p_name =  $_;
                $p_name =~ s/'//g;
            }
            elsif($in_line > 0)
            {
                if(!($_ =~ /<\/td/))
                {
                    <$infh>;
                }
                
                if($in_line == 2)
                {
                    #EC
                    $_ =~ s/.*link">//;
                    $_ =~ s/<.*//;
                    $_ =~ s/ //g;
                    $EC = $_;
                    $EC =~ s/'//g;
                }
                elsif($in_line == 3)
                {
                    #ORG
                    $_ =~ s/^   //;
                    $_ =~ s/<.*//;
                    $_ =~ s/^ $//;
                    $org = $_;
                    $org =~ s/'//g;
                }
                elsif($in_line == 4)
                {
                    #GENBANK
                    my @fields = split / /, $_;
                    foreach my $field (@fields)
                    {
                        next if not $field =~ /protein/;
                        $field =~ s/.*val=//;
                        chop $field;
                        $field =~ s/'//g;
                        $gbs{$field} = 1;
                        print $chomped_fh '';
                        print $chomped_fh join "'\t'", ($domain, $p_name, $org, $EC, $field);
                        print $chomped_fh "'\n";
                    }
                }
                elsif($_ =~ /^    <\/tr/)
                {
                    $in_line = 0;
                    $p_name = "";
                    $EC = "";
                    $org = "";
                }
                $in_line ++;
            }
        }
    }
    elsif(1 == $in_meta)
    {
        print $meta_fh $_;
        if($_ =~ /<\/table>/)
        {
            # out of the meta
            $in_seqs = 1;
            print $meta_fh $html_tail;
            close $meta_fh;
        }
    }
}

close $infh;
close $chomped_fh;
foreach my $acc (keys %gbs)
{
    my $gb = Bio::DB::GenBank->new();
    print "Getting $acc.gb\n";
    my $seq = $gb->get_Seq_by_id($acc);
    my $seqout = Bio::SeqIO->new(
                             -file   => ">../GB/$acc.gb",
                             -format => 'genbank',
                             );
    $seqout->write_seq($seq);
    print "written $acc.gb\n";
}

