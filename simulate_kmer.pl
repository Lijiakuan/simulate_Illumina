#!usr/bin/perl
# Program:
#       get kmer
# 2017/2/9      linchao

use strict;
use warnings;

die "perl $0 <kmer>, <outputfile>, <inputFileList>\n" if(@ARGV==0);
my($kmer, $outputfile,@inputfileList) = @ARGV;


die("inputFileList not input") if(!defined($inputfileList[0]));
die("kmer not input") if(!defined($kmer));
die("outputfile not input") if(!defined($outputfile));

# ========================= Main Program ============================

open(OUTPUT,">$outputfile") or die "can't open $outputfile:$!";

my %DNAtimes;
for(my $i=0;$i<=$#inputfileList;$i++)
{
    open(INPUT, "<$inputfileList[$i]") or die "can't open $inputfileList[$i]:$!";
    print "$inputfileList[$i]========================================\n";
    $/ = ">";
    while(<INPUT>)
    {
        chomp;
        next if(length $_ == 0);
        
        my @array = split(/\n/, $_);
        my $seq = join("", @array[1..$#array]);

        &k_frequen($seq, $kmer,\%DNAtimes);

    }
    close INPUT;
}

print OUTPUT "output_1:\n";

#delete the DNA or the oppotion DNA by the appearance times.
my @keys = keys(%DNAtimes);
for(my $i=0; $i<@keys; $i++)
{
    my $oppotion = &revcom($keys[$i]);
    if(defined($DNAtimes{$oppotion})&&defined($DNAtimes{$keys[$i]}))
    {
        if($DNAtimes{$keys[$i]} > $DNAtimes{$oppotion}){
            delete($DNAtimes{$keys[$i]});
        }else{
            delete($DNAtimes{$oppotion});
        }
    }
}

# output DNA and the times of the same DNA
@keys = keys(%DNAtimes);
for(my $i=0; $i<@keys; $i++)
{
    print OUTPUT "$keys[$i], $DNAtimes{$keys[$i]}\n";
}
my %kmer = &kmer(values(%DNAtimes));

# output the times of DNA,and the kind of DNA.1
@keys = keys(%kmer);
@keys = sort({$a <=> $b} @keys);
print OUTPUT "output_2:\n";
for(my $i=0; $i<@keys; $i++)
{
    print OUTPUT "$keys[$i], $kmer{$keys[$i]}\n";
}
close OUTPUT;


# ======================== Son Function ===========================
# get frequence that length is k
sub k_frequen
{
    my ($seq, $kmer,$modelHash) = @_;
    my $model;
    for(my $j=0; $j < length($seq)-$kmer; $j++)
    {
        $model = substr($seq, $j, $kmer);
        if(!defined($modelHash->{$model})){
            $modelHash->{$model} = 0;
        }
        $modelHash->{$model}++;
    }
}


# get kmer
sub kmer
{
    my @modelArr = @_;
    my %kmers;
    for(my $i=0; $i<@modelArr; $i++)
    {
        if(!defined($kmers{$modelArr[$i]})){
            $kmers{$modelArr[$i]} = 0;
        }
        $kmers{$modelArr[$i]}++;
    }
    return %kmers;
}

# oppositon DNA
sub revcom{
    my ($seq) = @_;
    $seq = reverse $seq;
    $seq =~ tr/[ACGTacgt]/[TGCAtgca]/;
    return $seq;
}