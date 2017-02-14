#! usr/bin/perl

use strict;
use feature "state";

die "perl $0 <inputfile> <SNP_rate> <Indel_rate> <insertSize_mean> <insertSize_sd> <error_rate> <read_len> <coverage>\n" if(@ARGV==0);

my($inputfile, $SNP_rate, $Indel_rate, $insertSize_mean, $insertSize_sd, $error_rate, $read_len, $coverage) = @ARGV;

# ========================= set default value ============================
$read_len = 100 if(!defined($read_len));
$insertSize_mean = 500 if(!defined($insertSize_mean));
$insertSize_sd = 25 if(!defined($insertSize_sd));
$error_rate = 0.00 if(!defined($error_rate));
$SNP_rate = 0 if(!defined($SNP_rate));
$Indel_rate = 0 if(!defined($Indel_rate));
$coverage = 10 if(!defined($coverage));
$inputfile = "ref_sequence.fa" if(!defined($inputfile));

# ###################### Main Program #########################
open INPUT, "<$inputfile";
open READ1, ">read1.fa";
open READ2, ">read2.fa";

$/ = ">";
while(<INPUT>) {
    next if(length $_ == 0);
    
    my @array = split(/\n/, $_);
    my $seq = join("", @array[1..$#array]);

    my $ref1 = $seq;
    my $ref2 = &revcom($seq);
    my $length = length($ref1);
    if($SNP_rate>0 || $Indel_rate>0){
        $ref2 = &SNP_Indel($ref1, $SNP_rate, $Indel_rate);
    }
    my $PCR;
    for(my $i=0; $i<$length*$coverage; $i += 2*$read_len)
    {
        if(int rand 2 == 1) {
            $PCR = &cutSeq($ref1, $insertSize_mean, $insertSize_sd);
        }
        else {
            $PCR = &cutSeq($ref2, $insertSize_mean, $insertSize_sd);
            #$PCR = &revcom($PCR);
        }
        &getReads($PCR, $error_rate, $read_len);
    }
}

close INPUT;
close READ1;
close READ2;

# ###################### Son Function #########################
# simulate the SNP and Indel problom
sub SNP_Indel{
    my ($ref, $sRate, $iRate) = @_;
    my $len = length $ref;
    my $sNum = int $sRate * $len;
    my $iNum = int $iRate * $len;
    my @array = ("A", "G", "C", "T");
    for(my $i=0; $i<$sNum; $i++) {
        substr($ref, int rand $sNum, 1) = $array[int rand 4];
    }
    for(my $i=0; $i<$iNum; $i++) {
        if(int rand 2 == 1) {
            substr($ref, int rand $iNum, int rand 3) = "";
        }
        else {
            my $position = int rand $iNum;
            for(my $i=0; $i<rand 3;$i++) {
                substr($ref, $position , 1) = @array[int rand 4];
            }
        }
    }
    return $ref;
}

# simulate the option of PCR
sub cutSeq{
    my ($ref, $mean, $range) = @_;
    my $len = int($mean + rand 2 * $range - $range);
    my $PCR = substr $ref, int rand(length($ref) -$len), $len;
    return $PCR;
}



# start read the DNA
sub getReads{
    my ($PCR, $rate, $read_len) = @_;
    my @array = ("A", "G", "C", "T");
    state $n = 0;
    my $len = length $PCR;
    my $read1 = substr($PCR, 0, $read_len);
    my $read2 = substr($PCR, -$read_len, $read_len);
    if($len > 1000){
        $read1 = &revcom($read1);
    }else{
        $read2 = &revcom($read2);
    }
    for(my $i=0; $i<$read_len; $i++) {
        if(rand(1) < $rate)
        {
            substr($read1, $i, 1) = @array[int rand 4];
        }
        if(rand(1)<$rate)
        {
            substr($read2, $i, 1) = @array[int rand 4];
        }
    }
    print READ1 ">reads_$n\_100 1\n$read1\n";
    print READ2 ">reads_$n\_100 2\n$read2\n";
    $n++;
}

# oppositon DNA
sub revcom{
    my ($seq) = @_;
    $seq = reverse $seq;
    $seq =~ tr/[ACGTacgt]/[TGCAtgca]/;
    return $seq;
}