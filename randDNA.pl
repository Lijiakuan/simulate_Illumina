#!/usr/bin/perl
# Program:
#       用于产生虚拟的DNA序列
# 2017/2/7      林超
use strict;
use warnings;

die "perl $0 <len> <outfile>\n" if(@ARGV==0);

my($len, $outfile) = @ARGV;


# ========================= Main Program =========================
open FILE, ">$outfile";

my @DNA;
my @Nucleotide = ("A", "G", "C", "T");
for(my $i = 0;$i <$len; $i++)
{	
	my $rand = int(rand(4));
    $DNA[$i] = $Nucleotide[$rand];

}
my $str = join("", @DNA);

print FILE ">test\n";
print FILE "$str\n";
close FILE;
