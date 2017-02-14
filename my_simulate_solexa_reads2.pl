#! usr/bin/perl
# Program:
#       simulate the program of reading DNA
# 2017/2/9      linchao

use strict;
use warnings;
use Getopt::Long;

=head1 Option
    -i      <string>        the road of DNA file
    -r      <int>           set the length of max reading
    -im     <int>           set the mean of insertsize distribution
    -id     <int>           set the sd of insertsize distribution
    -e      <float>         set the error rate
    -sr     <float>         set the SNP rate
    -ir     <float>         set the Indel rate
    -o      <string>        set the road of output file
=cut

my ($inputref, $read_len, $coverage, $insertsize_mean, $insertsize_sd, $error_rate);
my ($heterSNP_rate, $heterIndel_rate, $outputref);

GetOptions(
    "i:s"=>\$inputref,
    "r:i"=>\$read_len,
    "c:i"=>\$coverage,
    "im:i"=>\$insertsize_mean,
    "id:i"=>\$insertsize_sd,
    "e:f"=>\$error_rate,
    "sr:f"=>\$heterSNP_rate,
    "ir:f"=>\$heterIndel_rate,
    "o:s"=>\$outputref,
);

# ========================= set default value ============================
$read_len = 100 if(!defined($read_len));
$coverage = 40 if(!defined($coverage));
$insertsize_mean = 500 if(!defined($insertsize_mean));
$insertsize_sd = 25 if(!defined($insertsize_sd));
$error_rate = 0.01 if(!defined($error_rate));
$heterSNP_rate = 0 if(!defined($heterSNP_rate));
$heterIndel_rate = 0 if(!defined($heterIndel_rate));
$outputref = "read" if(!defined($outputref));
$inputref = "ref_sequence.fa" if(!defined($inputref));

die("arg contains the num that is no bigger than 0\n") if($read_len<=0||$insertsize_mean<=0||$insertsize_sd<=0||$error_rate<0||$heterSNP_rate<0||$heterIndel_rate<0);

# ========================== Main Program ================================
open (INPUT, "<$inputref") || die ("can't open file $inputref\n");
open (OUTPUT, ">$outputref\1.fa") || die ("can't open file $outputref\n");
open (OUTPUT2,">$outputref\2.fa") || die ("can't open file $outputref\n");
my $count = 0;
while(<INPUT>)
{
    chomp;
    # filter blank row
    next if(($_ eq ""));
    if(/^[agctAGCT]*$/)
    {
        my $fa_sqe = $_;
        my $mo_sqe = $_;
        # deal the fa_sqe
        if($heterSNP_rate > 0)
        {
            $fa_sqe = &simulate_SNP($fa_sqe, $heterSNP_rate);
        }
        if($heterIndel_rate > 0)
        {
            $fa_sqe = &simulate_Indel($fa_sqe, $heterIndel_rate);
        }
        &simulate_read($fa_sqe, $mo_sqe, $insertsize_mean, $coverage, $insertsize_sd, $error_rate, $read_len);
    }
}
close(OUTPUT);
close(OUTPUT2);
close(INPUT);


# ============================= Son Function =============================
# read DNA
sub simulate_read{
    my ($sqe, $coverage, $insertsize_mean, $insertsize_sd, $error_rate, $read_len)
    my $len = length($sqe);
    # preview the num after breaking;
    my $read_num = $len / $insertsize_mean * $coverage;
    # get the length of every breaking DNA
    my %Isize = &simulate_insertsize($insertsize_sd, $insertsize_mean, $read_num);
    # get the position of every reading error
    my %Esize = &simulate_err($error_rate, $read_len, $read_num, $coverage);
    # start read the DNA
    &startReads($refsqe, $read_len, \%Isize, \%Esize);
}
# get the length of every breaking DNA
sub simulate_insertsize{
    my ($sd, $mean, $num) = @_;
    my $pi = 3.1415926;
    my (%Isize,$total, $total1);
    $total = 0;
    for(my $i = $mean-5*$sd; $i < $mean+5*$sd; $i++)
    {
        $Isize{$i} = 1/$sd/sqrt(2*$pi)/exp((($i-$mean)/$sd)**2);
        $total += $Isize{$i};
    }
    $total1 = 0;
    for(my $i = $mean-5*$sd; $i < $mean+5*$sd; $i++)
    {
        $Isize{$i} = ($Isize{$i}/$total)*$num;
        $total1 += $Isize{$i};
    }
    if($total1 < $num)
    {
        $Isize{$mean} = $num - $total1;
    }
    return %Isize;
}

# get the position of the read Error;
sub simulate_err{
    my ($rate, $len, $num) = @_;
    my (%Esize, $sizes);
    for(my $i = 0; $i < $len; $i++)
    {
        $sizes="";
        my $enum = int($rate * $i *10 * $num);
        # print("$enum -=-=-=-=-=-=- $num\n");
        while($enum > 0)
        {
            #print("$sizes-----------------\n");
            my $size = int(rand($num));
            if(!($sizes =~ m/<$size>/))
            {
                #print("x-----------------\n");
                $Esize{$size} = $i;
                $enum--;
                $sizes .= "<$size>";
            }
        }
    }
    return %Esize;
}

# simulate the SNP problem
sub simulate_SNP{
    my ($sqe, $rate) = @_;
    my $length = length($sqe);
    my $num = $rate * $length;
    my %table = (
        "A" => ["G", "C", "T"],
        "G" => ["A", "C", "T"],
        "C" => ["G", "A", "T"],
        "T" => ["G", "C", "A"]
    );
    my $sizes="";
    while($num > 0)
    {
        my $size = rand($length);
        if(!($sizes =~ m/<$size>/))
        {
            my $n = int(rand(3));
            # 改进用正则
            $sqe =~ s/(.{$size})([agctAGCT])/$1$table{$2}[$n]/;
            $num--;
            $sizes .= "<$size>";
        }
    }
    return $sqe;
}

# simulate the Indel problem
sub simulate_Indel{
    my ($rate, $sqe) = @_;
    my $length = length($sqe);
    my $num = int($rate * $length / 2);
    my $i = $num;
    my $sizes="";
    # delete
    while($i > 0)
    {
        my $size = rand($length);
        if(!($sizes =~ m/<$size>/))
        {
            my $deletelen = int(rand(3));
            $sqe =~ s/(.{$size})[agctAGCT]{$deletelen}/$1/;
            $i--;
            $sizes .= "<$size>";
        }
    }
    # add 
    $sizes = "";
    my @table = ("A", "G", "C", "T");
    while($i > 0)
    {
        my $size = rand($length);
        if(!($sizes =~ m/<$size>/))
        {
            my $n = int(rand(3));
            my $addDNA = "";
            while($n>0)
            {
                $addDNA .= $table[int(rand(4))];
                $n--;

            }
            $sqe =~ s/(.{$size})/$1$addDNA/;
            $i--;
            $sizes .= "<$size>";
        }
    }
    return $sqe;
}

# start reads the DNA, accoding the data
sub startReads{
    my ($sqe, $sqe2, $read_len, $Isize, $Esize) = @_;
    my $length = length($sqe);
    my $length2 = length($sqe2);
    my @keys = keys(%{$Isize});
    my state $n=1;
    for(my $i=0; $i<@keys; $i++)
    {
        for(my $j=0; $j < $Isize->{$keys[$i]}; $j++)
        {
            if(int rand 2 ==1 ){
                my $pos = int(rand($length-$keys[$i]));
                # breaking DNA
                my $origin = substr($sqe, $pos, $keys[$i]);
            }else{
                my $pos = int(rand($length2-$keys[$i]));
                # breaking DNA
                my $origin = substr($sqe2, $pos, $keys[$i]);
            }
            # opposion DNA
            my $copy = $origin;
            $copy =~ tr/ACGTacgt/TGCAtgca/;
            # start
            # a reads the origin
            my $a = substr($origin, 0, $read_len);
            # b reads the opposion
            my $b = substr($copy, -$read_len, $read_len);

            if($length > 1000)
            {
                $a = reverse($a);
                $a =~ tr/ACGTacgt/TGCAtgca/;
            }
            else
            {
                $b = reverse($b);
                $b =~ tr/ACGTacgt/TGCAtgca/;
            }
            # output
            my $d=int ((rand 10)%2);
			if ($d) {
					print(OUTPUT ">reads_$n\_100 1\n$a\n");
                    print(OUTPUT2 ">reads_$n\_100 2\n$b\n");
			}else{
					print(OUTPUT ">reads_$n\_100 1\n$b\n");
                    print(OUTPUT2 ">reads_$n\_100 2\n$a\n");
			}
            $n++;
        }
    }
}