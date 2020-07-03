#!/usr/bin/perl
# Nazar, 24/09/2011 

($file, $threshold) = @ARGV;

# Removing the existing matrices 
# *******************************************
system "del matrix.s";
system "rm matrix.s";

# Combine the interacting protein in one list
# *******************************************
open (File, "<$file");
while (<File>)
{
	if (/(\S+)\s+(\S+)/)
	{
		push @a1, "$1";
		push @a1, "$2";
	}
}
$n1 = @a1;
print "The numver of proteins is $n1\n"; 

# Remove duplicated proteins in the list
# **************************************
open (OUT, ">unique-pro.txt");
my %hash   = map { $_ => 1 } @a1;
my @unique = keys %hash;                

foreach my $uni (@unique)
{
	print OUT "$uni\n";
	push @a2, "$uni";
}
close OUT;

$n2 =  @a2;
print "The numver unique proteins is $n2\n";

# Retrieve the corresponding proteins
# ***********************************
$path = '';
$Proteins = 'unique-pro.txt';
$Sequences = 'all_seq.fasta';
$outfile = 'sequences.txt';
open(PROT, $path.$Proteins) || die ("Could not open file <br> $!");
open(SEQ, $path.$Sequences) || die ("Could not open file <br> $!");
open(OUT,">$path$outfile") || die ("Could not open file \n$!");

$line = <PROT>;
@sequence = <SEQ>;
$size = @sequence;

while ($line)
{
 	$p = $line;
 	chomp($p);
  $i=0; $flag = 0;
  while($i<$size && !$flag){
     if($sequence[$i] =~ />/){
     	 $sequence[$i] =~ />(\S+)\s/;
       if($1 eq $p){
       	 $flag = 1;
       	 $flag2 = 0;
       	 print OUT $sequence[$i];
         $i++;
       	 while($i<$size && !$flag2){
           if($sequence[$i] =~ />/){
           	 $flag2 = 1;
           }else {
           	print OUT $sequence[$i];
            }
           $i++;
         }
       }
     }
     $i++;
  }
  $line = <PROT>;
}
close(PROT);
close(SEQ);
close(OUT);  

# Replace the protein lables by integers
# **************************************
open (FILE,"<sequences.txt");
open (OUT,">data.s");

$i=1;
while (<FILE>)
{
	if (/>/)
	{
		s/.*/>$i /;
		print OUT $_;
		$i++;
		chomp  $_;
	}else{
		print OUT $_;
	}
}
close (FILE);
close (OUT);

# Returning one protein at a time
# *********************************
open IN,"<data.s";
$line = "0";
while(<IN> )
{
       chomp;
       open OUT,">testfile.s" if $line == 0;
       print OUT "$_\n";
       $line++;
       if($line == 2){ close OUT;

# Running FASTA
#**************
print ".";
system "fasta36.exe -q -s BL62 -f -11 -g -1 -d 0 -H -b $n2 testfile.s data.s > fastaoutput.s";

# Returning the sequence label and score lines
# ********************************************
@arrayone=();
open (File,"<fastaoutput.s");
while(<File>)
{
	if (/(\d.* \s+\(.*)/)
	{
		push @arrayone, "$1\n";
	}
}
close File;

# Sorting
# *******
@arraytwo=();
for (@arrayone)
{
	chomp; push @array, $_;
	if (@array == $n2) 
	{
		@AB = sort { $a <=> $b; } @array;
		foreach (@AB)
		{
			push @arraytwo, "$_\n";
		}
		@array = ();
	}
}

# Returning the scores only and push them to an array
# ***************************************************
@score=();
for (@arraytwo)
{
		/\d.*\s+\(*.*\)\s+(\S+)\s+(\S+)\s+(\S+)/;
	push @score,"$1";
}
close File;

# save the score of the substrings in columns
# *******************************************
@new = @score;
chomp @new;
$file = "matrix.s";
@prev = ();
open(fhOut, "$file");
while ($Line = <fhOut>)
{
	chomp($Line);
	push (@prev,$Line);
}
close fhOut;

open(fhOut, ">$file") or die "can't open $file";
$i=0;
if ($#prev <= 0)
{
	$spc = ",";
}else{
	$spc = ",";
}
foreach (@new)
{
	$Line = $prev[$i++];
	$Line = $Line.$spc.$_."\n";
	print fhOut $Line;
}
$line = "0";
}
}
print "\n";
close (fhOut);
close (IN);

# prepare the matrix file
#***************************
open(F,"<matrix.s");
open(F1,"<unique-pro.txt");
open(OUT,">similaritymatrix.txt");

@a = <F>;
chomp @a;
@b = <F1>;
chomp @b;

for $i(0..$n2-1)
{
	print OUT "$b[$i],$a[$i]\n";
}
close OUT;
close F;
close F1;
system "del *.s";
