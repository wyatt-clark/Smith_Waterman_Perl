#!/usr/bin/perl -w
use strict;

#this will have to be the absolute path on your machine
require "/Users/wtclark/Desktop/RNA/www/cgi-bin/wtclark/smith_waterman/sw.pl";
require "/Users/wtclark/Desktop/RNA/www/cgi-bin/wtclark/smith_waterman/is_protein.pl";



#sequence 1
my $seq1 = "four score and seven years ago our fathers brought forth on this continent a new nation";
#sequence 2
my $seq2 = "the score was seven goals go fathers forth inning grand slam nation";
#type should be 
  #<OPTION> protein
   # <OPTION> DNA
    #<OPTION> Text split on word
    #<OPTION> text split on letter
my $type = "Text split on word";
    #<OPTION> alignment (html table)
    #<OPTION> alignment (tabb delineated text)
my $output_form = "alignment (tabb delineated text)";
#-2
my $gapop = 0;
#-1
my $gapext = -1;
# "pam250_easy" or "blosum62_easy" 
my $matrix = "pam250_easy";
#2
my $match = 2;
#-
my $missmatch = -1;



if($type eq "DNA" || $type eq "protein"){
	my @seq1_array = split "\n", $seq1;

	if($seq1_array[0] =~ /^>/){
		shift @seq1_array;
	}
	$seq1 = join '', @seq1_array;
	$seq1=~ s/\W//g;
	$seq1=~ s/\d//g;

	my @seq2_array = split "\n", $seq2;
	if($seq2_array[0] =~ /^>/){
		shift @seq2_array;
	}
	$seq2 = join '', @seq2_array;

	$seq2=~ s/\W//g;
	$seq2 =~ s/\d//g;

	#print "$seq1<br>$seq2<br>";

}
#$seq1=~ s/\W//g;

#$seq2=~ s/\W//g;
if($output_form ne "alignment (tabb delineated text)"){
	print start_html("Smith Waterman alignment");
}
if($type eq "protein"){
	my $a = is_protein($seq1);
	my $b = is_protein($seq2);
	#print "$a $b<br>";
	if(is_protein($seq1) == 0 || is_protein($seq2) == 0){
		print "It seems that you have illegal characters in your sequence<br>";	
		exit;
	}
}
initialize_sw($type, $output_form, $gapop, $gapext, $matrix, $match, $missmatch);

#print "$seq1\n";
#print "<br>\n";
#print "$seq2\n";
#print "<br>\n";

sw($seq1, $seq2);


