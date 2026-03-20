#!/usr/bin/perl -w
use strict;


use CGI qw(:standard);
use CGI::Carp qw(warningsToBrowser fatalsToBrowser);
use CGI::Pretty qw( :html3 );
require "sw.pl";
require "is_protein.pl";


print header;



my %form;
foreach my $p(param()){
        $form{$p} = param($p);

}


my $seq1 = uc($form{"seq1"});
my $seq2 = uc($form{"seq2"});
my $type = $form{"data_type"};
my $output_form = $form{"output_type"};
my $gapop = $form{"gap_open"};
my $gapext = $form{"gap_ext"};
my $matrix = $form{"matrix"};
my $match = $form{"match"};
my $missmatch = $form{"missmatch"};

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


