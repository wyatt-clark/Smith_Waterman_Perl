#!/usr/bin/perl -w
use strict;

sub is_protein{
	my $string = shift;
	my $ans = 0;
	if($string !~ /[^ARNDCQEGHILKMFPSTWYVBZX\*]/){
		$ans = 1;
	}


	return $ans;
}1;
