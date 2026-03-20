#!/usr/bin/perl -w
use strict;
use List::Util;

#for now i match is 1, gap is -1 and gapopen is 0
my $match;
my $gap_open;
my $gap_ext;
my $missmatch;
my $matrix;
my %sub_mat;
my $type;
my $output_type;
#type defined in sw_wrapper.pl
#my $type = "protein";
sub initialize_sw{
	$type = shift;
	$output_type = shift;
	$gap_open = shift;
	$gap_ext = shift;
	$matrix = shift;
	$match = shift;
	$missmatch = shift;

	if($type eq "protein"){%sub_mat = load_scoring_matrix("$matrix")};

}


sub sw{
	my $s1 = shift;
	my $s2 = shift;
	my @s_1;
	my @s_2;
	$s1 =~ s/\s+/ /g;
	$s2 =~ s/\s+/ /g;
	#print "$gap_open<br>$gap_ext<br>";
	#print "type is $type\n";
	if($type eq "Text split on word"){
		
		@s_1 = split ' ', $s1;
		@s_2 = split ' ', $s2;
		
	}
	
	else{
		if($type eq "text split on letter"){
			$s1 =~ s/\s/ /g;
			$s2 =~ s/\s/ /g;
		
		}
		else{
			$s1 =~ s/\s//g;
			$s2 =~ s/\s//g;
		}
		@s_1 = split '', $s1;
		@s_2 = split '', $s2;
		
	}

	score_matrix(\@s_1, \@s_2);
	#my @m1 = @$m1_ref;
	#my @m2 = @$m2_ref;
	#my @best = @$best_ref;

	
	#(my $t_1, my $t_2) = trace_back(\@m1, \@m2, \@best);
	#my @t1 = @$t_1;
	#my @t2 = @$t_2;
	
	#if($output_type eq "alignment (tabb delineated text)"){
	#	print_results_text(\@t1, \@t2, \@s_1, \@s_2);
	#}
	#else{
	#	print_results_html(\@t1, \@t2, \@s_1, \@s_2);
	#}
	#open(FILE, ">alignment");
	
	
	#close FILE;	
	
}1;

sub print_out{
	(my $m1_ref, my $m2_ref, my $best_ref, my $s1_ref, my $s2_ref, my $go_to_ref) = @_;
	my @m1 = @$m1_ref;
	my @m2 = @$m2_ref;
	my @best = @$best_ref;
	my @s_2 = @$s2_ref;
	my @s_1 = @$s1_ref;
	my @go_to = @$go_to_ref;

	
	(my $t_1, my $t_2) = trace_back(\@m1, \@m2, \@best, \@go_to);
	my @t1 = @$t_1;
	my @t2 = @$t_2;
	
	if($output_type eq "alignment (tabb delineated text)"){
		print_results_text(\@t1, \@t2, \@s_1, \@s_2);
	}
	else{
		print_results_html(\@t1, \@t2, \@s_1, \@s_2);
	}


}1;
sub score_matrix{
	my $s_1 = shift;
	my $s_2 = shift;
	my @s1 = @$s_1;
	my @s2 = @$s_2;
	my $best = 0;
	my @best_loc;
	my @m1;
	my @m2;
	my @go_to;
	
	#we first go though the first column and row
	for my $i (0..$#s1){
		if($s1[$i] eq $s2[0]){
			$m1[$i][0] = $match;
			#arrow might have to be changed
			$m2[$i][0] = 2;
		}
		else{
			#arrow might have to be changed
			$m1[$i][0] = 0;
			$m2[$i][0] = 3
		}
	}
	
	for my $i (0..$#s2){
		if($s2[$i] eq $s1[0]){
			$m1[0][$i] = $match;
			#arrow might have to be changed
			$m2[0][$i] = 2;
		}
		else{
			#arrow might have to be changed
			$m1[0][$i] = 0;
			$m2[0][$i] = 1
		}
	}
	
	#now go through rest of matrix
	#sequence 1 is on dim 1, seq 2 is in dim 2
	for my $i (1..$#s1){
		for my $j(1..$#s2){
			#print "i = $i j = $j\n";
			
			
			my $sub_val = $sub_mat{$s1[$i]}{$s2[$j]};
			my @scores;
			$scores[0] = score_up($m1[$i][$j-1], $m2[$i][$j-1]);
			if($type eq "protein"){
				$scores[1] = score_diag($m1[$i-1][$j-1], $m2[$i-1][$j-1], $sub_val);
			}
			else{
				my $m = 0;
				if($s1[$i] eq $s2[$j]){$m = 1};
				$scores[1] = score_diag_normal($m, $m1[$i-1][$j-1], $m2[$i-1][$j-1]);
			
			}
			$scores[2] = score_left($m1[$i-1][$j], $m2[$i-1][$j]);
			
			#sort in descending
			my @sorted = sort{$scores[$b] <=> $scores[$a]} 0..$#scores;
			#print "best is $sorted[0] $scores[$sorted[0]]\n";
			$m1[$i][$j] = $scores[$sorted[0]];
			$m2[$i][$j] = $sorted[0] +1;
			
			if($scores[$sorted[0]] > $best){
				$best = $scores[$sorted[0]];
				($best_loc[0], $best_loc[1]) = ($i, $j);
			
			}
			if($scores[$sorted[0]] == 0 && $best > 0){
				if($best_loc[1] - $go_to[1] > 4){
					if($best_loc[0] - $go_to[0] > 4){
#			(my $m1_ref, my $m2_ref, my $best_ref, my $s1_ref, my $s2_ref, my $go_to_ref)
						print_out(\@m1, \@m2, \@best_loc, \@s1, \@s2, \@go_to);
						$go_to[0] = $i;
						$go_to[1] = $j;
						$best = 0;
					}
				}
			
			}
			
		}
	}
	#return \@m1, \@m2, \@best_loc;

}1;

sub score_diag{
	(my $score, my $dir_val, my $sub_val) = @_;
	#if match then returned score is 
	#print "$m\t$score\t$dir_val\n";
	my $n_score = $score + $sub_val;
	
	if($n_score < 0){$n_score = 0};
	return $n_score;


}1;

sub score_diag_normal{
	(my $m, my $score, my $dir_val) = @_;
	#if match then returned score is 
	#print "$m\t$score\t$dir_val\n";
	my $n_score;
	if($m == 1){
		$n_score = $score + $match;
	}
	else{
		$n_score = $score + $missmatch;
	}
	if($n_score < 0){$n_score = 0};
	return $n_score;


}1;

sub score_up{
	(my $score, my $dir_val) = @_;
	#print "$m\t$score\t$dir_val\n";
	my $n_score; 
	if($dir_val !=2){
		$n_score = $score + $gap_ext;
	}
	else{
		$n_score = $score + $gap_open;
	}
	if($n_score < 0){$n_score = 0};
	return $n_score;
	
}1;

sub score_left{
	(my $score, my $dir_val) = @_;
	#print "$m\t$score\t$dir_val\n";
	my $n_score; 
	if($dir_val !=2){
		$n_score = $score + $gap_ext;
	}
	else{
		$n_score = $score + $gap_open;
	}
	if($n_score < 0){$n_score = 0};
	return $n_score;
	
}1;

sub trace_back{
	(my $ref_1, my $ref_2, my $best_ref, my $go_to_ref) = @_;
	my @m1 = @$ref_1;
	my @m2 = @$ref_2;
	my @best = @$best_ref;
	my @go_to = @$go_to_ref;
	my @ts1;
	my @ts2;
	my $i = $best[0];
	my $j = $best[1];
	my $stop = 1;
	while($stop > 0){
		$stop = $m1[$i][$j];
		my $arrow = $m2[$i][$j];
		
		if($arrow == 2){
			
			unshift @ts1, $i;
			unshift @ts2, $j;
			$i--;
			$j--;
		}
		if($arrow == 1){
			
			unshift @ts1, -1;
			unshift @ts2, $j;
			$j--;
			
		}
		if($arrow == 3){
			
			unshift @ts1, $i;
			unshift @ts2, -1;
			$i--;
		}
		
		if($i < $go_to[0] || $j< $go_to[1]){
			$stop = -5;
		}
	
		
	}
	
	return \@ts1, \@ts2;

}1;

sub load_scoring_matrix{
	my $file = shift;
	open(FILE, "$file");
	my @data = <FILE>;
	close FILE;
	
	chomp(@data);
	my %hash;
	foreach my $line(@data){
		my @crap = split ' ', $line;
		$hash{$crap[0]}{$crap[1]} = $crap[2];
	
	}
	return %hash;
}1;

sub print_results_html{
	my $t_1 = shift;
	my $t_2 = shift;
	my $s1 = shift;
	my $s2 = shift;
	my @s_1 = @$s1;
	my @s_2 = @$s2;
	my @t1 = @$t_1;
	my @t2 = @$t_2;
	my @out1;
	my @out2;
	
	for my $p (@t1){
		if($p >= 0){
			push @out1, $s_1[$p];
		}
		else{
			push @out1, " ";
		}
	}
	
	
	for my $p (@t2){
		if($p >= 0){
			push @out2, $s_2[$p];
		}
		else{
			push @out2, " ";
		}
	}
	
	my @middle;
	for my $i (0..$#out1){
		if($out1[$i] eq $out2[$i]){
			$middle[$i] = "\|";
		}
		else{
			if($out1[$i] eq " " || $out2[$i] eq " "){
				$middle[$i] = " ";
			}
			else{
				$middle[$i] = "\*";
			}
		}
	
	
	}
	if($out1[0] ne $out2[0]){
		shift(@out1);
		shift(@out2);
		shift(@middle);
	
	}
	print table({-border=>''},
	TR({-align=>'center'},
	[
	td(\@out1),
	td(\@middle),
	td(\@out2)
	]
	)
	);
	print "<br>";
	print "<br>";

};

sub print_results_text{
	my $t_1 = shift;
	my $t_2 = shift;
	my $s1 = shift;
	my $s2 = shift;
	my @s_1 = @$s1;
	my @s_2 = @$s2;
	my @t1 = @$t_1;
	my @t2 = @$t_2;
	for my $p (@t1){
		if($p >= 0){
			print "$s_1[$p]\t";
		}
		else{
			print "\t";
		}
	}
	
	print "\n";
	for my $p (@t2){
		if($p >= 0){
			print "$s_2[$p]\t";
		}
		else{
			print "\t";
		}
	}
	print "\n";


};