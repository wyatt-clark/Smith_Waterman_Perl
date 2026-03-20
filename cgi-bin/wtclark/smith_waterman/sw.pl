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
			$s1 =~ s/\W//g;
			$s2 =~ s/\W//g;
		
		}
		else{
			$s1 =~ s/\s//g;
			$s2 =~ s/\s//g;
		}
		@s_1 = split '', $s1;
		@s_2 = split '', $s2;
		
	}

	(my $m1_ref, my $m2_ref, my $best_ref, my $score) = score_matrix(\@s_1, \@s_2);
	my @m1 = @$m1_ref;
	my @m2 = @$m2_ref;
	my @best = @$best_ref;

	
	#(my $t_1, my $t_2) = trace_back(\@m1, \@m2, \@best, \@s_1, \@s_2, $score);
	trace_back(\@m1, \@m2, \@best, \@s_1, \@s_2, $score);

	#my @t1 = @$t_1;
	#my @t2 = @$t_2;
	
	#if($output_type eq "alignment (tabb delineated text)"){
	#	print_results_text(\@t1, \@t2, \@s_1, \@s_2);
	#}
	#else{
	#	print_results_html(\@t1, \@t2, \@s_1, \@s_2);
	#}
		
	
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
			
		}
	}
	return \@m1, \@m2, \@best_loc, $best;

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
	(my $ref_1, my $ref_2, my $best_ref, my $s1_ref, my $s2_ref, my $score) = @_;
	my @m1 = @$ref_1;
	my @m2 = @$ref_2;
	my @best = @$best_ref;
	my @s2 = @$s2_ref;
	my @s1 = @$s1_ref;
	my @ts1;
	my @ts2;
	my $start1;
	my $start2;
	
	my @ret_1;
	my @ret_2;
	my @middle;
	
	my $i = $best[0];
	my $j = $best[1];
	
	my $match_count = 0;
	my $stop = 1;
	my $gap = 0;
	my $miss = 0;
	while($stop > 0){
		$stop = $m1[$i][$j];
		my $arrow = $m2[$i][$j];
		
		if($arrow == 2){
			unshift @ret_1, $s1[$i];
			unshift @ret_2, $s2[$j];
			if($s1[$i] eq $s2[$j]){
				unshift @middle, "\|";
				$match_count++;
			}
			else{
				unshift @middle, "\*";
				$miss++;
			}
			unshift @ts1, $i;
			unshift @ts2, $j;
			$i--;
			$j--;
		}
		if($arrow == 1){
			unshift @ret_1, " ";
			unshift @ret_2, $s2[$j];
			unshift @middle, " ";
			unshift @ts1, -1;
			unshift @ts2, $j;
			$j--;
			$gap++;
			
		}
		if($arrow == 3){
			unshift @ret_1, $s1[$i];
			unshift @ret_2, " ";
		
			unshift @middle, " ";
			unshift @ts1, $i;
			unshift @ts2, -1;
			$i--;
			$gap++;
		}
		
		if($i < 0 || $j< 0 || $m1[$i][$j] == 0){
			$stop = -5;
			
		}
		$start1 = $i + 2;
		$start2 = $j + 2;
	
	}
	my $pid = sprintf("%.2f", ($match_count / ($#middle + 1)) * 100);
	my $stop1 = $best[0]+1;
	my $stop2 = $best[1]+1;
	# if($output_type ne "alignment (tabb delineated text)"){
	# 		print table(
	# 			{-border=>''},
	# 			Tr(td[("start"),("stop")]),
	# 			Tr(td[($start1),($stop1)]),
	# 			Tr(td[($start2),($stop2)]),
	# 			Tr(td[("score"),($score)]),
	# 			Tr(td[("percent id"),($pid)]),
	# 			Tr(td[("gaps"),($gap)]),
	# 			Tr(td[("miss-matches"),($miss)])
			
	# 	);
	# 			print "<br><br>";
	
	# 	print table(
	# 			{-border=>''},
	# 			TR({-align=>'center'},
	# 			[
	# 				td(\@ret_1),
	# 				td(\@middle),
	# 				td(\@ret_2)
	# 			]
	# 		)
	# 	);
	# }
	# else{
		print map "$_\t", @ret_1;
		print "\n";
		print map "$_\t", @middle;
		print "\n";
		print map "$_\t", @ret_2;
		
	
	#}
	#return \@ts1, \@ts2;

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

