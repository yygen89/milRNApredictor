package Performance;
#!/usr/bin/perl
use strict;
use warnings;




####################################################

sub draw_roc_curve {
	
	my $hash_predict = shift or die;
	my $outfile_roc_data = shift or die;
	
	
	my $cutoff_step = 0.001; 
	
	open  OUT,">$outfile_roc_data" or die;
	print OUT "Cutoff\tTOT_P\tTOT_N\tTP\tTN\tFP\tFN\tPr\tAc\tGm\t1-Sp\tSn\tSp\tMcc\n";
	
	my $max_mcc = 0;
	my $threshod = 0.5;
	
	my @arrayX = ();
	my @arrayY = ();
	my $cutoff = 0;
	while ( $cutoff < 1 ) {	
		
		my %hash_count = ();
		&TP_TN_FP_FN($hash_predict,\%hash_count,$cutoff);
		
		my %hash_perf = ();
		&sn_sp_ac_mcc(\%hash_count,\%hash_perf);
	
		my $TP = $hash_count{TP};
		my $TN = $hash_count{TN};
		my $FP = $hash_count{FP};
		my $FN = $hash_count{FN};
		my $TOT_P = $hash_count{total_pos};
		my $TOT_N = $hash_count{total_neg};
		
		my $sn = $hash_perf{SN};
		my $sp = $hash_perf{SP};
		my $ac = $hash_perf{AC};
		my $pr = $hash_perf{PR};
		my $Gm = $hash_perf{GM};
		my $mcc = $hash_perf{MCC};
		
		if($mcc > $max_mcc){
			$max_mcc = $mcc;
			$threshod = $cutoff;
		}
			
		my $t = "";
		if( &is_numeric($sp) ){
			$t = &round(1-$sp,6);
		} else {
			$t = "NaN";
		}
		
		
		push(@arrayX, $t);
		push(@arrayY, $sn);
		
		$cutoff = &round($cutoff,4);
		print OUT "$cutoff\t$TOT_P\t$TOT_N\t$TP\t$TN\t$FP\t$FN\t$pr\t$ac\t$Gm\t$t\t$sn\t$sp\t$mcc\n";
		
		$cutoff += $cutoff_step;

	}
	
	close OUT;
	
	my $auc = &calculate_AUC(\@arrayX,\@arrayY);
	
	return ($threshod,$auc);
}



####################################################

sub TP_TN_FP_FN {
	
	my $hash_predict = shift or die;
	my $hash_count = shift or die;
	my $cutoff = shift;
	
	my $total_pos = 0;	
	my $total_neg = 0; 	

	foreach my $id ( sort keys %$hash_predict ){		
		if( $hash_predict->{$id}{real} == 1 ){ 		
			$total_pos++;		
			if ( $hash_predict->{$id}{proP} >= $cutoff ){ 
				unless ( exists $hash_count->{TP} ){
					$hash_count->{TP} = 1;
				} else {
					$hash_count->{TP}++;
				}		
			} else { 
				unless ( exists $hash_count->{FN} ){
					$hash_count->{FN} = 1;
				} else {
					$hash_count->{FN}++;
				}
			}				
		} elsif ( $hash_predict->{$id}{real} == -1 ){
			$total_neg++;
			if ( $hash_predict->{$id}{proP} >= $cutoff ){ 
				unless ( exists $hash_count->{FP} ){
					$hash_count->{FP} = 1;
				} else {
					$hash_count->{FP}++;
				}	
			} else { 
				unless ( exists $hash_count->{TN} ){
					$hash_count->{TN} = 1;
				} else {
					$hash_count->{TN}++;
				}
			}
		} else {
			print STDERR "Error: label must be 1 or -1 ($hash_predict->{$id}{real})\n";
			die;
		}				
	}
	
	$hash_count->{total_pos} = $total_pos;
	$hash_count->{total_neg} = $total_neg;
	

	unless ( exists $hash_count->{TP} ){
		$hash_count->{TP} = 0;
	}
	unless ( exists $hash_count->{TN} ){
		$hash_count->{TN} = 0;
	}
	unless ( exists $hash_count->{FP} ){
		$hash_count->{FP} = 0;
	}
	unless ( exists $hash_count->{FN} ){
		$hash_count->{FN} = 0;
	}
	
}



####################################################

sub sn_sp_ac_mcc {
	
	my $hash_count = shift or die;
	my $hash_perf  = shift or die;

	my $TP = $hash_count->{TP};
	my $TN = $hash_count->{TN};
	my $FP = $hash_count->{FP};
	my $FN = $hash_count->{FN};
	
	my $sn = "NaN";
	my $sp = "NaN";
	my $ac = "NaN";
	my $pr = "NaN";
	my $Gm = "NaN";
	my $mcc = "NaN";
	
	if ( $TP + $FN > 0 ) {
		$sn = $TP / ( $TP + $FN );
		$sn = &round($sn,6);
	}

	if ( $TN + $FP > 0  ){
		$sp = $TN / ( $TN + $FP );
		$sp = &round($sp,6);
	} 

	if ( $TP + $FP + $TN + $FN > 0 ) {
		$ac = ($TP + $TN) / ( $TP + $FP + $TN + $FN );
		$ac = &round($ac,6);
	}

 	if ( $TP + $FP > 0 ){
 		$pr = $TP / ( $TP + $FP );
 		$pr = &round($pr,6);
 	}
 	
	if( ($TP+$FN) * ($TN+$FP) * ($TP+$FP) * ($TN+$FN) > 0 ){
		$mcc = ( ($TP*$TN) - ($FN*$FP) ) / sqrt( ($TP+$FN) * ($TN+$FP) * ($TP+$FP) * ($TN+$FN) );
		$mcc = &round($mcc, 6);
	}
	
	if ( (&is_numeric($sn)) && (&is_numeric($sp)) && $sn * $sp >= 0 ){
		$Gm = sqrt( $sn * $sp );
		$Gm = &round($Gm,6);
	}
			
	
	$hash_perf->{SN}  = $sn;
	$hash_perf->{SP}  = $sp;
	$hash_perf->{AC}  = $ac;
	$hash_perf->{PR}  = $pr;
	$hash_perf->{GM}  = $Gm;
	$hash_perf->{MCC} = $mcc;
	
	
	my $t = "";
	if( &is_numeric($sp) ){
		$t = &round(1-$sp,6);
	} else {
		$t = "NaN";
	}
	
	
	return "TP:$TP\tTN:$TN\tFP:$FP\tFN:$FN\tPR:$pr\tAC:$ac\tGM:$Gm\t1-SP:$t\tSN:$sn\tSP:$sp\tMCC:$mcc";
	
}


####################################################

sub is_numeric {
	
	use Scalar::Util qw(looks_like_number);

	my $v = shift;

	if( looks_like_number( $v ) ){  	
		return 1;
	} else {
		return 0;
	}
}


####################################################

sub round {
	
	my $usage = "<  numerical value > < the numbers after the decimal point >";
	my $val = shift;
	my $col = shift;

	unless(defined $val){
		die $usage;
	}

	unless(defined $col){
		die $usage;
	}

	unless( &is_numeric($val) ){
		print STDERR "Error:$val not a numeric";
		die;
	}

	unless( &is_numeric($col) ){
		print STDERR "Error:$col not a numeric";
		die;
	}

	my $r = 10 ** $col;
	my $a = ($val > 0) ? 0.5 : -0.5;

	return int($val * $r + $a) / $r;

}

####################################################


sub calculate_AUC {
	
	my $arrayX = shift or die;
	my $arrayY = shift or die;
	
	my $sum_area = 0;
	for (my $i = 0; $i < scalar @{$arrayY} - 1; $i++ ){
		my $l = $arrayX->[$i+1] - $arrayX->[$i];
		my $h = ($arrayY->[$i+1] + $arrayY->[$i]);
		$sum_area += $l * $h / 2;
	}
	
	$sum_area = abs($sum_area);
	$sum_area = &round($sum_area,4);
	
	print "\n\n\n##### AUC: $sum_area #####\n\n\n";
	
	return $sum_area;
}

####################################################
1;