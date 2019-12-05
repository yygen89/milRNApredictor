package FeatureSelect;
#!/usr/bin/perl
use strict;
use warnings;


my $is_feat_selct = 0;

#############################################################


sub feature_selection {
	
	my $hash_feature_positive = shift or die;
	my $hash_feature_negative = shift or die;
	my $hash_feature_test = shift or die;
	my $selected_train_feature = shift or die;
	my $selected_test_feature = shift or die;
	my $dir = shift or die;


	my $flag_no_selected = 1;
	if( $is_feat_selct ){
		
		print "\n****** Feaction slection ******\n\n";
		
		my $train_feature = "$dir/train_all_feature.txt";
		open FEATUREOUTTRAIN,">$train_feature" or die;
		&print_feature($hash_feature_positive,\*FEATUREOUTTRAIN,'label'=>"positive",'title'=>1);
		&print_feature($hash_feature_negative,\*FEATUREOUTTRAIN,'label'=>"negative");
		close FEATUREOUTTRAIN;
		
		
		my $feature_model = "$dir/features.model";	
		my $featset_num = &rrf($train_feature,$selected_train_feature,$feature_model,$dir);
		
		
		if( $featset_num > 1 ){
			
			my $test_feature = "$dir/test_all_feature.txt";
			open FEATUREOUTTEST,">$test_feature" or die;
			&print_feature($hash_feature_test,\*FEATUREOUTTEST,'title'=>1);
			close FEATUREOUTTEST;
			
			&pick_feature_by_optimized_results($test_feature,$selected_test_feature,$feature_model);
			
			print "\nThe number of features: $featset_num\n";
			
			$flag_no_selected = 0;
		}
	}
	
	
	if( $flag_no_selected ){
		
		open FEATUREOUTTRAI,">$selected_train_feature" or die;
		my $featset_num = &print_feature($hash_feature_positive,\*FEATUREOUTTRAI,
		'label'=>"positive",'title'=>1);
		&print_feature($hash_feature_negative,\*FEATUREOUTTRAI,'label'=>"negative");
		close FEATUREOUTTRAI;
		
		open FEATUREOUTTEST,">$selected_test_feature" or die;
		&print_feature($hash_feature_test,\*FEATUREOUTTEST,'title'=>1);
		close FEATUREOUTTEST;
		
		print "\nThe number of features: $featset_num\n";
	}
}


#############################################################

sub print_feature {
	
	my $hash_feature = shift or die;
	local *FEATUREOUT = shift or die;
	
	
	my $label;
	my $title_flag = 0;
	my $fscore_file;
	my $fscore_threshold;
	while (@_) {
		my $argument = shift @_;	
		if ($argument=~/label/i) {$label=shift @_}
		if ($argument=~/title/i) {$title_flag=shift @_}
		if ($argument=~/fscoreFile/i) {$fscore_file=shift @_}
		if ($argument=~/threshold/i) {$fscore_threshold=shift @_}
	}
	
	
	my %hash_feat_name = ();
	for my $id ( keys %{$hash_feature} ){
		for my $feature ( keys %{$hash_feature->{$id}} ){
			$hash_feat_name{$feature} = 1;
		}
	}
	my @array_feats = sort keys %hash_feat_name;
	
		
	
	if( $title_flag ){	
		my $title = "Name";
		for my $feature ( @array_feats ){
			$title .= "\t$feature";
		}
		
		if(defined $label){
			print FEATUREOUT "$title\ttype\n";
		}else{
			print FEATUREOUT "$title\n";
		}
	}
	
	
	for my $id ( sort {$a cmp $b} keys %{$hash_feature} ){
		my $str = $id;
		for my $feature ( @array_feats ){
			my $value = $hash_feature->{$id}{$feature};
			$value = &round($value,4);
			$str .= "\t$value";
		}
		
		if(defined $label){
			print FEATUREOUT "$str\t$label\n";
		}else{
			print FEATUREOUT "$str\n";
		}	
	}
	
	my $featset_num = scalar @array_feats;
	return $featset_num;
}

#############################################################

sub rrf {

	my $feature_file= shift or die;
	my $selected_feature_file = shift or die;
	my $feature_model  = shift or die;
	my $dest_dir = shift or die;
	
	mkdir($dest_dir) unless -e $dest_dir;
	
	# ================= Training =================
	
	
	my @R_cmds = (
		'library(RRF)',	
		'set.seed(1)',
		
		'training_data <- read.table("'. $feature_file .'", sep = "\t", quote = "\"", header = TRUE)',
		'X<- training_data[1:nrow(training_data),3:ncol(training_data)-1]',
		'class <- training_data[1:nrow(training_data),ncol(training_data)]',
		
		'rf  <-  RRF(X,as.factor(class),  flagReg  =  0)',
		'impRF  <-  rf$importance',
		'impRF  <-  impRF[,"MeanDecreaseGini"]',
		'rf$feaSet',
		
		'imp=impRF/(max(impRF))',
		'coefReg=0.9*0.8+0.1*imp',
		'rrf <- RRF(X,as.factor(class),mtry=ncol(X),coefReg=coefReg)',
		'imp=rrf$importance',
		'imp=imp[,"MeanDecreaseGini"]',
		'chosed_feaSet <- rrf$feaSet',
		'chosed_feaSet',
		'write.table(chosed_feaSet, "' . $feature_model . '", quote = FALSE, sep = "\t")',
		'q("no")'
	);
	
			
	my $R_script = "$dest_dir/RRF.R";
	open RSCRIPT, ">$R_script" or die;
	for my $cmd ( @R_cmds ){
	    print RSCRIPT "$cmd\n";
	}
	close RSCRIPT;
	
	die if system("R CMD BATCH $R_script");
	
	
	# ================= Clean =================
	
	unlink(".RData") if -e ".RData";
	unlink("RRF.Rout") if -e "RRF.Rout";
	
	# ==================================
	
	my $featset_num = &pick_feature_by_optimized_results($feature_file,$selected_feature_file,$feature_model);
	
	# ==================================
	
	return $featset_num;
}



#############################################################

sub pick_feature_by_optimized_results {

	my $feature_file = shift or die;
	my $selected_feature_file = shift or die;
	my $feature_model = shift or die;

	
	# ============ Parse ===============
	
	# get feature ID
	my %hash_feat_ID = ();
	open RRFIN,"<$feature_model" or die;
	while(defined(my $line = <RRFIN> )){
		chomp($line);
		next if $line =~/^\s*$/;
		next if $line =~/^x$/;
		my @tmp = split/\t/,$line;
		my $feat_ID = pop(@tmp);
		$hash_feat_ID{$feat_ID} = 1;
	}
	close RRFIN;
	
	my $featset_num = scalar keys %hash_feat_ID;
	
	# =================================
	
	
	my %hash_row = ();
	my @classNames = &toMatrix($feature_file,\%hash_row);
	
	
	my @featsets = ();
	for my $feat_ID ( sort {$a<=>$b} keys %hash_feat_ID ){
		my $feat_name = $classNames[$feat_ID];
		push(@featsets,$feat_name);
	}
	
	
	
	my $title = "Name";
	for my $feat_name ( @featsets ){
		$title .= "\t$feat_name";
	}
	
	open RRFOUT,">$selected_feature_file" or die;
	
	if($classNames[-1] eq "type"){
		print RRFOUT "$title\ttype\n";
	}else{
		print RRFOUT "$title\n";
	}
	
	for my $name ( sort keys %hash_row ){
		my $str = $name;
		for my $feat_name ( @featsets ){
			my $value = $hash_row{$name}{$feat_name};
			$str .= "\t$value";
		}
		
		if(exists $hash_row{$name}{"type"}){
			my $type = $hash_row{$name}{"type"};
			print RRFOUT "$str\t$type\n";
		}else{
			print RRFOUT "$str\n";
		}
	}
	close RRFOUT;
	
	return $featset_num;
}

	
#############################################################

sub toMatrix {
	
	my $infile = shift or die;
	my $hash_row = shift or die;
	
	my $name;
	while (@_) {
		my $argument = shift @_;
		if ($argument=~/name/i) {$name=shift @_}
	}
	

	unless(-e $infile){
		print STDERR "Error:$infile not exists\n";
		die;
	}
	
	open INFEAT,"<$infile" or die;
	
	my $flag = 1;
	my @columns = ();
	my %hash_col_tmp = ();
	
	while( defined( my $line = <INFEAT> )){
		chomp($line);
		next if $line =~/^\s*$/;
		
		my @tmp = split/\t/,$line;
		my $row = shift @tmp;
		
		if ( $flag ) {
			if ( defined $name ){
				unshift(@tmp,$row);
			} else {	
				push(@columns,$row);
			}
			for my $i ( 0 .. $#tmp ){		
				unless( exists $hash_col_tmp{ $i } ){
					$hash_col_tmp{ $i } = $tmp[$i];
				}	
				push(@columns,$tmp[$i]);		
			}
			$flag = 0;
			next;
		} 
		
		for my $i ( 0 .. $#tmp ){
			my $col = $hash_col_tmp{ $i };
			unless( exists $hash_row->{ $row }{ $col } ){	
				$hash_row->{ $row }{ $col } = $tmp[$i];
			} else {
				print STDERR "Error:($row,$col) repeat\n";
				die;
			}
		}
	}
	close INFEAT;
	
	return @columns;
}


sub round {
	
	my $usage = "<  numerical value > < the numbers after the decimal point >";
	my $val = shift;
	my $col = shift;

	my $r = 10 ** $col;
	my $a = ($val > 0) ? 0.5 : -0.5;

	return int($val * $r + $a) / $r;
}


#############################################################
1;