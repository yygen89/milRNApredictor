package RandomForest;
#!/usr/bin/perl
use strict;
use warnings;



####################################################

sub training {

	my $feature_file = shift or die;
	my $model_file = shift or die;
	my $dest_dir = shift or die;
	
	
	unless(-e $feature_file){
		print STDERR "Error:$feature_file not exists\n";
		die;
	}
	
	
	unless(-e $dest_dir){
		mkdir($dest_dir);
	}
	
	
	# ================= Training =================
	
	my $R_bin = "R";
	my $seed = 17;

	my $cmd = "model <- randomForest(type ~ .,  data = training, importance = TRUE, do.trace = 100, proximity = TRUE)";

	my @R_cmds = (
		'library(randomForest)',	
		'training_data <- read.table("'. $feature_file .'", sep = "\t", quote = "\"", header = TRUE)',
		'training <- training_data[1:nrow(training_data),2:ncol(training_data)]',
		'type <- training_data[1:nrow(training_data),ncol(training_data)]',
		'set.seed("'. $seed  .'")',
		$cmd,
		'save(model, file = "'. $model_file  . '")',	
		'rm(model)',
		'q("no")'
	);
	
			
	my $R_script = "$dest_dir/training.R";
	open RSCRIPT, ">$R_script" or die;
	for my $cmd ( @R_cmds ){
	    print RSCRIPT "$cmd\n";
	}
	close RSCRIPT;
	
	
	
	die if system("$R_bin CMD BATCH $R_script");
	

	# ================= Clean =================
	
	unlink(".RData") if -e ".RData";
	unlink("training.Rout") if -e "training.Rout";
	
}
	


####################################################

sub predict {
	
	my $feature_file = shift or die;
	my $model_file = shift or die;
	my $outfile_out = shift or die;
	my $dest_dir = shift or die;	
	

	unless(-e $feature_file){
		print STDERR "Error:$feature_file not exists\n";
		die;
	}
	
	unless(-e $model_file){
		print STDERR "Error:$model_file not exists\n";
		die;
	}	
	
	unless(-e $dest_dir){
		mkdir($dest_dir);
	}
		
	
	# =========== Test ===============
	
	my $R_bin = "R";
	
	my $R_outfile = "$dest_dir/predict.out";
	
	my @R_cmds = (
		'library(randomForest)',
		'load("' . $model_file . '")',
		'test_data <- read.table("'. $feature_file .'", sep = "\t", quote = "\"", header = TRUE)',
		'test <- test_data[1:nrow(test_data),2:ncol(test_data)]',
		'print(model)',
		'pred <- predict(model, test, type = "prob")',
		'write.table(pred, "' . $R_outfile . '", quote = FALSE, sep = "\t")',
		'q("no")'
	);
	
		
	my $R_script = "$dest_dir/predict.R";
	open RSCRIPT, ">$R_script" or die;
	for my $cmd ( @R_cmds ){
	    print RSCRIPT "$cmd\n";
	}
	close RSCRIPT;


	die if system("$R_bin CMD BATCH $R_script ");
	
	unlink(".RData") if -e ".RData";
	unlink("predict.Rout") if -e "predict.Rout";
	
	
	# ============ Parse ===============
	
	# get name 
	my $m = 1;
	my %hash_ID = ();
	my %hash_predict = ();
	open IN,"<$feature_file" or die;
	while(defined(my $line = <IN> )){
		chomp($line);
		next if $line =~/^\s*$/;
		next if $line =~/^Name/i;
		my @tmp = split /\t/,$line;
		my $name = shift( @tmp );
		$hash_ID{ $m } = $name;
		$m++;
	}
	close IN;
	
	
	my %hash_row = ();
	my @classNames = &toMatrix($R_outfile,\%hash_row,'name'=>"Name");
	
		
	my $title = "Name\tpredictType";
	for my $class ( @classNames ){
		$title .= "\t$class";
	}
	
	open  RFOUT,">$outfile_out" or die;
	print RFOUT "$title\n";

	my $positive_score_threshold = 0.5;
	for my $row ( sort {$a <=> $b} keys %hash_row ){
		
		my $name = $hash_ID{$row};
		my $str = "$name";
		
		# ============== predicted type ================
			
		my $max_class;
		my $max_prob = 0;
		
		for my $class ( @classNames ){
			unless(exists $hash_row{$row}{$class}){
				print STDERR "Error:($row,$class)\n";
				die;
			}
			my $prob = $hash_row{$row}{$class}; 
			if ( $prob > $max_prob ){
				$max_class = $class;
				$max_prob = $prob;
			}
		}
		
		if ( $max_class eq "positive" ){
			if ( $max_prob >= $positive_score_threshold ){
				$max_class = 1;
			}else{
				$max_class = -1;
			}
		} elsif ( $max_class eq "negative" ){
			$max_class = -1;
		} else {
			print STDERR "Error:$max_class must be positive or negative\n";
			die;
		}
		$str .= "\t$max_class";
			
		# =============== output ===========
		
		for my $class ( @classNames ){
			my $prob = $hash_row{$row}{$class}; 
			$str .= "\t$prob";
		}
		
		print RFOUT "$str\n";
			
		unless(exists $hash_predict{$name}{proP}){
			unless(exists $hash_row{$row}{positive}){
				print STDERR "Error:($row,positive) not defined\n";
				die;
			}
			$hash_predict{$name}{proP} = $hash_row{$row}{positive};
		} else {
			print STDERR "Error:$name repeat\n";
			die;
		}
	}
	
	close RFOUT;

}



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
	
	open IN,"<$infile" or die;
	
	my $flag = 1;
	my @columns = ();
	my %hash_col_tmp = ();
	
	while( defined( my $line = <IN> )){
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
	close IN;
	
	return @columns;
	
}


####################################################
1;