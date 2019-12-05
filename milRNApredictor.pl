#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use Getopt::Long  qw(:config bundling);
use File::Basename;
use lib dirname(__FILE__);
##################################################
use FeatureSelect;
use RandomForest;
use Performance;
use DistancePotential;
use KmerS;
##################################################
my $is_KmerS = 1;
my $is_distance_potential = 1;
##################################################
my $dir;
my $help;
my $positive_fasta;
my $negative_fasta;
my $test_fasta;
##################################################
GetOptions( 	
	"help|h!" => \$help,
	"posFasta|P=s" => \$positive_fasta,
	"negFasta|N=s" => \$negative_fasta,
	"testFasta|T=s" => \$test_fasta,
	"dir|D=s" =>\$dir,
);
##################################################

if ( $help ) { &usage(); exit; }
unless( defined $positive_fasta ){&usage(); exit;}
unless( defined $negative_fasta ){&usage(); exit;}
unless( defined $test_fasta ){&usage(); exit;}

my $current_dir = getcwd;
unless(defined $dir){$dir = "$current_dir/result_dir";}
mkdir($dir) unless -e $dir;

my $model_dir = "$current_dir/model";
unless(-e $model_dir ){
	mkdir( $model_dir );
}

##################################################

my $start_time = &start();

&train_predict($positive_fasta,$negative_fasta,$test_fasta,$model_dir,$dir);

&end($start_time);

##################################################

sub train_predict {
	
	my $training_positive_fasta = shift or die;
	my $training_negative_fasta = shift or die;
	my $test_fasta = shift or die;
	my $model_dir = shift or die;
	my $dir = shift or die;

	
	# ==========  Feature extraction ==========

	my %hash_feature_positive = ();
	my %hash_feature_negative = ();
	my %hash_feature_test = ();
	
	if( $is_distance_potential ){
		my $potential_model_file = "$model_dir/potential.model";
		DistancePotential::extract_feature($training_positive_fasta,$training_negative_fasta,$test_fasta,
		\%hash_feature_positive,\%hash_feature_negative,\%hash_feature_test,$potential_model_file);
	}
	
	if( $is_KmerS ){
		KmerS::extract_feature($training_positive_fasta,$training_negative_fasta,$test_fasta,
		\%hash_feature_positive,\%hash_feature_negative,\%hash_feature_test);
	}
	
	# ==========  Feature selection ==========

	my $train_feature = "$dir/train_feature.txt";
	my $test_feature = "$dir/test_feature.txt";
	FeatureSelect::feature_selection(\%hash_feature_positive,\%hash_feature_negative,
	\%hash_feature_test,$train_feature,$test_feature,$dir);
	

	# ============= Traning ==============
	
	my $training_model = "$model_dir/RF_model.RData";
	RandomForest::training($train_feature,$training_model,$dir);

	# ============= Testing ==============
	
	my $outfile_out = "$dir/predict_results.txt";
	RandomForest::predict($test_feature,$training_model,$outfile_out,$dir);

}


####################################################

sub start {
	
	my $out = shift; 

	my($second, $minute, $hour, $dayOfMonth, $month, 
	$yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	$second = "0$second" if($second =~ /^\d$/);
	my $sTime = "$hour:$minute:$second";
	my $stime = time;
	my $start_time = localtime;

	print "\n\nstarted: $start_time\n\n";
	if( defined $out ){
		print $out "started: $start_time\n";
	}
    
	return $stime;
}


####################################################

sub end {
	
	my $stime = shift or die;
	my $out = shift;

	my $etime = time - $stime;
	my ($second, $minute, $hour, $dayOfMonth, $month, 
	$yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	$second = "0$second" if($second =~ /^\d$/);
	my $eTime = "$hour:$minute:$second";

	my $end_time = localtime;
    
	print "\n\nended: $end_time\n";
	print "total:", int($etime / 3600),"h:",
	int(($etime % 3600) / 60),"m:",int($etime % 60),"s\n\n";

	if( defined $out ){
		print $out "\n\nended: $end_time\n";
		print $out "total:", int($etime / 3600),"h:",
		int(($etime % 3600) / 60),"m:",int($etime % 60),"s\n\n";
	}

	print "\n\n============== End ============\n\n";
}


####################################################

sub usage {
	
my $usage = << "USAGE";

Program: $0
Contact: Yao Yuangen <yygen89\@163.com>

Usage:

	+++++++++++++++++++++++++++++++++++++++++
		
	Options:
	--posFasta	| -P [file] : fasta file containing positive sequences
	--negFasta	| -N [file] : fasta file containing negative sequences
	--testFasta	| -T [file] : fasta file containing test sequences
	--dir		| -D [file] : file containing calculated results (default: results)
	--help		| -h : help
	
	+++++++++++++++++++++++++++++++++++++++++
	
	$0 [options] -P <postive fasta> -N <negative fasta> -T <test fasta> -D [directory] 

	+++++++++++++++++++++++++++++++++++++++++ 

USAGE
print $usage;

}
