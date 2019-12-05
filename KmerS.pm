package KmerS;
#!/usr/bin/perl
use strict;
use warnings;

##############
use WordSet;
##############


my @array_kmers = ();
my @array_kmer_nmismatch = ("2#0","3#0");

###############################################

sub extract_feature {
	
	my $training_positive_fasta = shift or die;
	my $training_negative_fasta = shift or die;
	my $test_fasta = shift or die;
	
	my $hash_feature_positive = shift or die;
	my $hash_feature_negative = shift or die;
	my $hash_feature_test = shift or die;
	
	
	# ========== extract features =============
	
	print "====== extract k-mer features ========\n";
	&get_feature($training_positive_fasta,$hash_feature_positive);
	&get_feature($training_negative_fasta,$hash_feature_negative);
	&get_feature($test_fasta,$hash_feature_test);
	
}


###############################################

sub get_feature{
	
	my $fasta = shift or die;
	my $hash_feature = shift or die;
	
	
	unless( @array_kmers ){
		@array_kmers = WordSet::get_mismatch_profile(\@array_kmer_nmismatch);
	}
	
	my %hash_seq = ();
	&readfasta($fasta,\%hash_seq);
	
	for my $id ( sort keys %hash_seq ){
		for my $kmer ( @array_kmers ){
			my $nmismatch = $kmer =~ tr/*/*/;  
			my $hit_count = &matchPattern($nmismatch,$kmer,$hash_seq{$id});
			$hash_feature->{$id}{$kmer} = $hit_count;
		}
	}
	
}
			

####################################################

sub readfasta {
	
	my $usage = "< fasta file > < hash reference >";
	my $infile = shift or die $usage;
	my $hash_seq = shift or die $usage;
	
	unless(-e $infile){
		print STDERR "Error:$infile not exists";
		die;
	}
	
	open IN, $infile || die;
	
	my $c=0;
	my $seqId;
	while (defined (my $line=<IN>)) {
		chomp($line);
		next if $line=~/^\s*$/;
		if ($line =~/^>/) {
			$line =~s/\s*$//g;	
			$line =~s/^\s*//g;	
			$seqId = substr($line,1);
			$c++;
		} else {
			$line =~s/\s*//g;
			$hash_seq->{$seqId}.=$line;
		}
	}
	close IN;
	
	for my $id ( keys %{$hash_seq} ){
		my $seq = $hash_seq->{$id};
		$seq = uc($seq);
		$seq =~ tr/T/U/;
		$hash_seq->{$id} = $seq;
	}
	
	return $c;
	
}


####################################################

sub matchPattern {
	
	my $nmismatch = shift;
	my $pattern = shift or die;
	my $subjectSeq = shift or die;
	my $hash_pos = shift;
	
	my @pattern = split '', $pattern;
	my @subjectSeq = split '', $subjectSeq;

	my $lng = length( $pattern );
	my %err_count;
	for my $i ( 0 .. @subjectSeq - $lng ) {
		my $n_err = 0 ;
		for my $j ( 0 .. @pattern - 1 ) {
			$n_err++ if( $pattern[$j] ne $subjectSeq[$i+$j] );
		}
		next if ( $n_err > $nmismatch );
		$err_count{$i} = $n_err;
	}

	my $n = 0;
	foreach( keys %err_count ) {
		if( $err_count{$_} == $nmismatch ){
			my $value = substr( $subjectSeq, $_, $lng );
			if( defined $hash_pos ){
				$hash_pos->{ $_ } = $value;
			}	
			$n++;
		}
	}
	
	return $n;
}


# =============================
# %(G+C) = (|C|+|G|)/L * 100
# =============================

sub GC_content {	
	
	my $seq = shift or die;
	
	my @tmp = split ( '', $seq );
	
	my $n = 0;
	my $m = 0;	
	for my $base ( @tmp ) {
		if ( $base eq "G" or $base eq "g" ) {
			$m++;
		} elsif ( $base eq "C" or $base eq "c" ) {
			$m++;
		}
		$n++;
	}
	
	return ($m/$n) * 100;
}

####################################################
1;


