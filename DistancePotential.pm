package DistancePotential;
#!/usr/bin/perl
use strict;
use warnings;

#################
use WordSet;
#################

my $dr = 2;
my $nbin = 20;

my @array_kmers = ();
my @array_kmer_nmismatch = ("2#0","3#0","4#0");

my %hash_kmer_pair = ();
my %hash_kmer_pair_length = ();

#############################################################

sub extract_feature {
	
	my $training_positive_fasta = shift or die;
	my $training_negative_fasta = shift or die;
	my $test_fasta = shift or die;
	
	my $hash_feature_positive = shift or die;
	my $hash_feature_negative = shift or die;
	my $hash_feature_test = shift or die;

	my $save_potential_file = shift;
	
	# ========== Construct the set of kmers =============
	
	unless( @array_kmers ){
		print "\n======  Construct the set of kmers ======\n\n";
		@array_kmers = WordSet::get_mismatch_profile(\@array_kmer_nmismatch);
	}
	

	# ========== Calculate statistical potentials =============
	
	print "\n========  Calculate statistical potentials =======\n\n";
	my %hash_potential = ();
	&construct_potential($training_positive_fasta,$training_negative_fasta,\%hash_potential);
	
	
	# ================ Save potential =================
	
	if(defined $save_potential_file){
		&save_potential(\%hash_potential,$save_potential_file);
	}
	
	# ================ Extract features ================
	
	print "\n======  Extract statistical potential features ======\n\n";
	&get_feature($training_positive_fasta,\%hash_potential,$hash_feature_positive);
	&get_feature($training_negative_fasta,\%hash_potential,$hash_feature_negative);
	&get_feature($test_fasta,\%hash_potential,$hash_feature_test);

}


sub get_feature {
	
	my $fasta_file = shift or die;
	my $hash_potential = shift or die;
	my $hash_feature = shift or die;
	
	
	#assign initial values for global variables
	unless( @array_kmers ){
		@array_kmers = WordSet::get_mismatch_profile(\@array_kmer_nmismatch);
	}
	
	
	unless( scalar keys %hash_kmer_pair ){
		for my $d ( sort {$a<=>$b} keys %{ $hash_potential } ){
			for my $kp ( sort keys %{ $hash_potential->{$d} } ){
				$hash_kmer_pair{$kp} = 1;	
				my $lng = length($kp) - 1;
				$hash_kmer_pair_length{$lng} = 1;
			}
		}
	}
	
	my $dif_lng_num = 0;
	for my $lng ( keys %hash_kmer_pair_length ){
		$dif_lng_num++;
	}

	my %hash_seq = ();
	&readfasta($fasta_file,\%hash_seq);
	
	for my $id ( keys %hash_seq ){
	
		#obtain distance-dependent kmer pairs
		my %hash_dist = ();
		&get_kmer_pair_frequency($hash_seq{$id},\%hash_dist);	
		
		#extract features
		my $total_score = 0;
		my %hash_score_lng = ();
		my %hash_score_dist = ();
		for(my $d = 0; $d < $nbin; $d++) {
			my $score_d = 0;
			for my $kp ( keys %hash_kmer_pair ){			
				if( (exists $hash_potential->{$d}{$kp}) && (exists $hash_dist{$d}{$kp}) ){
					my $v = $hash_potential->{$d}{$kp};	
					$total_score += $v;
					$score_d += $v;
					
					my $lng = length($kp) - 1;
					unless(exists $hash_score_lng{$lng}){
						$hash_score_lng{$lng} = $v;
					}else{
						$hash_score_lng{$lng} += $v;
					}
				}
			}
			$hash_score_dist{$d} = $score_d;
		}
			
		#record features	
		$hash_feature->{$id}{"dp"} = $total_score;
		
		for(my $d = 0; $d < $nbin; $d++) {
			$hash_feature->{$id}{"dpD-$d"} = $hash_score_dist{$d};
		}
	
	
		if( $dif_lng_num > 1 ){
			for my $lng ( keys %hash_kmer_pair_length ){
				my $v = 0;
				if(exists $hash_score_lng{$lng}){
					$v = $hash_score_lng{$lng};
				}
				$hash_feature->{$id}{"dpL-$lng"} = $v;
			}
		}
		
	}
}


sub construct_potential {
	
	my $training_positive_fasta = shift or die;
	my $training_negative_fasta = shift or die;
	my $hash_potential = shift or die;

	
	# ================= Positive ==================

	my %hash_seq_positive = ();
	&readfasta($training_positive_fasta,\%hash_seq_positive);

	my %hash_dist_positive = ();
	for my $id ( keys %hash_seq_positive ){
		&get_kmer_pair_frequency($hash_seq_positive{$id},\%hash_dist_positive);
	}
	
	# ================= Negative ==================
	
	my %hash_seq_negative = ();
	&readfasta($training_negative_fasta,\%hash_seq_negative);
	
	my %hash_dist_negative = ();
	for my $id ( keys %hash_seq_negative ){
		&get_kmer_pair_frequency($hash_seq_negative{$id},\%hash_dist_negative);
	}
	
	#all possible kmer pairs
	my %hash_all_kmer_pairs = ();
	for my $kmer_i ( @array_kmers ){
		for my $kmer_j ( @array_kmers ){
			my ($k1,$k2) = sort ($kmer_i,$kmer_j);
			my $kp = "$k1-$k2";
			$hash_all_kmer_pairs{$kp} = 1;
		}
	}
	
	# ====== distance-dependent potential  ===========

	for(my $d = 0; $d < $nbin; $d++) {
		for my $kp ( keys %hash_all_kmer_pairs ){	
			if( (exists $hash_dist_positive{$d}{$kp}) && (exists $hash_dist_negative{$d}{$kp}) ){
				my $r = $hash_dist_positive{$d}{$kp}/$hash_dist_negative{$d}{$kp};
				my $v = -log($r);
				if( abs($v) > 0 ){
					$hash_potential->{$d}{$kp} = $v;
				}
			}
		}
	}
}


sub get_kmer_pair_frequency {
	
	my $seq = shift or die;
	my $hash_dist = shift or die;
	
	
	my %hash_freq = ();
	for my $kmer ( @array_kmers ){
		my $nmismatch = $kmer =~ tr/*/*/;  
		my %hash_coordinate = ();
		my $n_hit = &matchPattern($nmismatch,$kmer,$seq,\%hash_coordinate);
		for my $b ( keys %hash_coordinate ){
			$hash_freq{$b}{$kmer} = 1;
		}
	}
	
	#calculate distance between two kmers 
	my %hash_distribution = ();
	for(my $i = 0; $i < length($seq)-1; $i++) {
		for(my $j = $i + 1; $j < length($seq); $j++) {
			
			#calculate distance distribution
			my $d = $j - $i; 
			my $r = 0;
			while( 1 ){
				last if $r <= $d && $d < $r + $dr;
				$r += $dr;
			}
			
			for my $kmer_i ( keys %{$hash_freq{$i}} ){
				for my $kmer_j ( keys %{$hash_freq{$j}} ){
					my ($k1,$k2) = sort ($kmer_i,$kmer_j);
					my $kp = "$k1-$k2";
					unless(exists $hash_dist->{$r}{$kp}){
						$hash_dist->{$r}{$kp} = 1;
					}else{
						$hash_dist->{$r}{$kp} += 1;
					}
				}
			}
		}
	}
}


sub save_potential{
	
	my $hash_potential = shift or die;
	my $outfile = shift or die;
	
	open POT,">$outfile" or die;
	for my $d ( sort {$a<=>$b} keys %{ $hash_potential } ){
		for my $kp ( sort keys %{ $hash_potential->{$d} } ){
			my $v = $hash_potential->{$d}{$kp};
			print POT "$d\t$kp\t$v\n";
		}
	}
	close POT;
	
}


sub readin_potential{
	
	my $infile = shift or die;
	my $hash_potential = shift or die;
	
	open IN,"<$infile" or die;
	while(my $line = <IN>){
		chomp($line);
		next if $line =~ /^\s*$/;
		
		my @tmp = split/\t/,$line;
		
		my $d = $tmp[0];
		my $kp = $tmp[1];
		my $v = $tmp[2];
		
		unless(exists $hash_potential->{$d}{$kp}){
			$hash_potential->{$d}{$kp} = $v;
		}else{
			print "Error: ($d,$kp) repeat\n";
			die;
		}
	}
	close IN;
	
}




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
			my $matched_pattern = substr( $subjectSeq, $_, $lng );
			if( defined $hash_pos ){
				$hash_pos->{ $_ } = $matched_pattern;
			}	
			$n++;
		}
	}
	
	return $n;
}



sub logN { 
	
	my $v = shift;
	my $b = shift;

	if ( $b <= 0 ){
		print STDERR "Error:base number ($b) <= 0\n";
		die;
	}
	
	if ( $v <= 0 ){
		print STDERR "Error:input value ($v) <= 0\n";
		die;
	} else {
		return log($v)/log($b);    
	}
}

####################################################
1;