package WordSet;
#!/usr/bin/perl
use strict;
use warnings;


####################################################

sub get_mismatch_profile {
	
	my $array_kmers_nmismatchs = shift or die;
	my $array_character = shift;


	unless(defined $array_character){
		@{$array_character} = qw(A C G U);
	}


	my %hash_wordset = ();
	for my $t ( @{$array_kmers_nmismatchs} ){	
		my ($word_length,$nmismatch) = split/#/,$t;		
		if( $word_length > $nmismatch ){				
			my @words = ();
			&permutation($array_character,\@words,$word_length-1);
			for my $word ( sort @words ){
				&get_wordset($nmismatch,$word,\%hash_wordset);
			}
		}
	}	
	
	my $c = 0;
	my @word_set = ();
	for my $word ( sort keys %hash_wordset ){
		push(@word_set,$word);
		$c++;
	}
	
	print "\nNumber of words: $c\n";
	
	return @word_set;
}


sub get_wordset {
	
	my $nmismatch = shift;
	my $substring = shift or die;
	my $hash_wordset = shift or die;
	
	#get character set for input string.
	my %hash_character = ();
	$hash_character{"*"} = 1; # arbitrary Character
	my @arry_character = split '', $substring;
	for my $i ( 0 .. @arry_character - 1 ) {
		$hash_character{$arry_character[$i]} = 1;
	}
	
	my @wordsets = ();
	my @characters = sort keys %hash_character;
	&permutation(\@characters,\@wordsets,length($substring) - 1);

	for my $word ( sort @wordsets ){		
		if( &matchPattern($nmismatch,$substring,$word)){
			my $count = 0;
			while( $word =~ /\*/g ){
		 		$count ++;
			}	
			if( $count == $nmismatch ){
				$hash_wordset->{$word} = $count;
			}
		}
	}
}


sub permutation {
	
	my ($in_ref_array,$out_ref_array,$n,$prefix) = @_;
	
	unless( defined $prefix ){
		$prefix = "";
	}
	
	if( $n == 0 ){
		for my $r ( @$in_ref_array ){
			push (@$out_ref_array,$prefix.$r);
		}
		return;
   	}
   
   	for my $r ( @$in_ref_array ){
		&permutation($in_ref_array,$out_ref_array,$n-1,$prefix.$r);
	}
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
			my $value = substr( $subjectSeq, $_, $lng );
			if( defined $hash_pos ){
				$hash_pos->{ $_ } = $value;
			}	
			$n++;
		}
	}
	
	return $n;
}

####################################################
1;


