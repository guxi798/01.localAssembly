my $searchio = new Bio::SearchIO( -format => 'blast', -file => $blast_output);

# Go through each report (= each result from one query)
while(my $result = $searchio->next_result){
    # process the Bio::Search::Result::ResultI object
    my $queryName = $result->query_name;
    my $queryAcc = $result->query_accession;
    my $queryLength = $result->query_length;
    my $queryDesc = $result->query_description;
    my $dbName = $result->database_name;
    my $numSeqsInDb = $result->database_entries;
    my $numHits = $result->num_hits;


    # Go through each each matching sequence for this query
    while(my $hit = $result->next_hit){
	# process the Bio::Search::Hit::HitI object
        my $hitName = $hit->name;
	my $hitAcc = $hit->accession;
	my $hitDesc = $hit->description;
	my $hitEvalue = $hit->significance;
	my $hitBits = $hit->bits;
	my $numHsps = $hit->num_hsps;


	# Go through each each HSP for this sequence
	while (my $hsp = $hit->next_hsp){
	    # process the Bio::Search::HSP::HSPI object
	    my  $hspEvalue = $hsp->evalue, "\n";
	    my  $fracIdentical = $hsp->frac_identical;
	    my  $fracConserved = $hsp->frac_conserved;
	    my  $numGaps = $hsp->gaps;
	    my  $hspLength = $hsp->hsp_length;
	    my  $hspQuerySeq = $hsp->query_string;
	    my  $hspHitSeq = $hsp->hit_string;
	    my  $hspConsensus = $hsp->homology_string;
	    my  $hspLength = $hsp->hsp_length;
	    my  $hspRank = $hsp->rank;

	    my  $queryStrand = $hsp->strand('query');
	    my  $hitStrand = $hsp->strand('hit');
	    my  $queryStart = $hsp->start('query');
	    my  $queryEnd = $hsp->end('query');
	    my  $hspStart = $hsp->start('hit');
	    my  $hspEnd = $hsp->end('hit');
	    
	    my  $hspScore = $hsp->score;
	    my  $hspBits = $hsp->bits;
	}
    }
}
