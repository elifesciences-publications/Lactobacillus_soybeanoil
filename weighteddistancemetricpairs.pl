#pull out of the distance matrix the distance between PMA and nonPMA pairs for the same sample 
#create a table with the relevant information from the mapping file
open (weighted, "</output/open_ref_otus/wf_bdiv_40000_AP4only_HF/weighted_unifrac_dm.txt") or die "Can't open weighted data file\n";
open(weightedout, ">/output/open_ref_otus/wf_bdiv_40000_AP4only_HF/weighted_unifrac_dm_PMApairs.txt") or die "Can't open weighted data file\n";
#open (weighted, "</output/open_ref_otus/wf_bdiv_40000_AP4only_LF/weighted_unifrac_dm.txt") or die "Can't open weighted data file\n";
#open(weightedout, ">/output/open_ref_otus/wf_bdiv_40000_AP4only_LF/weighted_unifrac_dm_PMApairs.txt") or die "Can't open weighted data file\n";
open (mappingfile, "</mappingfile.txt") or die "Can't open mapping file\n";

while (my $line = <mappingfile>) {
	chomp $line;
	#print $line, "\n";	
	if ($line =~ "#SampleID") {
		$headerline = $line;
	}
	elsif ($line !~ "#SampleID") {
		my @line = split("\t", $line);
		my $SampleID = @line[0];
		#print $SampleID, "\n";
		my @SampleID = split('\.', $sample);
		my $samplename= @sample[0], "\n";
		my $FAMI = @line[5];
		my $Cage = @line[8];
		my $Diet = @line[11];
		my $Gavage = @line[25];
		my $PMA = @line[27];
		my $Description = @line[38];
		my $FAMIDiet = @line[36];
		my $SequenceCounts = @line[35];
		if ($FAMI eq "FAMIAP4" and $SequenceCounts >40000 ) {  #so only look at the FAMI samples
			$mappingfile{$SampleID} = [$FAMI,$Cage,$Diet,$Gavage,$PMA, $Description, $FAMIDiet];
			#print $SampleID;	
		}
	}
}



close(mappingfile);
		
#read in file
while (my $line = <weighted>) {
	#print $line;
	chomp $line;
	my @line = split("\t",$line);
	if ($line =~ "#") {
		#print $line;
		
		#print @line, "\n";
		#need loop here
		for ($index =1; $index<scalar(@line);$index++) {
			#print $index, "\n";
			my $sample = @line[$index];
			$samples{$sample} = [$index, "NA", @{$mappingfile{$sample}}];
			#print $samples{$sample}, "\n";
		}
	}
	else {
		my $sample2 = @line[0];
		my @sample2 = split('\.', $sample2);
		$sample2name= @sample2[0], "\n";
		$sample2PMA = @sample2[1];
		foreach $sample (keys %samples) {
			my @sample = split('\.', $sample);
			$samplename= @sample[0], "\n";
			$samplePMA = @sample[1];
			#print $samplename, "\t", $samplePMA, "\t", "\n";
			if ($samplename eq $sample2name and $samplePMA ne $sample2PMA) {
				$index = @{$samples{$sample}}[0];
				@{$samples{$sample}}[1]= @line[$index];
				#print $samplename, "\t", @line[$index], "\n";
				
				
			
			} 
		}
	}
}
	print weightedout "sample\tindex\tWeightedDistance\tFAMI\tcage\tDiet\tGavage\tignore\tFAMI_Diet_Gavage\tFAMI_Diet\n";

foreach $sample (keys %samples) {
	my @sample = split('\.', $sample);
	$samplename= @sample[0], "\n";
	$,="\t";
	print weightedout $samplename, "\t", @{$samples{$sample}}, "\n";
}

close(weighted);
close(weightedout);
