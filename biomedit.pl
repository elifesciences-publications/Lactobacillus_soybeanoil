#only use samples present with at least 100 counts total and in at least 25% of the samples and at least

#run either the LF samples or the HF samples

#open (biomfile, "</output/open_ref_otus/wf_bdiv_40000/otu_table_mc2_w_tax_even40000_AP4only_LF_PMA.biom.txt") or die "Can't open biom file\n";
open (biomfile, "</output/open_ref_otus/wf_bdiv_40000/otu_table_mc2_w_tax_even40000_AP4only_HF_PMA.biom.txt") or die "Can't open biom file\n";
#open (biomout100, ">/output/open_ref_otus/wf_bdiv_40000/otu_table_mc2_w_tax_even40000_AP4only_LF_PMA_25percent_100greater.biom.txt") or die "Can't create biom file\n";
open (biomout100, ">/output/open_ref_otus/wf_bdiv_40000/otu_table_mc2_w_tax_even40000_AP4only_HF_PMA_25percent_100greater.biom.txt") or die "Can't create biom file\n";


while(my $line = <biomfile>) {
	chomp $line;
	#print $line,"\n";
	if ($line =~ "#" and $line =~"OTU") {
		print biomout100 $line, "\n";
		my @line = split("\t", $line);
		#print @line[0], "\t", @line[(scalar(@line)-2)], "\n";
	}
	else {
		my @line = split("\t", $line);
		#print @line, "\n";
		#$sumcolumn = scalar(@line)-1;
		#print $sumcolumn;
		my $samplepresencesum = 0;
		my $hundredgreater = "F";
		my $linesum = 0;
		for ($column= 1; $column<scalar(@line)-2; $column++) {
			$linesum = @line[$column] +$linesum;
			if (@line[$column] > 0) {
				$samplepresencesum++;
			}
			if (@line[$column] > 100) {
				$hundredgreater = "T";
			}
		}
		print @line[0], "\t", @line[(scalar(@line)-2)], "\n";
		#print $linesum, "\n";
		if ($samplepresencesum>($sumcolumn-2)/4 and $hundredgreater eq "T"){
			print biomout100 $line, "\n";  
		}
	}
}

close(biomfile);
close(biomout100);

