@defFiles=<data/*HaploCoV_V2_5_100*>;
$num=0;

foreach $d (@defFiles)
{
	$segm=(split(/\./,$d))[0];
	$segm=~s/data\///;
	open(IN,$d);
	open(OUT,">data/$segm.phenetic_V2.csv");
	
	$header=<IN>;
	@labs=();
	%see=();
	%track=();
	while(<IN>)
	{
		chomp();
		($lab,@muts)=(split());
		push(@labs,$lab);
		foreach $mut (@muts)
		{
			$track{$mut}{$lab}++;
		}
	}
	print OUT " @labs\n";
	foreach $mut (sort keys %track)
	{
		print OUT "$mut";
		foreach $gr (@labs)
		{
			$val= $track{$mut}{$gr} ? 1 : 0;
			print OUT " $val";
		}
		print OUT "\n";
	}
}
