$defining=shift;
open(IN,$defining);
$header=<IN>;
$totG=0;

while(<IN>)
{
	chomp();
	($lin,$muts)=(split(/\t/))[9,10]; #rimettere 10
	@muts=(split(/\,/,$muts));
	$tot{$lin}++;
	foreach $m (@muts)
	{
		$tl{$lin}{$m}++ if $lin ne "NA";
		$tt{$m}++;
	}
	$totG++;

}

open(OUT,">$defining.frequentVariants.txt");
print OUT "mutation\tcount\n";
foreach $m (keys %tt)
{
	print OUT "$m\t$tt{$m}\n" if $tt{$m}/$totG>=0.01;
}


open(OUT,">$defining.definingVariants.txt");
print OUT "designation\tgenomicVariants\n";
foreach $lin (sort keys %tl)
{
	next unless $tot{$lin}>=10;
	next if $lin eq "NA";
	print OUT "$lin";
	$cmut="\t";	
	foreach $m (keys %{$tl{$lin}})
	{
		$cmut.="$m " if $tl{$lin}{$m}/$tot{$lin}>=0.5;
	}
	chop($cmut);
	print OUT "$cmut\n";
}
