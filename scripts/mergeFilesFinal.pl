$checkfile="./scripts/internalCheckFile.csv";
open(IN,$checkfile);
while(<IN>)
{
	chomp();
	($id,@segms)=(split());
	$hc=pop(@segms);
	$s=join('',@segms);
	$checkData{$id}=[$hc,$s];
}

$dataFolder="./data";

@addF=qw(HA.HaploCoV.newD_v2_100  NA.HaploCoV.newD_v2_100  NS.HaploCoV.newD_v2_100  PB1.HaploCoV.newD_v2_100 M.HaploCoV.newD_v2_100  NP.HaploCoV.newD_v2_100  PA.HaploCoV.newD_v2_100  PB2.HaploCoV.newD_v2_100);

foreach $add (@addF)
{
	$nameSeq=(split(/\./,$add))[0];
	#$nameSeq="N0" unless $nameSeq;
	open(AD,"$dataFolder/$add");
	<AD>;
	while(<AD>)
	{
		chomp();
		($id,$segm)=(split(/\t/))[0,9];
		$dataCompl{$id}{$nameSeq}=$segm;
	}
}

$metadata="MetadataTable.tsv";
open(IN,"$dataFolder/$metadata");
$headerM=<IN>;
@segments=qw(PB2 PB1 PA HA NP NA M NS);
print "Virus name\tGenotype\tClade\tLineage\tContinent\tCountry\tregion\tHost\tCollection date\tOffet collection\tSubmission date\tOffset submission\tPB2\tPB1\tPA\tHA\tNP\tNA\tM\tNS\tHC_comb\n";
while(<IN>)
{
	chomp();
	($name,$genotype,$clade,$lineage,$location,$host,$collection,$submission)=(split(/\t/));
	($continent,$country,$region)=(split(/\//,$location));
	$continent ="NA" if $continent eq "";
	$country = "NA" if $country eq "";
	$region = "NA" if $region eq "";
	$ofsC=diff_d($collection);
	$ofsS=diff_d($submission);
	$segString="";
	$outString="$name\t$genotype\t$clade\t$lineage\t$host\t$continent\t$country\t$region\t$collection\t$ofsC\t$submission\t$ofsS";
	foreach $segment (@segments)
	{
		$val=$dataCompl{$name}{$segment} ? $dataCompl{$name}{$segment} : "notAvail";
		$outString.="\t$val";
		$segString.="$val";
	}
	$HC_comb=$checkData{$name}[0];
	$segStringCompare=$checkData{$name}[1];
	$isEq= $seqString eq $seqStringCompare ? 1 : 0;
	$HC_comb= "NotConsistent" unless $isEq;
	$outString.="\t$HC_comb";
	print "$outString\n";
}

sub diff_d
{
        my $diff="NA";
        my $date=$_[0];
        my @vl=(split(/\-/,$date));
        if ($date eq "NA")
        {
                return($diff);
        }else{
                $diff=0;
                my $diffY=($vl[0]-2019)*365;
                my $diffM=(int($vl[1]-1)*30.42); #Months
                my $diffD=$vl[2]-1;
                $diff=$diffY+$diffM+$diffD;
                return(int($diff));
        }
}

