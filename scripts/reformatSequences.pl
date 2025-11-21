$dataF="./data";

%valid=("HA"=>1,
	"PB1"=>1,
	"PB2"=>1,
	"PA"=>1,
	"NP"=>1,
	"MP"=>1,
	"NS"=>1,
	"NA"=>1,
	"M"=>1
);

foreach $segment (keys %valid)
{
	if (-e "$dataF/$segment.segm.fa")
	{
		system("rm $dataF/$segment.segm.fa");
	}
}

while(<>)
{
	if ($_=~/^>(.*)/)
	{
		chomp();
		$id=$1;
		($id,$segment)=(split(/\|/,$id))[1,2];
		unless ($valid{$segment})
		{
			print "$_\n";
		}
		$segment="M" if $segment eq "MP";
		open(OUT,">>$dataF/$segment.segm.fa");
		print OUT ">$id\n";
	}else{
		$l=uc $_;
		print OUT $l;
	}
	
}
