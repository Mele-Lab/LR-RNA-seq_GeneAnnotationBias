#!usr/bin/perl
#
open(REM, $ARGV[0]);
open(IN, $ARGV[1]);
open(REMOVE, '>',$ARGV[2]);
open(KEEP, '>',$ARGV[3]);

while(<REM>)
{	chomp;
	if($_=~/^@/)
	{
		$id=substr($_,1);
	}
	else{$id=$_;}
	@idOnly=split("_",$id);
	#print "$idOnly[0]\n";
	$multiSplitIDs{$idOnly[0]}="";
#	print "$_\t$id\n";
}
#print %multiSplitIDs;

$line=0;
$removed=0;$kept=0;
while(<IN>)
{	chomp;
	if(($line % 4)==0)
	{	$match=0;
		@ID_info=split(" ",$_);
		$ID1=substr($_,1);
		#print "$ID_info\n";
		@currID=split(" ",$ID1);		
		if($currID[0]=~/_/)
		{	#print "$currID[0]\n";
			@currID=split("_",$currID[0]);
			#print "&&$currID[0]\t@currID\n";
		}
		chomp $currID[0];
		#print "#####$currID[0]\t$multiSplitIDs{$currID[0]}\n";
#
		if(exists($multiSplitIDs{$currID[0]}))
		{	$match=1;#print "MATCHES$_\n";
		}
	}#print "$match\n";
	if($match == 1)
	{	if(($line % 4)==0)
		{	print REMOVE "$ID_info[0]\n";
		}else{	print REMOVE "$_\n";
		}
		$removed++;
	}
	else
	{	if(($line % 4)==0)
		{	print KEEP "$ID_info[0]\n";
		}else{	print KEEP "$_\n";
		}
		$kept++;
	}
	$line++;

}
print "$removed removed and $kept lines kept\n";

