use warnings;
open F1, "strain.txt";
open F2, "/home/mlk442/FIND/MicroMGx_Pilot/amplicons/429_16s.fasta";
open FH, ">16s_strain.fasta";
while ($line=<F1>){
	chomp($line);
	$hash{$line}=0;
}
#read= ">";
$flag =0;
$count =0;
while ($fa = <F2>){
	chomp($fa);
	if ($fa =~ /^>\w+_(\w+)/){
		$strain = $1;
		if (exists $hash{$strain}){
			print FH $fa, "\n";
			$count ++;
			$flag =1;
		}
		else{
			$flag =0;
		}
	}
	else{
		if ($flag == 1){
			print FH $fa, "\n";
		}
	} 
}
print $count, "\n";
