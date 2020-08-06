use warnings;
open F1, "/home/mlk442/FIND/mlk442/correlation_table/taxonomy.txt";
open F2, "/home/mlk442/FIND/MicroMGx_Pilot/amplicons/429_16s.fasta";
open FH, ">16s_strain.fasta";
while (<F1>){
	chomp;
	@line = split/\t/,$_;
	if (scalar(@line) > 3){
		$hash{$line[0]} = join (" ", $line[0],$line[2],$line[3]);
	#	$hash{$line[0]} = $line[0]."_".$line[2]."_".$line[3];
	}
}
#print scalar(keys @hash);
#read= ">";
$flag =0;
$count =0;
while ($fa = <F2>){
	chomp($fa);
	if ($fa =~ /^>\w+_(\w+)/){
		$strain = $1;
		if (exists $hash{$strain}){
			print FH ">", $hash{$strain},"\n";
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
