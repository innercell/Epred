#|===================================================================================|	
#|                __________  _____                             _                    |	
#|               /_  __/ __ \/ ___/____ _____ ___  ____  ____ _(_)___                |	
#|                / / / /_/ /\__ \/ __ `/ __ `__ \/ __ \/ __ `/ / __ \               |	
#|               / / / _, _/___/ / /_/ / / / / / / /_/ / /_/ / / /_/ /               |	
#|              /_/ /_/ |_|/____/\__,_/_/ /_/ /_/ .___/\__,_/_/\____/                |	
#|                                             /_/                                   |	
#|                            B-HidroPred TRSampaio v1.0.5                           |	
#|        Prediction of B-cell Linear Epitopes with Parker Hidrophicility Scale      |
#| ================================================================================= |
#|   Developed by Tiago Roberti Sampaio | Federal University of Mato Grosso, Brazil  |  
#| ================================================================================= |
#|                trsampaio@outlook.com | tiagosampaio@outlook.com                   |
#| ================================================================================= |
#| v1.0.5 : [ADD] Automatic multiple windows prediction.	                         |
#| v1.0.5 : [ADD] New major function Run(), controls the execution order.            |
#| v1.0.5 : [ADD] Log file created.							                         |
#| v1.0.4 : [FIX] Unused variables and functions were excluded.                      |
#| v1.0.4 : [FIX] Translated and Restructured code for better understanding.         |
#| v1.0.3 : [ADD] Epigroup output filename changed to w"window_value"-epigroup.      |
#| v1.0.3 : [ADD] Storage in array the Start and End of the epitopes grouped.        |
#| v1.0.2 : [ADD] New option to save the epitopes grouped in file (EPIGROUP_FILE).   |
#| v1.0.2 : [ADD] Minimum window is set to 3.                                        |
#| v1.0.2 : [ADD] Search window range is set to 3 (i-3,i+3).                         |
#| v1.0.2 : [ADD] Epitopes with sequences in common are grouped and displayed.       |
#| v1.0.1 : [ADD] If the threshold value equals zero, then threshold equals average. |
#| v1.0.1 : [ADD] Auxiliary TEMP file to keep all sequences without epitope flag.    |
#| v1.0.1 : [ADD] Epitope flag storage.                                              |
#| v1.0.1 : [FIX] Output file wrong spaces.                                          |
#| v1.0.1 : [ADD] Average, minimum and maximum epitope values.                       |
#| v1.0.0 : Version Developed.                                                       |
#|===================================================================================|

####MODIFIABLE VARIABLES####
$window_size=3; #Size of the search window  (WINDOW) 
$threshold=0; #threshold tax (0 - calculate the average of all sequences)
$test_all_windows=1; #Test all windows values, 3 to 7 | 0 = false | 1 = true
$max_window_search=15;
####SYSTEM USE VARIABLES#### 
$range=3; #number of extra A-A to evade false negatives, 3 extra in the begin of chain and 3 extra in the end of chain
$maximum=0;
$minimum=0;
$average=0;
$group=0; #number of epigroups

#obs. start + (range+1) = A-A number on sequence


##########################################
open(LOG_FILE,">LOGFILE.log");
		print LOG_FILE "LOG FILE \n";


if($test_all_windows==0) #Execute one time
{
	
	run();
}else{
	for($w=3;$w<=$max_window_search;$w++)
	{
		$window_size=$w;
		print"\n\n++++++++++++WINDOW SIZE = $window_size++++++++++\n";
		print LOG_FILE "\n\n++++++++++++WINDOW SIZE = $window_size++++++++++\n";
		run();
		read_and_write_epigroups();
		
	}
}

sub run
{
	first_step();
	make_combinations();
	parker();
	write_result_file();
	group_epitopes();	
}

sub first_step
{
		if($archive_name eq "")
		{
		#Read the source file (*.fsa,*.seq,etc)
		print "Enter the name of the file to be read:\n"; 
		$archive_name=<STDIN>; #Capture the archive name
		chomp $archive_name; #Remove the \n of the end of the file
		}
		
		
		
		#If the archive can't open display a error message
		unless(open(SOURCE_FILE,$archive_name))
		{
			print LOG_FILE "Error!\"$archive_name\" not found!\n\n";
			print "Error!\"$archive_name\" not found!\n\n";
			exit;
		}
		#Storage the FILE in a array
		@SF_content=<SOURCE_FILE>;
		close SOURCE_FILE;
		#Convert the data of the array @SF_content to string (Removes \n)
		$SF_content_string=join('',@SF_content);
		#Remove all white space
		$SF_content_string=~ s/\s//g;

		#Check if the source_file is valid
		check_fasta();

		#Transform the cleaned string to array again (The source file is prepared to be used)
		@SF_content=split('',$SF_content_string); #Now every positions is a A-A

		#Calculate the properties of the source file
		$tam=@SF_content; #How many elements has in the file
		$original_size_SF=$tam;
		#Number of the complete peptides
		$num_compep=$tam/$window_size; #Calcule how many time the chain will be splited - Ex: 23/5 = 4.6
		$num_compep= int($num_compep); #Return the integer part
		#Number of the incomplete peptides
		$num_incompep=$tam % $window_size; #Return how many elements will remain after the split_chain - Ex: 23%5 = 3

		
}

sub check_fasta
		{
			#INPUT ALPHABET: ACDEFGHIKLMNPQRSTVWYXZ
			$count=($SF_content_string=~tr/ACDEFGHIKLMNPQRSTVWYXZacdefghiklmnpqrstvwyxz//); #Count how many elements are in the input alphabet
			if($count!=length $SF_content_string){ #Compare with the chain size
			print "Input file isn't a protein!\n";
			exit;
			}
		} 
		
sub make_combinations
{
	#Get the size of the chain
	$size=length $SF_content_string;
	#If has incomplete peptides it'll need to be completed by 'X' to replace the NULL SPACE  
	#Ex: Window = 5 | Combination created = 'ATE  ' (put the X) = 'ATEXX'
	for($i=1;$i<=$window_size-1;$i++)
	{
		@SF_content[$size+$i]="X";
	}
	#Get the new source file content (with the 'X' values)
	$SF_content_string=join('',@SF_content);
	
	#Split the chain with the $window_size lenght, and storage into a array to be written later
	for($i=0;$i<$size;$i++)
	{
		@write[$i]=substr($SF_content_string,$i,$window_size); #Ex: Sequence start in 0 and ends in $window_size (5) -> 0~5, after, 1~5, 2~5, ...
	}
	
	#Open the file to be written
	open(W_COMB_FILE,">$archive_name.comb");
	
	$size=@write; #Get how many chains was created after the split combination
	for($k=0;$k<$size;$k++)
	{
		print W_COMB_FILE "@write[$k]\n"; #Write to combination file
	}
	close (W_COMB_FILE);
	print "Generated sequences: $size\n";
	print LOG_FILE "Generated sequences: $size\n";

}

		
sub parker
{
	#Declare the hidrophilicty scale values
	my %parker_scale=('R' =>4.2,
						'D' =>10.0,
						'E' =>7.8,
						'K' =>5.7,
						'S' =>6.5,
						'N' =>7.0,
						'Q' =>6.0,
						'G' =>5.7,
						'P' =>2.1,
						'T' =>5.2,
						'A' =>2.1,
						'H' =>2.1,
						'C' =>1.4,
						'M' =>-4.2,
						'V' =>-3.7,
						'I' =>-8.0,
						'L' =>-9.2,
						'Y' =>-1.9,
						'F' =>-9.2,
						'W' =>-10.00,
						'X' =>0.0,
						'Z' =>0.0);
						
	#If the archive can't open display a error message		
	unless(open(R_COMB_FILE,"$archive_name.comb"))
	{
		print "Error!\"$archive_name.comb\" not found!\n\n";
		exit;
	}
	#Open the file to be written
	open(W_TEMP_FILE,">$archive_name.temp"); #Will storage all the sequences and values without the EPITOPE FLAG
	
	$lines_count=0; #Index do vetor final que ira armazenar o resultado
	while($sequence=<R_COMB_FILE>) #Read line per line all the sequences value
	{
		chomp $sequence; #Removes the \n of the end of line
		@data=split('',$sequence); #Convert the sequence string into array - Every position will be a A-A
		$total=0; #Used to calculate the average of the sum of the all A-A scale values
		
		printf W_TEMP_FILE "%s ",$sequence; #Print the sequence 
		for($i=0;$i<$window_size;$i++) #Travels horizontally and replace the A-A letter to your respect value in the scale
		{
			#This if/else its only to a better view format
			if(($parker_scale{@data[$i]}>=10)||($parker_scale{@data[$i]}*-1>=10)){ #If value its bigger than 9, then use one space Ex: "10 "
					printf W_TEMP_FILE "%+.2f ",$parker_scale{@data[$i]};
			}else{ #If value its smaller than 10, then use two spaces Ex: "7  "
					printf W_TEMP_FILE "%+.2f  ",$parker_scale{@data[$i]};
			}
			$total+=$parker_scale{@data[$i]}; #Sum all the A-A scale values
		}
		$total/=$window_size; #Calculate the average of the sum of all the A-A scale values
		if($total>=$maximum){ $maximum=$total; } #Get the maximum value
		if($total<=$minimum){ $minimum=$total; } #Get the minimum value
		
		#Do you remember of the 'X' that i used to complete the incomplet peptides? Okay...
		#The number of the 'X' incremented in the chain is $window_size-1
		#i'll not count in the average the sequences that were filled with the 'X'
		if($lines_count<$original_size_SF-($window_size-1)){ $average+=$total;}
		$lines_count++; #With this i know what line i'm reading
		
		printf W_TEMP_FILE "= %+.3f \n",$total;	#Write in the file the average of the sum of the all A-A scale values		 		
	}
	close (R_COMB_FILE);
	
	$avg=($average/($original_size_SF-($window_size-1))); #Calculate the average of the individual sequence results (without the sequences were filled with 'X')
	close (W_TEMP_FILE);
}

sub write_result_file
{
		if($threshold==0){ $threshold=$avg; }
			
		unless(open(R_TEMP_FILE,"$archive_name.temp"))
		{
		print "Error!\"$archive_name.temp\" not found!\n"
		
		}
		open(W_RESULT_FILE,">$archive_name.out");		
		
		$epi_count=0; #Count the possible epitopes
		@epitopes; #Will storage all the possible epitopes
		$row_number=0; #Point to current line
		while($data=<R_TEMP_FILE>) #Read line per line
		{
			chomp $data; #Remove the \n of the end of the line			
			@arr=split(' ',$data); #Convert the line string into array
			$last=@arr; #Get the lenght of the array (Every position is a column of the line)			
			$last--; #Point to the last position of the array @arr
			
			#Check if its epitope, i.e., if the last position of the line (the average of the sum of all A-A scale values)>= threshold(fit)
			#And exclude all the sequences that were filled with 'X', automatically isn't a epitope.
			if((@arr[$last]>=$threshold) && ($row_number<$original_size_SF-($window_size-1))) 
			{	
				print "$data E\n"; #Visual Add to the line another column, the epitope flag with the result (POSITIVE='E' or NEGATIVE='-')
				print LOG_FILE "$data E\n"; #Visual Add to the line another column, the epitope flag with the result (POSITIVE='E' or NEGATIVE='-')
				print W_RESULT_FILE "$data E\n"; #File Add to the line another column, the epitope flag with the result (POSITIVE='E' or NEGATIVE='-')
				@epitopes[$epi_count]=@arr[0];
				$epi_count++;
				
			}else{
				print "$data -\n"; #Visual Add to the line another column, the epitope flag with the result (POSITIVE='E' or NEGATIVE='-')
				print LOG_FILE "$data -\n"; #Visual Add to the line another column, the epitope flag with the result (POSITIVE='E' or NEGATIVE='-')
				print W_RESULT_FILE "$data -\n"; #FILE Add to the line another column, the epitope flag with the result (POSITIVE='E' or NEGATIVE='-')
			}	
			$row_number++; #Point to the next line
		}
		close (R_TEMP_FILE);
		
		#Display and write the global average, minimum, maximum and threshold values
		printf "\nAverage %.3f|Minimum %.3f|Maximum %.3f|Threshold %.3f\n",$avg,$minimum,$maximum,$threshold;
		printf W_RESULT_FILE "\nAverage %.3f|Minimum %.3f|Maximum %.3f|Threshold %.3f\n",$avg,$minimum,$maximum,$threshold;
		printf LOG_FILE "\nAverage %.3f|Minimum %.3f|Maximum %.3f|Threshold %.3f\n",$avg,$minimum,$maximum,$threshold;
		
		###Display all the possible epitopes
		print "\n----------------Possible Epitopes----------------\n\n";
		print LOG_FILE "\n----------------Possible Epitopes----------------\n\n";
		print W_RESULT_FILE "\n----------------Possible Epitopes----------------\n\n";
		for($i=0;$i<$epi_count;$i++) #The same than a foreach, displays all @epitopes positions
		{
			print "		Epitope -> $epitopes[$i]\n";
			print LOG_FILE "		Epitope -> $epitopes[$i]\n";
			print W_RESULT_FILE "		Epitope -> $epitopes[$i]\n";
		}
		print "\n--------------------------------------------------\n";
		print LOG_FILE "\n--------------------------------------------------\n";
		print W_RESULT_FILE "\n--------------------------------------------------\n";
		###
		close (W_RESULT_FILE);
}


sub group_epitopes
{
	unless(open(R_RESULT_FILE,"$archive_name.out"))
	{
		print "Não foi possivel abrir o arquivo \"$archive_name.out\" para escrever!\n"
	}
	
	@content=<R_RESULT_FILE>; #Storage the file content into a array
	close (R_RESULT_FILE);
	$group=0; #number of the current group
	$tam=@content; #Get how many lines has into the file
	for($i=0;$i<$tam;$i++) #Read line per line
	{
		@line=split(' ',$content[$i]); #Transform the single string inside the array in a array (every position is a column of the line)
		$last=@line; #Get how many columns has in line
		$last--; #Get the last position of the line
		
		if($line[$last]eq"E") #The last position is the epitope flag
		{	
			$count_range=0; #Storage how many sequences is a false epitope flag
			if($i>=3)
			{			
				if($i-3<$start){$start=$i-3;} #Define the new start of the macro sequence
				if($i+3>$end){$end=$i+3;} #Define the new end of the macro sequence
			}else{ #if the current position is less than 3, then to prevent bugs
				$start=0; #Start is smaller than its possible
				$end=$i+3; 
			}
		}elsif($line[$last]eq"-"){
			$count_range++; #When count_range is equal 3 a macro sequence will be generated
		}
		if($window_size==3){$range=2;}
		if($count_range==$range)
		{
			$count_range=0; #Set 0 to discover the next macro sequence
			#if($end!=-1){print "************epitope grouped: start=$start|end=$end \n";}
			#Storage the start and the end of the macro sequence			
			@start_epigroup[$group]=$start; 
			@end_epigroup[$group]=$end;
			if($window_size==3)
			{ 
				$range=2; 
				if($start!=0){ $start--;}
				
			}else{ $range=3; }
			
				$for_start=$start+3;
				if($start==0)
				{ 
					$for_start=$start; 
				}
				@row=split(' ',$content[$for_start]); #Transform the current line into a array (Every position will be a column)					
				@epigroup[$group]=substr($row[0],0,$range);
				for($k=$for_start;$k<=$end;$k++)
				{
					@row=split(' ',$content[$k]); #Transform the current line into a array (Every position will be a column)					
					@letter=split('',$row[0]);#Transform the selected A-A sequence into a array (Every position will be a A-A)											
					$arrstr=join('',$epigroup[$group]); #Recover the sequence already storaged					
					if($window_size>3)
					{
						@epigroup[$group]=$arrstr.$letter[3]; #Concatenate the sequence already storaged if the new A-A letter	
					}elsif($window_size==3){
						@epigroup[$group]=$arrstr.$letter[2]; #Concatenate the sequence already storaged if the new A-A letter	
					}else{
						print "Search window is too small!\n";
						exit;
					}
					#printf("Add to epigroup n:$group <- $letter[3]\n");
					#printf("Epigroup n:$group = $epigroup[$group]\n");				
				}
							
			if($end!=-1){$group++;} #If the $end was changed, then point to the new group ID
			$end=-1; #Reset end value
			$start=9999; #Reset start value			 
		}	
	}
	
	open(W_EPIGROUP_FILE,">$archive_name.w$window_size-epigroup");	
	
	#Display and write all the epitopes grouped storaged
	for($i=0;$i<$group;$i++)
	{
		print W_EPIGROUP_FILE "EpiGroup n$i: $epigroup[$i] start $start_epigroup[$i] | end $end_epigroup[$i]\n";
		print "EpiGroup n$i: $epigroup[$i] start $start_epigroup[$i] | end $end_epigroup[$i]\n";
		print LOG_FILE "EpiGroup n$i: $epigroup[$i] start $start_epigroup[$i] | end $end_epigroup[$i]\n";
	}
	close (W_EPIGROUP_FILE);		
}


sub read_and_write_epigroups
{
	$s_index=0; #index of the array that will storage all sequences
	for($i=3;$i<=$max_window_search;$i++)
	{
			open(R_EPIGROUP_FILE,"$archive_name.w$i-epigroup");
			print "====FROM $archive_name.w$i-epigroup====\n";
			print LOG_FILE "====FROM $archive_name.w$i-epigroup====\n";
			#@stock="";
			while($data=<R_EPIGROUP_FILE>) #Read line per line
			{
				chomp $data; #Remove the \n of the end of the line			
				@arr=split(' ',$data); #Convert the line string into array
				@stock[$s_index]=$arr[2];
				print "$arr[2]\n";
				print LOG_FILE "$arr[2]\n";
				$s_index++;
			}
			#print "estoq= @stock\n";
			close (R_EPIGROUP_FILE);
	}
	
	
}
printf LOG_FILE "\nAverage %.3f|Minimum %.3f|Maximum %.3f|Threshold %.3f\n",$avg,$minimum,$maximum,$threshold;
printf "\nAverage %.3f|Minimum %.3f|Maximum %.3f|Threshold %.3f\n",$avg,$minimum,$maximum,$threshold;
close (LOG_FILE);
