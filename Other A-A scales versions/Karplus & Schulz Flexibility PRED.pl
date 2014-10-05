############################################################
#Este aplicativo faz parte do pacote 'TRSampaio Prediction'#
# 														   #
# 						VERSÃO BETA						   #
#Predição de Epitopos Lineares de Células-B, Flexibility   #
#														   #
#Desenvolvido por Tiago Roberti Sampaio					   #
#Graduando em Ciência da Computação na UFMT Brasil		   #
############################################################

####VARIAVEIS CONFIGURAVEIS DO SISTEMA####
my $max_pepsize=7; #Tamanho da janela de predição (WINDOW)
my $threshold=0; #Taxa de precisao para verificar se é epitopo
$maximum=0;
$minimum=0;
$average=0;
##########################################


print "Entre com o nome do arquivo a ser lido\n"; #Cabecalho
my $archive_name=<STDIN>; #Le o nome do arquivo
chomp $archive_name; #Retira quebra de linha do nome do arquivo

#Imprime mensagem de erro caso não consiga abrir o arquivo
unless(open(FILE,$archive_name))
{
	print "Erro ao abrir \"$archive_name\"\n\n";
	exit;
}

#Arquivo é armazenado em forma de array
my @SEQ=<FILE>;
close FILE;
#Converte os dados do array para string
my $SEQS=join('',@SEQ);
#Retira possiveis espaços em branco
$SEQS=~ s/\s//g;

########Arquivo de entrada limpo!#####
#Verifica se o arquivo de entrada é do formato fasta
check_fasta();

#Transforma a string limpa em array novamente (retirando todas as quebras de linha)
@SEQ=split('',$SEQS);


#Cálculos para efetuar a divisão da cadeia
$tam=@SEQ; #Armazena quantos elementos continha no arquivo
$num_compep=$tam/$max_pepsize; #Calcula quantas vezes a cadeia será dividida - Ex: 23/5 = 4.6
$num_compep= int($num_compep); #Retorna a parte inteira
$num_incompep=$tam % $max_pepsize; #Retorna quantos elementos sobrarão das cadeias divididas - Ex: 23%5 = 3

sub check_fasta
{
	#ALFABETO DE ENTRADA: ACDEFGHIKLMNPQRSTVWY
	$count=($SEQS=~tr/ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy//);
	if($count!=length $SEQS){ 
	print "Arquivo de entrada nao e uma proteina \n";
	exit;
	}
}
sub cria_combinacoes
{
	$tamseq=length $SEQS;
	for($i=1;$i<=$max_pepsize-1;$i++)
	{
		@SEQ[$tamseq+$i]="X";
	}
	#print "OLD STRING: $SEQS \n @SEQ \n";
	$SEQS=join('',@SEQ);
	#print "NEW STRING: $SEQS \n";
	for($i=0;$i<$tamseq;$i++)
	{
		@write[$i]=substr($SEQS,$i,$max_pepsize);
	}
	unless(open(NEWFILE,">$archive_name.comb"))
	{
		print "Não foi possivel abrir o arquivo $archive_name.comb para escrever!\n"
	}
	$size=@write;
	for($k=0;$k<$size;$k++)
	{
		#print "Gravando: @write[$k] \n";
		print NEWFILE "@write[$k]\n";
	}
	close (NEWFILE);
	print "$size sequencias geradas\n";

}
cria_combinacoes();
flex();


		
sub flex
{
	my %flex_scale=('R' =>2.79,
						'D' =>1.42,
						'E' =>1.60,
						'K' =>2.88,
						'S' =>3.0,
						'N' =>1.77,
						'Q' =>1.18,
						'G' =>1.86,
						'P' =>0.52,
						'T' =>1.18,
						'A' =>-1.27,
						'H' =>-0.82,
						'C' =>-1.09,
						'M' =>-1.84,
						'V' =>-1.75,
						'I' =>-2.89,
						'L' =>-2.29,
						'Y' =>-3.3,
						'F' =>-2.14,
						'W' =>-3.78,
						'X' =>0.0);
						
			
	unless(open(FILE,"$archive_name.comb"))
	{
		print "Não foi possivel abrir o arquivo \"$archive_name.out\" para escrever!\n"
	}
	unless(open(TEMP,">$archive_name.temp"))
	{
		print "Não foi possivel abrir o arquivo \"$archive_name.out\" para escrever!\n"
	}
	$out=0; #Index do vetor final que ira armazenar o resultado
	@arq_final; #Index = $out
	while($data=<FILE>)
	{
		chomp $data;
		#print "$data ";
		@data=split('',$data);
		#print "@data \n";
		$total=0;
		#@flex_out[0]=$data;
		printf TEMP "%s ",$data;
		for($i=0;$i<$max_pepsize;$i++)
		{
			#print $flex_scale{@data[$i]},"  ";
			if(($flex_scale{@data[$i]}>=10)||($flex_scale{@data[$i]}*-1>=10)){
			#printf "%+.2f ",$flex_scale{@data[$i]};
			#@flex_out[$i+1]="$flex_scale{@data[$i]} ";
			printf TEMP "%+.2f ",$flex_scale{@data[$i]};
			}else{
			#printf "%+.2f  ",$flex_scale{@data[$i]};
			#@flex_out[$i+1]="$flex_scale{@data[$i]}  ";
			printf TEMP "%+.2f  ",$flex_scale{@data[$i]};
			}
			#printf "%+.2f  ",$flex_scale{@data[$i]};
			$total+=$flex_scale{@data[$i]};
		}
		#print "parkeeeer \n @flex_out \n parkeeer";
		$total/=$max_pepsize;
		if($total>=$maximum){ $maximum=$total; }
		if($total<=$minimum){ $minimum=$total; }
		if($out<$max_pepsize-1){ $average+=$total; }
		#print "average agora é $average\n";
		@flex_out[$max_pepsize+1]="= $total";
		$out++;
		#print "= $total \n";		
		#printf "= %+.3f \n",$total;
		printf TEMP "= %+.3f \n",$total;			
		
		#@arq_final[$out]=join(" ",@flex_out);
		 
		
	}
	close (FILE);
	$avg=($average/($max_pepsize-1));
		
	close (TEMP);
	
	####Le arquivo .out e verifica se é epitopo
	if($threshold==0){ $threshold=$avg; }
			
		unless(open(TEMP,"$archive_name.temp"))
		{
		print "Não foi possivel abrir o arquivo \"$archive_name.temp\" para escrever!\n"
		}
		unless(open(FINALOUT,">$archive_name.out"))
		{
		print "Não foi possivel abrir o arquivo \"$archive_name.out\" para escrever!\n"
		}
		
		#print "Lendo arquivo out\n";
		$epi_count=0;
		@epitopes;
		$row_numer=0;
		while($data=<TEMP>)
		{
			chomp $data;			
			@arr=split(' ',$data);
			$len=@arr;			
			#print "out @arr ";
			#print "size=$size $row_number\n";
			if((@arr[$len-1]>=$threshold) && ($row_number<$size-($max_pepsize-1)))
			{	
				#print "size=$size $row_number\n";
				print "$data E\n";
				print FINALOUT "$data E\n";
				@epitopes[$epi_count]=@arr[0];
				$epi_count++;
				
			}else{
				print "$data -\n";
				print FINALOUT "$data -\n";
			}	
			$row_number++;	
		}
		printf "Average %.3f|Minimum %.3f|Maximum %.3f|Threshold %.3f\n",$avg,$minimum,$maximum,$threshold;
		printf FINALOUT "Average %.3f|Minimum %.3f|Maximum %.3f|Threshold %.3f\n",$avg,$minimum,$maximum,$threshold;
		##Apresenta possiveis epitopos
		print "\n----------------Possiveis Epitopos----------------\n\n";
		print FINALOUT "\n----------------Possiveis Epitopos----------------\n\n";
		for($i=0;$i<$epi_count;$i++)
		{
			print "		Epitope -> $epitopes[$i]\n";
			print FINALOUT "		Epitope -> $epitopes[$i]\n";
		}
		print "\n--------------------------------------------------\n";
		print FINALOUT "\n--------------------------------------------------\n";
		###
	

}


sub divide_cadeia
{
	$i=0;
	print "Numero de elementos da cadeia: $tam \n";
	while($i<$num_compep)
	{
		@write[$i]=substr($SEQS,$i*$max_pepsize,$max_pepsize); #Armazena a cadeia dividida
		$i++;
	}
	@write[$i]=substr($SEQS,($i*$max_pepsize),$num_incompep); #Armazena o que sobrou das cadeias divididas
	
	#Apresenta mensagem de erro caso arquivo de saida não seja encontrado
	unless(open(NEWFILE,">$archive_name.cadeias"))
	{
		print "Não foi possivel abrir o arquivo \"$arquivosaida\" para escrever!\n"
		
	}
	#Grava a cadeia dividida no arquivo de saida
	for($k=0;$k<=$i;$k++)
	{
		print "Gravando: @write[$k] \n";
		print NEWFILE "@write[$k]\n";
	}
	close (NEWFILE);
}	