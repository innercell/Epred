############################################################
#Este aplicativo faz parte do pacote 'TRSampaio Prediction'#
# 														   #
# 						VERSÃO BETA	(With Weight Position) #
#Predição de Epitopos Lineares de Células-B, Hydrophilicity#
#														   #
#Desenvolvido por Tiago Roberti Sampaio					   #
#Graduando em Ciência da Computação na UFMT Brasil		   #
############################################################

####VARIAVEIS CONFIGURAVEIS DO SISTEMA####
my $max_pepsize=5; #Tamanho da janela de predição (WINDOW)
my $threshold=3; #Taxa de precisao para verificar se é epitopo
$maior=0;
$menor=0;
##########################################

#Escala de Parker com o Peso de Posicao
#  A       R       N       D       C       Q       E       G       H       I       L       K       M       F       P       S       T       W      Y      V  
#0.115  -0.087   0.154   0.327  -0.702   0.178   0.243   0.467  -0.159  -0.492  -0.577   0.068  -0.517  -0.462   1.102   0.186   0.155   0.002 -0.176 -0.340
#0.131  -0.091   0.165   0.331  -0.672   0.194   0.248   0.451  -0.182  -0.491  -0.584   0.067  -0.508  -0.494   1.069   0.212   0.192  -0.006 -0.210 -0.343
#0.120  -0.088   0.183   0.350  -0.716   0.185   0.268   0.445  -0.138  -0.508  -0.598   0.080  -0.507  -0.509   1.056   0.209   0.194   0.112 -0.177 -0.368
#0.134  -0.079   0.191   0.368  -0.700   0.176   0.277   0.407  -0.126  -0.500  -0.578   0.095  -0.508  -0.499   1.006   0.215   0.175   0.023 -0.156 -0.369
#0.132  -0.097   0.174   0.375  -0.697   0.167   0.275   0.364  -0.132  -0.497  -0.578   0.093  -0.487  -0.477   1.035   0.213   0.195   0.068 -0.133 -0.362

%pos5=('R' =>-0.097,
		'D' =>0.375,
		'E' =>0.275,
		'K' =>0.093,
		'S' =>0.213,
		'N' =>0.174,
		'Q' =>0.167,
		'G' =>0.364,
		'P' =>1.035,
		'T' =>0.195,
		'A' =>0.132,
		'H' =>-0.132,
		'C' =>-0.697,
		'M' =>-0.487,
		'V' =>-0.362,
		'I' =>-0.497,
		'L' =>-0.578,
		'Y' =>-0.133,
		'F' =>-0.477,
		'W' =>0.068,
		'X' =>0.0);

  %pos4=('R' =>-0.079,
		'D' =>0.368,
		'E' =>0.277,
		'K' =>0.095,
		'S' =>0.215,
		'N' =>0.191,
		'Q' =>0.176,
		'G' =>0.407,
		'P' =>1.006,
		'T' =>0.175,
		'A' =>0.134,
		'H' =>-0.126,
		'C' =>-0.700,
		'M' =>-0.508,
		'V' =>-0.369,
		'I' =>-0.500,
		'L' =>-0.578,
		'Y' =>-0.156,
		'F' =>-0.499,
		'W' =>0.023,
		'X' =>0.0);
 
 %pos3=('R' =>-0.088,
		'D' =>0.350,
		'E' =>0.268,
		'K' =>0.080,
		'S' =>0.209,
		'N' =>0.183,
		'Q' =>0.185,
		'G' =>0.445,
		'P' =>1.056,
		'T' =>0.194,
		'A' =>0.120,
		'H' =>-0.138,
		'C' =>-0.716,
		'M' =>-0.507,
		'V' =>-0.368,
		'I' =>-0.508,
		'L' =>-0.598,
		'Y' =>-0.177,
		'F' =>-0.509,
		'W' =>0.112,
		'X' =>0.0);
						
 %pos2=('R' =>-0.091,
		'D' =>0.331,
		'E' =>0.248,
		'K' =>0.067,
		'S' =>0.212,
		'N' =>0.165,
		'Q' =>0.194,
		'G' =>0.451,
		'P' =>1.069,
		'T' =>0.192,
		'A' =>0.131,
		'H' =>-0.182,
		'C' =>-0.672,
		'M' =>-0.508,
		'V' =>-0.343,
		'I' =>-0.491,
		'L' =>-0.584,
		'Y' =>-0.210,
		'F' =>-0.494,
		'W' =>-0.006,
		'X' =>0.0);
 
 
 %pos1=('R' =>-0.087,
		'D' =>0.327,
		'E' =>0.243,
		'K' =>0.068,
		'S' =>0.186,
		'N' =>0.154,
		'Q' =>0.178,
		'G' =>0.467,
		'P' =>1.102,
		'T' =>0.155,
		'A' =>0.115,
		'H' =>-0.159,
		'C' =>-0.702,
		'M' =>-0.517,
		'V' =>-0.340,
		'I' =>-0.492,
		'L' =>-0.577,
		'Y' =>-0.176,
		'F' =>-0.492,
		'W' =>0.002,
		'X' =>0.0);				

@wparker_scale=(\%pos1,\%pos2,\%pos3,\%pos4,\%pos5);
		
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
		print "Gravando: @write[$k] \n";
		print NEWFILE "@write[$k]\n";
	}
	close (NEWFILE);
	print "$size sequencias geradas\n";

}
cria_combinacoes();
parker();

		
sub parker
{
						
		
	unless(open(FILE,"$archive_name.comb"))
	{
		print "Não foi possivel abrir o arquivo \"$archive_name.out\" para escrever!\n"
	}

	$out=0; #Index do vetor final que ira armazenar o resultado
	@arq_final; #Index = $out
	while($data=<FILE>)
	{
		chomp $data;
		print "$data ";
		@data=split('',$data);
		#print "@data \n";
		$total=0;
		@parker_out[0]=$data;
		for($i=0;$i<$max_pepsize;$i++)
		{
			#print $parker_scale{@data[$i]},"  ";
			if(($wparker_scale[$i]{@data[$i]}>=10)||($wparker_scale[%i]{@data[$i]}*-1>=10)){
			printf "%+.3f ",$wparker_scale[$i]{@data[$i]};
			}else{
			printf "%+.3f  ",$wparker_scale[$i]{@data[$i]};
			}
			#printf "%+.3f  ",$parker_scale{@data[$i]};
			@parker_out[$i+1]=$wparker_scale[$i]{@data[$i]};
			$total+=$wparker_scale[$i]{@data[$i]};
		}
		#$total/=$max_pepsize;
		@parker_out[$max_pepsize+1]=$total;
		$out++;
		if($total>=$maior){ $maior=$total; }
		if($total<=$menor){ $menor=$total; }
		#print "= $total ";		
		printf "= %+.3f ",$total;
		if($total>$threshold)
		{
			@parker_out[$max_pepsize+2]="E";
			print "E \n";
		}else{
			@parker_out[$max_pepsize+2]="-";
			print "- \n";
		}		
		@arq_final[$out]=join(" ",@parker_out);
		
	}
	close (FILE);
	#Escrever arquivo final
	unless(open(OUT,">$archive_name.out"))
	{
		print "Não foi possivel abrir o arquivo \"$archive_name.out\" para escrever!\n"
	}
	
	foreach $gravar(@arq_final)
	{
		print OUT "$gravar\n";
	}
	close (OUT);
	
	print "Maximun $maior|Minimun $menor \n";

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