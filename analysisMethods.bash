if [[ "$@" =~ "--debug" ]]; then
	set -ex
else
	set -e
fi


##################PARSE PARAMETERS#############################

cfileband=0
localband=0
statusband=0
realdataband=0
simdataband=0

function parsefastaFunction {
	echo 'BEGIN{FS="|"}
	{
	if($1~">"){
		split($1,array,">")
		$1=array[2]
		for (i=1;i<100;i++){
		band=0;
			if($i != ""){
					if($i==ID){
						ID=$(i+1);
						print ID
						exit 0
					}
			}
					
		}
	}
	}'
}

function groupingFunction {
	echo 'args<-commandArgs()
	file<-c(args[6])
	headr<-c(args[7])
	if(headr=="T"){
			df<-read.csv(file, header = T, check.names = F)
			newdf<-aggregate(. ~ Species, df, FUN = sum)

	}else{
			df<-read.csv(file, header = F)
			colnames(df)<-c("COL1","COL2")
			newdf<-aggregate(. ~ COL1, df, FUN = sum)
	}
	write.table(newdf,file,row.names = F,quote = F)' > grp.R
}


for i in "$@"
do
	case $i in
	"--cfile")
		cfileband=1
	;;
	"--simdata")
		simdataband=1
	;;
	"--realdata")
		realdataband=1
	;;
	"--help")
		echo "Usage: bash analysisMethods.bash --cfile [config file] --simdata [simulation table csv format (',' separated)] --realdata [mconf from metasim]"
		echo -e"\nOptions explain:"
		echo "--cfile configuration file"
		echo "--simdata csv with your simulation data obtain from Parse Module provide by SEPA, or your own csv file (e.g. pathoscope_table.csv)"
		echo "--realdata a file that contain the real values of reads distribution (.mprf of metasim, or some file with format ti-reads [or gi-reads])"
		echo "Methods Aviables: RRMSE (normalized relative RMS error), ROC, SIMOBS. Specify in the config file using the flag 'ANALYSISTYPE'"
		echo -e "For example: ANALYSISTYPE=ROC,RRMSE,SIMOBS\n"
		echo "For ROC_CURVE you must specify a TOTALGENOMES=[integer] flag"
		echo "See the README for more information"
		exit
	;;
	*)
		
		if [ $((cfileband)) -eq 1 ];then
			
			if ! [ -f $i ];then
				echo "$i file no exist"
				exit
			fi

			for parameter in $(awk '{print}' $i)
			do
				Pname=$(echo "$parameter" |awk -F"=" '{print $1}')		
				case $Pname in
					"TOTALGENOMES")
						TOTALGENOMES=$(echo "$parameter" |awk -F"=" '{print $2}' |sed "s/,/ /g")
					;;
					"ANALYSISTYPE")
						ANALYSISTYPE=$(echo "$parameter" |awk -F"=" '{print $2}' |sed "s/,/ /g")
					;;
				esac
			done
			statusband=$((statusband+1))
			cfileband=0
		fi
		
		if [ $((realdataband)) -eq 1 ];then
			if [ -f $i ];then
				statusband=$((statusband+1))
				REALDATAFILE=$i
				realdataband=0
			fi
				#####################	FETCH ID REAL DATA	########################
				echo "checking your real data file (remember the format [abundance gi/ti id])"
				fields=$(awk '{print NF;exit}' $REALDATAFILE)
				case $fields in
					"2")
						echo "file already formated"
					;;
					"3")
						echo "ti reads" > rtmp
						parsefastaFunction > parsefasta.awk
						while read line
						do
							ID=$(echo "$line" |awk '{print $2}')
							#parsed file
							case $ID in
								"ti")
									echo "$line" |awk '{print $3, $1}' >> rtmp				
								;;
								"gi")
									abu=$(echo "$line" |awk '{print $1}')
									gi=$(echo "$line" |awk '{print $3}')
									ti=""
									echo "fetching ti by gi: $gi"
									
									while [ "$ti" == "" ]
									do
										ti=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nuccore&db=taxonomy&id=$gi" |grep "<Id>"|tail -n1 |awk '{print $1}' |cut -d '>' -f 2 |cut -d '<' -f 1)
									done
									echo "$ti $abu" >> rtmp				
								;;
							esac
						done < <(grep "" $REALDATAFILE)
						rm parsefasta.awk
						REALDATAFILE="rtmp"

						echo "Done"
					;;
					*)
						echo "unrecognized file"
						exit
					;;
				esac

				echo "convertion gixti for real data file"

				#fetch ti for species name
				firstline=0
				while read line
				do
					if [ $((firstline)) -eq 0 ]; then
						firstline=1
						echo "name reads" > newrealtmp
					else
						ti=$(echo "$line" |awk '{print $1}')
						abu=$(echo "$line" |awk '{print $2}')
						echo "fetching lineage from ti: $ti"
						#fetch the ti by ncbi api
						nofetch=""
						while [ "$nofetch" == "" ] || [[ "$nofetch" =~ "Connection refused" ]]
						do
							if curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=$ti" > tmp.xml ;then
								touch tmp.xml
								nofetch=$(cat tmp.xml)
							else
								echo "curl error fetch, internet connection?, retrying"
							fi
						done
						name=$(awk 'BEGIN{FS="[<|>]"}{if($2=="ScientificName"){printf "%s\n", $3;exit}}' tmp.xml) #be careful with \n
						spctoawk=$(awk 'BEGIN{FS="[<|>]"}{if($2=="ScientificName"){printf "%s\n", $3;exit}}' tmp.xml |awk '{print $2}')
						lineage=$(awk -v emergencyname=$spctoawk 'BEGIN{FS="[<|>]";prev="";superk="";phy="";class="";order="";fam="";gen="";spc=""}{if($2=="ScientificName"){prev=$3}if($3=="superkingdom"){superk=prev}if($3=="phylum"){phy=prev}if($3=="class"){class=prev}if($3=="order"){order=prev}if($3=="family"){fam=prev}if($3=="genus"){gen=prev}if($3=="species"){spc=prev}}
						END{if(superk==""){printf "unknown,"}else{printf "%s,",superk};if(phy==""){printf "unknow,"}else{printf "%s,",phy}; if(class==""){printf "unknow,"}else{printf "%s,",class}; if(order==""){printf "unknow,"}else{printf "%s,",order}; if(fam==""){printf "unknow,"}else{printf "%s,",fam}; if(gen==""){printf "unknow,"}else{printf "%s,",gen}; if(spc==""){if(emergencyname==""){print "unknow,"}else{printf "%s,",emergencyname}}else{printf "%s,",spc}}' tmp.xml)
						cand=$(echo "$lineage" |awk '{if($0 ~ "Candidatus"){print "YES"}else{print "NO"}}')
						unknowNum=$(echo "$lineage" |grep -o 'unknow' | wc -l)

						if [ $((unknowNum)) -ge 6 ] || [ "$lineage" == "" ];then
							echo "error to take lineage from ti $ti from taxonomy NCBI, using nuccore NCBI"
							#emergency step
							nofetch=""
							while [ "$nofetch" == "" ] || [[ "$nofetch" =~ "Connection refused" ]]
							do
								if curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$ti" > tmp.txt ;then
									touch tmp.txt
									nofetch=$(cat tmp.txt)
								else
									echo "curl (emergency step), error fetch, internet connection?, retrying"
								fi
							done
							name=$(grep "taxname" tmp.txt |awk -F"\"" '{print $2}' |awk '{print $1, $2}')
							lineageLine=$(grep -n "lineage" tmp.txt |cut -d " " -f 1 |sed "s/://g" |awk '{print $1+1}')
							lineage=$(head -n $lineageLine tmp.txt |tail -n2 |awk '{if($NF!=","){printf "%s",$0}else{printf "%s\n",$0}}' |awk -F"\"" '{gsub("; ",",");print $2}')
							#lineage only get until genus level, so, we added the species
							lineage=$(echo $lineage,$name)
						fi

						if [ "$cand" == "YES" ]; then
							name=$(echo "$name" |sed "s/Candidatus //g")
							echo "$name $abu" >> newrealtmp
						else
							name=$(echo "$name" |awk '{gsub("\\[|\\]","");print $1, $2}')
							echo "$name $abu" >> newrealtmp
						fi
						rm -f tmp.xml tmp.txt
					fi
				done < <(grep "" $REALDATAFILE)
				
				awk '{if(NR>1){print "\""$1" "$2"\""","$3}}' newrealtmp >rtmp
				rm -f newrealtmp
				REALDATAFILE="rtmp"
				
				#Grouping same names
				groupingFunction
				Rscript grp.R rtmp F

		fi
		
		if [ $((simdataband)) -eq 1 ];then
			if [ -f $i ];then
				#this is csv file
				statusband=$((statusband+1))
				SIMDATAFILE=$i
				BACKUPS=$i
				ORIGINALSNAME=$(echo "$SIMDATAFILE" |rev |cut -d "/" -f 1 |rev)

				#if [ $((notiWorkband)) -eq 1 ]; then
					echo "metaphlan,constrains,sigma convertion working" #now convertions apply to all softwares, check parseMethods for know the changes
					colti=$(awk -F"," '{if(NR==1){for (i=1;i<=9;i++){if($i ~ "Name"){print i;exit}}}}' $SIMDATAFILE)		
				#else
				#	colti=$(awk -F"," '{if(NR==1){for (i=1;i<=9;i++){if($i ~ "ti"){print i;exit}}}}' $SIMDATAFILE)
				#fi

				if [ "$colti" == "" ];then
					echo "header 'ti' not found in $SIMDATAFILE, check your csv"
					exit
				else
					#make sim data only with ti numbers
					awk -v initial=$colti -F"," '{for (i=initial;i<=NF;i++){printf "%s ",$i}printf "\n"}' $SIMDATAFILE > stmp
					SIMDATAFILE="stmp"

				fi

				simdataband=0

			else
				echo "ERROR: $i doesn't exist"
				exit
			fi
		fi
	;;
	esac
done

function getR2Function {
	
	echo 'args <-commandArgs()
	myfile<-args[6]
	
	correlaciones<-read.table(myfile,header=FALSE)
	
	summary(lm(V1~V2, data=correlaciones))'

}
function RMSfunction {
	folder=$(pwd)
	echo "RRMSE function called in $folder"

	totalcol=$(awk '{if(NR>1){print NF;exit}}' $SIMDATAFILE)
	echo -e "Analysis,RRMSE,NAVGRE" > RMS.$ORIGINALSNAME
	for coli in $(seq 3 1 $totalcol)	#col 1 and 2 always be the name
	do
		awk -v coli=$coli '{if(NR>1){print $1, $2, $coli}else{print $(coli-1) > "htmp"}}' $SIMDATAFILE > name_reads_tmp
		#in the next line, we find the ti that match in simulation data file.
		suma=0
		sumprom=0
		firstline=0

		while read line 
		do
			if [ $((firstline)) -eq 0 ];then
				firstline=1
			else
				name1=$(echo $line |awk '{print $1}')
				name2=$(echo $line |awk '{print $2}')
				abu=$(echo $line |awk '{print $3}')
				backupsuma=$suma
				backupprom=$sumprom

				suma=$(awk -v realname1=$name1 -v realname2=$name2 -v abu=$abu -v suma=$suma 'function abs(v) {return v < 0 ? -v : v}{
																		if($1==realname1 && $2==realname2){
																			if(abu==0){
																				exit #to avoid the 0/0
																			}
																			re=abs(($3-abu)/abu);
																			print (suma+(re*re));										
																			exit
																		}
																	}' name_reads_tmp)

				sumprom=$(awk -v realname1=$name1 -v realname2=$name2 -v abu=$abu -v suma=$sumprom 'function abs(v) {return v < 0 ? -v : v}{
																		if($1==realname1 && $2==realname2){
																				if(abu==0){
																					exit
																				}
																				re=abs(($3-abu)/abu);
																				print (suma+re);																						
																				exit
																			}
																		}' name_reads_tmp)
				if [ "$suma" == "" ];then
					suma=$(echo "$backupsuma")
				fi
		
				if [ "$sumprom" == "" ];then
					sumprom=$(echo "$backupprom")
				fi

			fi
		done < <(grep "" $REALDATAFILE)	#mprf parsed file have 'ti abundance' format
		#exit

		awk -v suma=$suma 'END{print sqrt(suma/(NR-1))}' $REALDATAFILE > rrmsetmp
		awk -v suma=$sumprom 'END{print suma/(NR-1)}' $REALDATAFILE > avg

		paste -d "," rrmsetmp avg > 2rms
		paste -d "," htmp 2rms > rrmse
		cat rrmse
		cat RMS.$ORIGINALSNAME rrmse > rmsvalues
		mv rmsvalues RMS.$ORIGINALSNAME
	done

	rm -f name_reads_tmp htmp rrmsetmp rrmse avg 2rms

}
function ROCfunction {
	
	re='^[0-9]+$'
	#only number no spaces
	if ! [[ "$TOTALGENOMES" =~ $re ]] ; then
		echo "you need to specify TOTALGENOMES flag with a integer number to calculate roc curves"
	else
		folder=$(pwd)
		echo "ROCfunction called in $folder"
	
		echo "fpr,tpr" > ROCtmp.dat
		echo "file" > filerocname
		totalcol=$(awk '{if(NR>1){print NF;exit}}' $SIMDATAFILE)
		for coli in $(seq 3 1 $totalcol)	#col 1 and 2 always be name in metaphlan,constrains and kraken, we begin in reads cols >=2
		do
			awk -v coli=$coli '{if(NR>1){print $1, $2, $coli}else{print $(coli-1) > "htmp"}}' $SIMDATAFILE > name_reads_tmp
			filename=$(cat htmp)
			TP=0	#true positive
			TN=0	#true negative
			FP=0	#false positive
			FN=0	#false negative
				
			###########TRUE POSITIVE AND FALSE POSITIVE###########
			while read name1 name2 reads
			do
				linetir=$(awk -v name1=$name1 -v name2=$name2 '{if($1==name1 && $2==name2){print $1, $2, $3;exit}}' $REALDATAFILE)
	
				if [ "$linetir" != "" ]; then
					name1=$(echo "$linetir" |awk '{print $1}')
					name2=$(echo "$linetir" |awk '{print $2}')
					readr=$(echo "$linetir" |awk '{print $3}')
					resultado=$(echo "$reads $readr" |awk '{if($1 >= $2/2){print "1"}else{print "0"}}') #bash doesn't work with float
					
					if [ "$resultado" == "1" ]; then
						TP=$(echo  "$TP+1" |bc)
					else
						band=$(echo "$reads $readr" |awk '{if($1==0 && $2>0){print "fp";exit}}')
						if [ "$band" == "fp" ]; then
							FP=$(echo "$FP+1" |bc)
						fi
					fi
					
				else
					FP=$(echo "$FP+1" |bc)
				fi
				
			done < <(grep "" name_reads_tmp)	#mprf parsed file have 'ti abundance' format
			
			###########FALSE NEGATIVE AND TRUE NEGATIVE###########
	
			while read line
			do	
				name1=$(echo "$line" |awk '{print $1}')
				name2=$(echo "$line" |awk '{print $2}')
				names=$(awk -v name1=$name1 -v name2=$name2 '{if($1==name1 && $2==name2){print $1, $2;exit}}' name_reads_tmp)
				if [ "$names" == "" ]; then
					FN=$(echo "$FN+1" |bc)
				fi								
			done < <(grep "" $REALDATAFILE)	#mprf parsed file have 'ti abundance' format
			TN=$(wc -l name_reads_tmp |awk '{print $1}' |awk -v total=$TOTALGENOMES -v tn=$TN '{print tn+(total-$1)}')
	
			######################################################
			
			echo -e "\nTP: $TP FP:$FP FN: $FN TN:$TN"
			tpr=$(echo "$TP $FN" | awk '{print $1/($1+$2)}')
			fpr=$(echo "$FP $TN" | awk '{print ($1/($1+$2))}')
			echo "fpr: $fpr tpr: $tpr"
			echo "$filename" >> filerocname
			echo "$fpr,$tpr" >> ROCtmp.dat
		done
		
		paste -d "," filerocname ROCtmp.dat > ROC.$ORIGINALSNAME
	#	sed "2d" ROC.$ORIGINALSNAME > tmp
	#	rm ROC.$ORIGINALSNAME
	#	mv tmp ROC.$ORIGINALSNAME
		rm -f name_reads_tmp htmp filerocname ROCtmp.dat
	fi

}
function SIMOBSfunction {
	folder=$(pwd)
	echo "Simulation vs Observed function called in: $folder"

	totalcol=$(awk '{if(NR>1){print NF;exit}}' $SIMDATAFILE)
	awk '{if(NR>1){print $1" "$2","$3}else{print "OTU,REAL"}}' $REALDATAFILE > SIMOBS.$ORIGINALSNAME

	for coli in $(seq 3 1 $totalcol)	#col 1 and 2 always be the name we begin in reads cols >=3
	do
		awk -v coli=$coli '{if(NR>1){print $1, $2, $coli}else{print $(coli-1) > "htmp"}}' $SIMDATAFILE > name_reads_tmp

		firstline=0
		while read line
		do
			if [ $((firstline)) -eq 0 ];then
				firstline=1
			else
				name1=$(echo "$line" |awk '{print $1}')
				name2=$(echo "$line" |awk '{print $2}')
				readr=$(echo "$line" |awk '{print $3}')
				reads=$(awk -v name1=$name1 -v name2=$name2 '{if($1==name1 && $2==name2){print $3;exit}}' name_reads_tmp)
				if [ "$reads" == "" ]; then
					echo "0" >> simobstmp
				else
					echo "$reads" >> simobstmp
				fi
			fi
		done < <(grep "" $REALDATAFILE)	#mprf parsed file have 'ti abundance' format
		cat htmp simobstmp > header
		paste -d "," SIMOBS.$ORIGINALSNAME header > filetmp
		mv filetmp SIMOBS.$ORIGINALSNAME
		rm simobstmp header
	done
	rm -f htmp name_reads_tmp

}

if [ $((statusband)) -ge 3 ]; then

	if ! [[ "$@" =~ "--cfile" ]];then
		echo "* No --cfile was specified"
		exit
	fi

	for a in $ANALYSISTYPE
	do
		case $a in
			"R2")
				R2function
			;;
			"RRMSE")
				RMSfunction
	   		;;
	   		"SIMOBS")
				SIMOBSfunction
			;;
	   		"ROC")
	   			ROCfunction
	   		;;
	   		*)
	   			echo "no method aviable for $a"
	   			exit
	   		;;
		esac
	done
			
	rm -f stmp rtmp
else
	echo "Invalid or Missing Parameters, print --help to see the options"
	echo "Usage: bash analysisMethods.bash --cfile [config file] --simdata [table csv] --realdata [mconf from metasim maybe]"
	exit
fi
