#!/bin/sh


####################################################################################################
#  strainSpecific-iPCR.sh
####################################################################################################
#
#
#  Created by Jorge Omar Yanez Cuna on 14/03/2017.
#     NKI, Amsterdam
#
#   This program receives as input the alignment output in bed file for the mate1 and mate2 of a pair-end sequence.
#   Each mate has to be aligned separately. It will identify the strain for each of the reads.
#
#   This script received as input:
#     1. File with the list of all SNPs from both strains
#     2. Alignments bed file from bwbble, bowtie or the program of your choice for Mate1 and Mate2
#     3. fastq.gz file of the reads (for Mate1 and Mate2)
#     4. Directory of the output folder
#
#
#  Parameters
#     -s      File with SNPs
#     -d      Alignment file for reads from Mate1 (R1)
#     -e      Fastq file with all reads from Mate1 (R1)
#     -l      Alignment file for reads from Mate2 (R2)
#     -m      Fastq file with all reads from Mate2 (R2)
#     -o      output folder
#     -n      Pretty name for the output files
#
#  Input
#
#
#
#  OutPut
#
#
#  Example
#
#
#
####################################################################################################

snpFile="-"
mateOneBed="-"
mateOneFastq="-"
mateTwoBed="-"
mateTwoFastq="-"
outPutFolder=outFolder_strainSpecific-iPCR
gatcFile="-"
extansionBinArray=10000
shortFinalName="NA"

gapFile="/home/NFS/users/j.yanez/data/genomes/mm10/gaps_mm10.bed"
chrSizeFile="/home/NFS/users/j.yanez/data/genomes/mm10/mm10.chrom.sizes"

while getopts s:d:e:l:m:o:n: g
  do   case "$g" in
  s)  snpFile=$OPTARG;;
  d)  mateOneBed=$OPTARG;;
  e)  mateOneFastq=$OPTARG;;
  l)  mateTwoBed=$OPTARG;;
  m)  mateTwoFastq=$OPTARG;;
  o)  outPutFolder=$OPTARG;;
  n)  shortFinalName=$OPTARG;;
  esac
done

#####
# Checking parameters
#####

if [ ! -s $snpFile ]; then
 echo "-- Error: SNP file: ${gatcFile} do not exist. Please check parameter -s"
 exit
fi

if [ ! -s $mateOneBed -o ! -s $mateOneFastq ]; then
 echo "-- Error: Files for mate1 ${mateOneBed} and/or $mateOneFastq does not exist. Please check paramter -d and - e"
 exit
fi

if [ ! -s $mateTwoBed -o ! -s $mateTwoFastq ]; then
 echo "-- Error: Files for mate2 ${mateTwoBed} and/or $mateTwoFastq does not exist. Please check paramter -l and - m"
 exit
fi

#####
# Preparing variables and doing other stuff in order to start the analysis
#####

# Creating the output folder
mkdir -p $outPutFolder

# Obtaining base names for the output file
mateOneBedName=$(basename $mateOneBed | awk '{gsub(".bed","",$1); print $1}')
mateTwoBedName=$(basename $mateTwoBed | awk '{gsub(".bed","",$1); print $1}')

if [ $shortFinalName = "NA" ]; then
  shortFinalName=${mateOneBedName}_${mateTwoBedName}
fi

echo "${mateOneBedName}\t${mateTwoBedName}\t${shortFinalName}"

# unzipping fastq file
gunzip -c $mateOneFastq > ${outPutFolder}/${mateOneBedName}.fastq
gunzip -c $mateTwoFastq > ${outPutFolder}/${mateTwoBedName}.fastq

echo "-- 1. Identifying reads that overlap SNPs"

##############################
# Identifying reads that overlap SNPs
##############################

##########
# Identify the reads that mapped to
# 	129S1_SvImJ	CAST_EiJ
# Table with the SNP Analysis.
# Note that the bowtie are in 0-based coordinates and the SNP are in 1-based coordinates
# OutPut is in 1-based coordinates.
#   Format:
##########
 for line in "${mateOneBed}>${mateOneBedName}" "${mateTwoBed}>${mateTwoBedName}"; do
  alignmentFile=$(echo $line | awk '{split($0,arr,">"); print arr[1]}')
  shortName=$(echo $line | awk '{split($0,arr,">"); print arr[2]}')
  for pair in R1 R2; do
   cat $snpFile | awk -vOFS='\t' '{print $1,$2,$3,$4">"$5">"$6">"$7">"$8">"$9">"$10">"$11}' | \
    bedtools intersect -wo -a $alignmentFile -b - | \
     awk -vFile=${outPutFolder}/${shortName}.fastq -vOFS='\t' 'BEGIN{
      # Reading each of the reads
      while((getline<File)>0){ cRow++;
       if(cRow%4==1){ gsub("@","",$1); SeqID=$1;
       }else{ if(cRow%4==2){ hSeq[SeqID]=$1; } }
      }
     }{ split($4,aID,">"); readID=aID[1]; read=hSeq[readID]; snpStart=$8-$2;
      vGeno=$10; gsub("A/A","A",vGeno); gsub("C/C","C",vGeno); gsub("G/G","G",vGeno); gsub("T/T","T",vGeno); split(vGeno,aGeno,">"); strain="-";
      # Obtaining the sequence and if it is "-" reverse complement it
      if($6=="+"){ read=read;
      }else if($6=="-"){ tmpRead=read; read="";
       # Reversing the sequence
       gsub("A","t",tmpRead); gsub("T","a",tmpRead); gsub("G","c",tmpRead); gsub("C","g",tmpRead); for(i=length(tmpRead);i!=0;i--){ read=toupper(read substr(tmpRead,i,1))}
      }
      queryBP=substr(read,snpStart,1);
      # Checking the Genotypes
      if(aGeno[4]==queryBP){ if(strain=="-"){ strain="129S1"; }else{ strain=strain"_129S1" }}
      if(aGeno[7]==queryBP){ if(strain=="-"){ strain="CAST"; }else{ strain=strain"_CAST" }}
      print $1,$2,$3,readID,$5,$6,$7,$8,$9,vGeno,hSeq[readID],substr(read,1,snpStart-1)"-"queryBP"-"substr(read,snpStart+1),queryBP,strain;
     }' | gzip - > ${outPutFolder}/analysis_${shortName}_${pair}.txt.gz
  done
  rm ${outPutFolder}/${shortName}.fastq
 done

##########
# Assigning reads to strain & generating a bedpe file
#   Input Format:   chr, start, end, seqID, score, strand, Tot SNPs 129S1, Tot SNPs CAST, Tot SNPs 129 or CAST, Tot SNPs undefined;
#   OutPut Format:  ChrR1, StartR1, EndR1, ChrR2, StartR2, EndR2, SeqID, Barcode, Strand1, Strand2, Strain
##########
 for line in "${mateOneBed}>${mateOneBedName}" "${mateTwoBed}>${mateTwoBedName}"; do
  alignmentFile=$(echo $line | awk '{split($0,arr,">"); print arr[1]}')
  shortName=$(echo $line | awk '{split($0,arr,">"); print arr[2]}')
  for pair in R1 R2; do
   #   OutPut Format:  ReadID, Barcode, chr, start, end, score, strand, strain
   gunzip -c ${outPutFolder}/analysis_${shortName}_${pair}.txt.gz | \
    awk -vOFS='\t' '{ k=$1">"$2">"$3">"$4">"$5">"$6;
      if(!(k in hDid)){ c++; hReadIdx[c]=k; hDid[k]=1;
      } hReadStrain[k"_"$14]++;
     }END{ for(i=1;i<=c;i++){ k=hReadIdx[i]; tmpStrain=hReadStrain[k]; split(k,aPos,">"); print aPos[1],aPos[2],aPos[3],aPos[4],aPos[5],aPos[6],hReadStrain[k"_129S1"]+0,hReadStrain[k"_CAST"]+0,hReadStrain[k"_129S1_CAST"]+0,hReadStrain[k"_-"]+0; }
     }' | gzip - > ${outPutFolder}/assigningStrains_${shortFinalName}_${pair}.txt.gz
   # Separating the reads per strain
   gunzip -c ${outPutFolder}/assigningStrains_${shortFinalName}_${pair}.txt.gz | awk -vOFS='\t' '{ strain="unknown"
      if($7>0 && $8==0 && $7>=$10){ strain="129S1"; }
      if($7==0 && $8>0 && $8>=$10){ strain="CAST"; }
      # Additionally I want to add those reads where more than 2 times is assigned to one of the strains
      # Column 7 = Total of SNPs for 129S1; Col. 8= Tot SNPs for CAST; Col. 9= Tot SNPs undistinguishable bewteen 129S1 & CAST; Col. 10= Tot of SNPs undefined
      if( ($7>0 && $8>0) && ( ($8<=(($7/2)) && $8>=$10) ) ){ strain="129S1"; }
      if( ($7>0 && $8>0) && ( ($7<=(($8/2)) && $7>=$10) ) ){ strain="CAST"; }
      split($4,arr,"_"); print arr[1],arr[2],$1,$2,$3,$5,$6,strain
    }' > ${outPutFolder}/tmp_${shortFinalName}_${pair}_strainSpecific.txt
   # Reads that were not able to discriminate between SNPs
   cat $alignmentFile | awk -vFile=${outPutFolder}/tmp_${shortFinalName}_${pair}_strainSpecific.txt -vOFS='\t' 'BEGIN{
    while((getline<File)>0){ hDid[$1]=1; print $0
    }
   }{ split($4,arr,"_");
    if(!(arr[1] in hDid)){ print arr[1],arr[2],$1,$2,$3,$5,$6,"unknown" }
   }' | sort -k1,1 -k2,2 -k3,3g > ${outPutFolder}/strainSpecific_${shortFinalName}_${pair}.txt
   rm ${outPutFolder}/tmp_${shortFinalName}_${pair}_strainSpecific.txt ${outPutFolder}/analysis_${shortName}_${pair}.txt.gz
   rm ${outPutFolder}/assigningStrains_${shortFinalName}_${pair}.txt.gz
  done
 done

echo "-- 2. Generating a bedPE file by assembling the pairs R1, R2"

##########
#   Generating a bedPE file by assembling the pairs R1, R2
#   OutPut Format:  ChrR1, StartR1, EndR1, ChrR2, StartR2, EndR2, SeqID, Barcode, Strand1, Strand2, Strain
##########
 cat ${outPutFolder}/strainSpecific_${shortFinalName}_R1.txt | awk -vFile=${outPutFolder}/strainSpecific_${shortFinalName}_R2.txt -vOFS='\t' 'BEGIN{
    while((getline<File)>0){hR2[$1]=$0;}
   }{if($1 in hR2){ split(hR2[$1],arr,"\t");
      if($8==arr[8] || (arr[8]=="unknown" && $8!="unknown")){finalStrain=$8;}else if($8=="unknown" && arr[8]!="unknown"){ finalStrain=arr[8];}else{finalStrain=$8":"arr[8]}
      print $3,$4,$5,arr[3],arr[4],arr[5],$1,$2,$7,arr[7],finalStrain;
    }else{ print $3,$4,$5,".",".",".",$1,$2,$7,".",$8;
    }
   }' > ${outPutFolder}/tmp_strainSpecific_${shortFinalName}.bedpe
  cat ${outPutFolder}/strainSpecific_${shortFinalName}_R2.txt | awk -vFile=${outPutFolder}/tmp_strainSpecific_${shortFinalName}.bedpe -vOFS='\t' 'BEGIN{
     while((getline<File)>0){hDid[$7]==1; print $0;}
   }{
     if(!($1 in hDid)){ print ".",".",".",$3,$4,$5,$1,$2,".",$7,$8;}
   }' > ${outPutFolder}/strainSpecific_${shortFinalName}.bedpe
  rm ${outPutFolder}/strainSpecific_${shortFinalName}_R1.txt ${outPutFolder}/strainSpecific_${shortFinalName}_R2.txt
  rm ${outPutFolder}/tmp_strainSpecific_${shortFinalName}.bedpe


##########
# Generating a bedpe file with the unique list of integrations
#   Note. This is based on the script from the 2016.07.01 at the HybridMouseESC-V2.sh
##########
  # unique list of integrations for R1
  cat ${outPutFolder}/strainSpecific_${shortFinalName}.bedpe | awk -vOFS='\t' '{if($1!="." && $2!="."){print $1,$2,$3,$8}}' | sort -k4,4 -k1,1 -k2,2g | \
   awk -vOFS='\t' 'BEGIN{ fStart=1;
    # This time changing the barcode (line 4) should be considered as to change chromosome
    }{ k=$1"_"$2"_"$3"_"$4; fNotAdd=0; if(chr!=$1 || barcode!=$4){ if(fStart==0){ print chr,start,end,barcode,pos; } chr=$1; start=$2; end=$3; barcode=$4; fStart=0; pos=k;
    }else{
     if($2<=end){ split(pos,arr,";"); for(i=1;i<=length(arr);i++){ if(k==arr[i]){ fNotAdd=1; } } if(fNotAdd==0){ pos=pos";"k } if($3>end){ end=$3; }
     }else{ print chr,start,end,barcode,pos; chr=$1; start=$2; end=$3; barcode=$4; pos=k; } }
   }END{ print chr,start,end,barcode,pos; }' > ${outPutFolder}/uniqReadList_${shortFinalName}_R1.txt
  # unique list of integrations for R2
  cat ${outPutFolder}/strainSpecific_${shortFinalName}.bedpe | awk -vOFS='\t' '{if($4!="." && $5!="."){print $4,$5,$6,$8}}' | sort -k4,4 -k1,1 -k2,2g | \
   awk -vOFS='\t' 'BEGIN{ fStart=1;
    # This time changing the barcode (line 4) should be considered as to change chromosome
    }{ k=$1"_"$2"_"$3"_"$4; fNotAdd=0; if(chr!=$1 || barcode!=$4){ if(fStart==0){ print chr,start,end,barcode,pos; } chr=$1; start=$2; end=$3; barcode=$4; fStart=0; pos=k;
    }else{
     if($2<=end){ split(pos,arr,";"); for(i=1;i<=length(arr);i++){ if(k==arr[i]){ fNotAdd=1; } } if(fNotAdd==0){ pos=pos";"k } if($3>end){ end=$3; }
     }else{ print chr,start,end,barcode,pos; chr=$1; start=$2; end=$3; barcode=$4; pos=k; } }
   }END{ print chr,start,end,barcode,pos; }' > ${outPutFolder}/uniqReadList_${shortFinalName}_R2.txt
  # With the list from above, I updated the positions, make a bedpe file with the unique reads
  cat ${outPutFolder}/strainSpecific_${shortFinalName}.bedpe | \
   awk -vFileR1=${outPutFolder}/uniqReadList_${shortFinalName}_R1.txt -vFileR2=${outPutFolder}/uniqReadList_${shortFinalName}_R2.txt -vOFS='\t' 'BEGIN{
     while((getline<FileR1)>0){ split($5,arr,";"); for(i=1;i<=length(arr);i++){ hNewPosR1[arr[i]]=$1"\t"$2"\t"$3; }
     }
     while((getline<FileR2)>0){ split($5,arr,";"); for(i=1;i<=length(arr);i++){ hNewPosR2[arr[i]]=$1"\t"$2"\t"$3; }
     }
    }{ if($1!="."){ k=$1"_"$2"_"$3"_"$8; newPosMateOne=hNewPosR1[k]; }else{ newPosMateOne=$1"\t"$2"\t"$3 }
       if($4!="."){ k=$4"_"$5"_"$6"_"$8; newPosMateTwo=hNewPosR2[k]; }else{ newPosMateTwo=$4"\t"$5"\t"$6 }
       print newPosMateOne,newPosMateTwo,$7,$8,$9,$10,$11
    }' | sort -k8,8 -k1,1 -k2,2g -k4,4 -k5,5g > ${outPutFolder}/tmpUnique_strainSpecific_${shortFinalName}.bedpe
    cat ${outPutFolder}/tmpUnique_strainSpecific_${shortFinalName}.bedpe | awk -vOFS='\t' '{
      #printing unique integrations
      k=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$8"\t"$9"\t"$10"\t"$11; if(!(k in hDid)){ hDid[k]=1; cIntegration++; hInt[cIntegration]=k } hCountRead[k]++;
    }END{ for(i=1;i<=cIntegration;i++){ k=hInt[i]; split(k,arr,"\t"); print arr[1],arr[2],arr[3],arr[4],arr[5],arr[6],arr[7],arr[8],arr[9],hCountRead[k],arr[10]; }
    }' | sort -k1,1 -k2,2 -k4,4 -k5,5g > ${outPutFolder}/unique_strainSpecific_${shortFinalName}.bedpe

rm ${outPutFolder}/uniqReadList_${shortFinalName}_R1.txt ${outPutFolder}/uniqReadList_${shortFinalName}_R2.txt ${outPutFolder}/tmpUnique_strainSpecific_${shortFinalName}.bedpe








