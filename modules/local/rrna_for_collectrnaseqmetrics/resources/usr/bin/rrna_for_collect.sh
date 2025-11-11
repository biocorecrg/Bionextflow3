#!/usr/bin/bash

if [ x"$1" == x ]; then

        echo "please specify a faidx file"
        exit 1

fi


if [ x"$2" == x ]; then
        echo "please specify a GTF file"
        exit 1
fi

if [[ $2 =~ ".gz" ]]; then
    catcmd="zcat"
else
    catcmd="cat"
fi

awk 'BEGIN{OFS="\t"} {print "@SQ","SN:"$1,"LN:"$2}' $1 > rRNA.list

$catcmd $2 | grep -v "##"| \
awk -F"\t" '{ OFS="\t"; if ($3=="transcript" && $9 ~ "transcript_type \"rRNA") {
     split($9,array,";"); for (i in array) {
	 	   if (array[i]~"transcript_id") id=array[i]; 
		   gsub(/transcript_id /, "", id);
		   gsub(/"/, "", id);
		}
    print $1,$4,$5,$7,id 
    } 
}' | sort -k1V -k2n -k3n  >> rRNA.list

