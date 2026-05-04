#!/usr/bin/bash

if [ x"$1" == x ]; then

        echo "please specify a GTF file"

        exit 1

fi


if [ x"$2" == x ]; then

        echo "Searching for feature 'gene'"
        search="gene"

else
	echo "Searching for feature '${2}'"
	search=$2
fi

if [[ $1 =~ ".gz" ]]; then
    catcmd="zcat"
else
    catcmd="cat"
fi
echo "gene_id	gene_name	gene_type" > gene_desc.txt
dict=()
$catcmd $1 |\
awk -v feat=$search -F"\t" '{ \
	if ($3==feat) {split($9,array,";"); \
		dict["1"]=""; \
                dict["2"]=""; \
		dict["3"]="NONE"; \
		for (i in array) { \
			if (array[i]~"gene_id") dict["1"]=array[i]; \
			if (array[i]~"gene_name") dict["2"]=array[i]; \
			if (array[i]~"gene_type") dict["3"]=array[i]; \
			if (array[i]~"gene_biotype") dict["3"]=array[i]; \
			gsub(/gene_id\ /, "", dict["1"]);\
			gsub(/gene_name\ /, "", dict["2"]);\
			gsub(/gene_type\ /, "", dict["3"]);\
			gsub(/gene_biotype\ /, "", dict["3"]);\
		} \
		print dict[1]"\t"dict[2]"\t"dict[3];
	} \
}' | sed s/' '//g | sort | uniq  >> gene_desc.txt


### creating gene_desc_coordinates.txt
echo "gene_id	gene_name	gene_type	chr	start	end	strand" > gene_desc_coordinates.txt
$catcmd $1 | awk -v feat=$search '{ if ($3==feat) { gsub(";",""); gsub("\"",""); print $10"\t"$14"\t"$12"\t"$1"\t"$4"\t"$5"\t"$7} }' | sort | uniq  >> gene_desc_coordinates.txt


