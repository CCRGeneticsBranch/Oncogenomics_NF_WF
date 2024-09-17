#!/bin/bash

set +e

allele=`cut -f1 $1 |grep -v Allele|tr '\n' ',' |sed -e 's/,$//g'`


pvacseq run --iedb-install-directory /opt/iedb -e 8,9,10,11,12,13,14 -t $5  --fasta-size=400 $2 $3 ${allele} {NNalign,NetMHC,NetMHCIIpan,NetMHCcons,NetMHCpan,PickPocket,SMM,SMMPMBEC,SMMalign} $4    2>$4/$3_err


is_empty=$(grep "TSV file is empty" $4/$3_err|wc -l)

if [ "$is_empty" -eq "1" ]
then
exit 0
fi
