#!/usr/bin/env bash

##########
## Vars ##
##########
shopt -s extglob


##########
## Main ##
##########
## Make directories
rm -rf msa/
rm -rf hmms/
rm -rf models/
rm -rf ../HMMs
rm mmarc_model_annotations.tsv

mkdir msa
mkdir msa/{groupI,groupII,groupIII}

mkdir models
mkdir models/{groupI,groupII,groupIII}

mkdir hmms
mkdir hmms/{groupI,groupII,groupIII}

mkdir ../HMMs

./usearch -cluster_fast ../mmarc_groupI_seqs.fasta -id 0.8 -sort length -msaout msaI_ -uc ../cluster_info/groupI_clusters.uc
mv msaI_* msa/groupI/

./usearch -cluster_fast ../mmarc_groupII_seqs.fasta -id 0.8 -sort length -msaout msaII_ -uc ../cluster_info/groupII_clusters.uc
mv msaII_* msa/groupII/

./usearch -cluster_fast ../mmarc_groupIII_seqs.fasta -id 0.8 -sort length -msaout msaIII_ -uc ../cluster_info/groupIII_clusters.uc
mv msaIII_* msa/groupIII/

#./usearch -cluster_fast ../mmarc_groupIII_seqs.fasta -id 0.8 -sort length -msaout msa_blast_ -uc ../cluster_info/blast_augmented_clusters.uc
#mv msa_blast_* msa/blast_augmented/

./mmarc_annotate_models.py msa ../mmarc_sequence_annotations.csv models

(head -n 1 mmarc_model_annotations.tsv && tail -n +2 mmarc_model_annotations.tsv | sort -k1,1 -k2,2) > ../mmarc_model_annotations.tsv


for i in $( echo models/* ); do
	for j in $( echo ${i}/* ); do
		prefix=$( echo $i | sed 's/models/hmms/' )
		name=$( basename $j | sed 's/.fasta//' )
		hmmbuild -n $name --dna --cpu 7 ${prefix}/${name}.hmm $j
	done
done

for i in $( echo hmms/* ); do
	model_name=$( echo ${i} | sed 's/hmms\///' )
	cat $( echo ${i}/*) > ../HMMs/mmarc_${model_name}.hmm
		hmmpress ../HMMs/mmarc_${model_name}.hmm
done

exit 0

