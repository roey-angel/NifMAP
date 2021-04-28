#!/bin/bash -
#title          :NifMAP.sh
#description    :NifMAP - a bioinformatics pipline for analysing nifH amplicon data
#author         :Roey Angel
#date           :20210309
#version        :1.2
#usage          :./nifmap.sh <inputDir>
#notes          : Unlike V1.0, this pipeline is now independent and performs quality-filtering and OTU generation.
#               See Angel, R., Nepel, M., Panhölzl, C., Schmidt, H., Herbold, C. W., Eichorst, S. A., et al. (2018). Evaluation of primers targeting the diazotroph functional gene and development of NifMAP – a bioinformatics pipeline for analyzing nifH amplicon data. Front Microbiol 9, 703. doi:10.3389/fmicb.2018.00703.
#dependencies   : 1. HMMER - http://hmmer.org/
#                 2. FrameBot - https://github.com/rdpstaff/Framebot
#                 3. seqmagick - https://fhcrc.github.io/seqmagick/
#                 4. CART - https://wwwzehr.pmc.ucsc.edu/CART_model_public/
#                 5. R - https://www.r-project.org/
#bash_version   :4.3.48(1)-release
#============================================================================

eval "$(conda shell.bash hook)" # this is needed to be able to use conda activate

# Define variables
inputDir="./Data" #$1 # input fasta file containing merged MiSeq reads
HOMEFOLDER=`pwd` # base library
FRAMEBOTPATH="/usr/local/lib/Bioinformatics/rdptools/2.0.2/"
WORKFOLDER=${HOMEFOLDER}/nifH_work
RESULTSFOLDER=${HOMEFOLDER}/nifH_results
RESOURCEFOLDER=${HOMEFOLDER}/Resources
INPUTREFERENCEALIGNMENT=cluster.rep.nifH-chL-bchX.0.9_noDot.fasta
INPUTREFERENCETREE=RAxML_bipartitions.MultipleOriginal.tree
id=`echo 0.94 | bc` # OTU clustering, 1.0 = zOTUS, 0.97 = 97%
cores=4 # number of cores to use
MINLEN=200 # merged reads with a smaller size will be removed

if [ -d ${WORKFOLDER} ]
then
 rm -r ${WORKFOLDER}
fi
if [ -d ${RESULTSFOLDER} ]
then
 rm -r ${RESULTSFOLDER}
fi
mkdir ${WORKFOLDER}
mkdir ${RESULTSFOLDER}

echo "NifMAP1.2"> ${RESULTSFOLDER}/NifMAP.log

# 1. Merge MiSeq aplicon reads
echo "usearch -fastq_mergepairs">> ${RESULTSFOLDER}/NifMAP.log
cd $inputDir
#for gzfile in *.fastq.gz; do; gunzip $gzfile; done
usearch -fastq_mergepairs *_R1.fastq -fastqout merged_reads.fq  -relabel @  -threads $cores --log ${RESULTSFOLDER}/tmp.log
cd ../
cat ${RESULTSFOLDER}/tmp.log >>  ${RESULTSFOLDER}/NifMAP.log
echo "--">> ${RESULTSFOLDER}/NifMAP.log

# 2. Quality filter merged reads
echo "usearch -fastq_filter">> ${RESULTSFOLDER}/NifMAP.log
#for inputReads in ${inputDir}/*.fastq; do
#    usearch -fastq_filter ${inputReads} -fastq_maxee 1.0 -fastq_minlen 200 -relabel @ -fastaout ${WORKFOLDER}/${inputReads%.*}.fasta --log ${RESULTSFOLDER}/tmp.log;
#    cat ${RESULTSFOLDER}/tmp.log >>  ${RESULTSFOLDER}/NifMAP.log
#done
usearch -fastq_filter ${inputDir}/merged_reads.fq -fastq_maxee 1.0 -fastq_minlen $MINLEN -fastaout ${WORKFOLDER}/merged_reads.fa -threads $cores --log ${RESULTSFOLDER}/tmp.log
cat ${RESULTSFOLDER}/tmp.log >>  ${RESULTSFOLDER}/NifMAP.log
echo "--" >> ${RESULTSFOLDER}/NifMAP.log

# 3. Filter sequences using nifH HMM
# Screen merged reads using HMMer
echo "Screen merged reads using HMMer">> ${RESULTSFOLDER}/NifMAP.log
hmmsearch --cpu $cores --domtblout ${WORKFOLDER}/hmmOut1.out -o ${WORKFOLDER}/junk ${RESOURCEFOLDER}/hmm_nuc_1160_nifH.hmm ${WORKFOLDER}/merged_reads.fa
awk '{print $1}' ${WORKFOLDER}/hmmOut1.out | grep -v "#" > ${WORKFOLDER}/acceptable_hits
grep ">" ${WORKFOLDER}/merged_reads.fa | grep -v -F -f ${WORKFOLDER}/acceptable_hits >${WORKFOLDER}/shitHits
totalUnique=`grep ">" ${WORKFOLDER}/merged_reads.fa | wc -l`
totalAccepted=`cat ${WORKFOLDER}/acceptable_hits | wc -l`
totalRemoved=`cat ${WORKFOLDER}/shitHits | wc -l`
echo "hmmscreen of nifH removed ${totalRemoved} sequences out of ${totalUnique} sequences. ${totalAccepted} unique sequences retained" >> ${RESULTSFOLDER}/NifMAP.log
mv ${WORKFOLDER}/merged_reads.fa ${WORKFOLDER}/merged_reads_prehmm.fa
awk 'BEGIN{FS="\n";RS=">"};NR>1{print(">"$1);for(i=2;i<=NF;i++){printf($i)};print("")}' ${WORKFOLDER}//merged_reads_prehmm.fa | grep -A 1 -F -f ${WORKFOLDER}/acceptable_hits | grep -v "^\-\-$" >${RESULTSFOLDER}/merged_reads.fa
echo "--">> ${RESULTSFOLDER}/NifMAP.log

# 4. Generate OTUs
# Dereplicate sequences for OTU table
echo "usearch -fastx_uniques">> ${RESULTSFOLDER}/NifMAP.log
echo "Singletons will not be allowed to generate OTUs">> ${RESULTSFOLDER}/NifMAP.log
usearch -fastx_uniques ${RESULTSFOLDER}/merged_reads.fa -minuniquesize 2 -fastaout ${RESULTSFOLDER}/unique.fa -threads $cores --log ${RESULTSFOLDER}/tmp.log # -sizeout
cat ${RESULTSFOLDER}/tmp.log >>  ${RESULTSFOLDER}/NifMAP.log
echo "--">> ${RESULTSFOLDER}/NifMAP.log

# Make zOTUs or OTUs
if [[ (( $id == 1.0 )) ]]; then
    echo "usearch -unoise3 ">> ${RESULTSFOLDER}/NifMAP.log
    usearch -sortbylength ${RESULTSFOLDER}/unique.fa -fastaout ${WORKFOLDER}/unique_sorted.fa
    usearch -unoise3 ${WORKFOLDER}/unique_sorted.fa -zotus ${WORKFOLDER}/otus.fa -tabbedout ${RESULTSFOLDER}/unoise3.txt -threads $cores --log ${RESULTSFOLDER}/tmp.log
    cat ${RESULTSFOLDER}/tmp.log >>  ${RESULTSFOLDER}/NifMAP.log
    echo "--">> ${RESULTSFOLDER}/NifMAP.log
else
    echo "usearch -cluster_fast">> ${RESULTSFOLDER}/NifMAP.log
    usearch -sortbylength ${RESULTSFOLDER}/unique.fa -fastaout ${WORKFOLDER}/unique_sorted.fa
    usearch -cluster_fast ${WORKFOLDER}/unique_sorted.fa -id $id -centroids ${WORKFOLDER}/otus.fa -uc ${WORKFOLDER}/clusters.uc -threads $cores --log ${RESULTSFOLDER}/tmp.log
    cat ${RESULTSFOLDER}/tmp.log >>  ${RESULTSFOLDER}/NifMAP.log
    echo "--">> ${RESULTSFOLDER}/NifMAP.log
fi


# 5.Translate and correct OTU representatives using Framebot
#eVal_chL=1e-50
#score_chL=150
#eVal_bchX=1e-50
#score_bchX=150

java -jar ${FRAMEBOTPATH}/FrameBot.jar framebot -N -l 30 -i 0.4 -o ${WORKFOLDER}/nifH ${FRAMEBOTPATH}/Framebot/refset/nifh_prot_ref.fasta ${WORKFOLDER}/otus.fa
cp ${WORKFOLDER}/nifH_corr_nucl.fasta ${RESULTSFOLDER}


# 6. Make an OTU table based on corrected seqs
#echo "usearch -usearch_global ">> ${RESULTSFOLDER}/NifMAP.log
#usearch -usearch_global ${RESULTSFOLDER}/merged_reads.fa -db ${RESULTSFOLDER}/nifH_corr_nucl.fasta  -strand plus -id 0.97 -threads $cores -otutabout ${RESULTSFOLDER}/otu_table.txt --log ${RESULTSFOLDER}/tmp.log
#cat ${RESULTSFOLDER}/tmp.log >>  ${RESULTSFOLDER}/NifMAP.log
#echo "--">> ${RESULTSFOLDER}/NifMAP.log

# 7. Filter out homologous genes (bchL; chlL; bchX; parA) using HMM:
# Corrected AA sequences are in nifH_corr_prot.fasta
# Screen with hmm to identify all hits
hmmscan --cpu $cores  --domtblout ${WORKFOLDER}/hmmOut2.out -o ${WORKFOLDER}/full_output_NifH_ChlL_bchX ${RESOURCEFOLDER}/NifH_ChlL_bchX.hmm ${WORKFOLDER}/nifH_corr_prot.fasta
#Rscript --vanilla < ${HOMEFOLDER}/nifH_bch_hmmEvaluation.R ${WORKFOLDER}/hmmOut2.out
#mv nifH_bch_hmmEvaluation.pdf ${RESULTSFOLDER}
cat ${WORKFOLDER}/hmmOut2.out | awk 'NR>3{if($8>bitarray[$4]){bitarray[$4]=$8;outArray[$4]=$1"\t"$4}}END{for(entry in outArray){print outArray[entry]}}' > ${WORKFOLDER}/assignments.txt
cp ${WORKFOLDER}/hmmOut2.out ${RESULTSFOLDER}/nifH_bch_hmmEvaluation.hmm.out
grep "nifH" ${WORKFOLDER}/assignments.txt | awk '{print $2}' | sort > ${WORKFOLDER}/acceptable_hits
grep ">" ${RESULTSFOLDER}/nifH_corr_nucl.fasta | awk '{print $1}' | grep -v -F -f ${WORKFOLDER}/acceptable_hits | sed 's/>//'>${WORKFOLDER}/shitHits
totalOTUs=`grep ">" ${RESULTSFOLDER}/nifH_corr_nucl.fasta | wc -l`
totalAccepted=`cat ${WORKFOLDER}/acceptable_hits | wc -l`
totalRemoved=`cat ${WORKFOLDER}/shitHits | wc -l`
echo "FRAMEBOT and hmmscreen of nifH removed ${totalRemoved} sequences out of ${totalOTUs} ASVs. ${totalAccepted} ASV reps retained">>${RESULTSFOLDER}/NifMAP.log

#mv ${RESULTSFOLDER}/${inputReads} ${WORKFOLDER}/${inputReads}
cat ${RESULTSFOLDER}/nifH_corr_nucl.fasta | sed 's/ //g' | awk 'BEGIN{RS=">";FS="\n"};NR>1{printf $1"_\t";for(i=2;i<=NF;i++){printf($i)}print("")}' | grep -F -f ${WORKFOLDER}/acceptable_hits | awk '{gsub("_$","",$1);print(">"$1"\n"$2)}' > ${RESULTSFOLDER}/nifH_corr_nucl_only_nifH.fasta
cat ${WORKFOLDER}/nifH_corr_prot.fasta | sed 's/ //g' |  awk 'BEGIN{RS=">";FS="\n"};NR>1{printf $1"_\t";for(i=2;i<=NF;i++){printf($i)}print("")}' | grep -F -f ${WORKFOLDER}/acceptable_hits | awk '{gsub("_$","",$1);print(">"$1"\n"$2)}' > ${RESULTSFOLDER}/nifH_corr_prot.fasta
cat ${WORKFOLDER}/nifH_corr_prot.fasta | sed 's/ //g' |  awk 'BEGIN{RS=">";FS="\n"};NR>1{printf $1"_\t";for(i=2;i<=NF;i++){printf($i)}print("")}' | grep -F -f ${WORKFOLDER}/shitHits | awk '{gsub("_$","",$1);print(">"$1"\n"$2)}' > ${RESULTSFOLDER}/nifH_rej_prot.fasta

# Make an OTU table based on filter corrected seqs
echo "usearch -usearch_global ">> ${RESULTSFOLDER}/NifMAP.log
#usearch -usearch_global ${RESULTSFOLDER}/merged_reads.fa -db ${RESULTSFOLDER}/nifH_corr_nucl_only_nifH.fasta -id 1.0 -strand plus -threads $cores -otutabout ${RESULTSFOLDER}/otu_table.txt
usearch -otutab ${RESULTSFOLDER}/merged_reads.fa -otus ${RESULTSFOLDER}/nifH_corr_nucl_only_nifH.fasta -otutabout ${RESULTSFOLDER}/otu_table.txt -mapout ${WORKFOLDER}/map.txt -threads $cores

# 7. Place OTUs on a reference tree using RAxML-EPA
#RAxML-EPA
mkdir ${RESULTSFOLDER}/EPA
cp ${RESULTSFOLDER}/nifH_corr_prot.fasta ${RESULTSFOLDER}/EPA/
cp ${RESOURCEFOLDER}/${INPUTREFERENCEALIGNMENT} ${RESULTSFOLDER}/EPA/
cp ${RESOURCEFOLDER}/${INPUTREFERENCETREE} ${RESULTSFOLDER}/EPA/

cd ${RESULTSFOLDER}/EPA
mafft --add nifH_corr_prot.fasta --thread 16 ${INPUTREFERENCEALIGNMENT} > completeAlignment.fa
cat completeAlignment.fa | awk '{print $1}' | sed 's/|$//'> RAxML_compatible.fa
raxmlHPC-PTHREADS-SSE3 -f v -T $cores -s RAxML_compatible.fa -m PROTCATJTT -t ${INPUTREFERENCETREE} -n EPAplaced -p 123

cd ../../

# 8. Classify sequences using CART
cat ${RESOURCEFOLDER}/AztVine8_AA.fasta ${RESULTSFOLDER}/nifH_corr_prot.fasta > ${WORKFOLDER}/nifH_corr_prot_4classification.fa
hmmalign --amino ${RESOURCEFOLDER}/Zehr_2014_1812genomes_nifH_AA_noGaps.hmm ${WORKFOLDER}/nifH_corr_prot_4classification.fa > ${WORKFOLDER}/nifH_corr_prot_hmm.sth
conda activate base
seqmagick convert ${WORKFOLDER}/nifH_corr_prot_hmm.sth  ${WORKFOLDER}/nifH_corr_prot_hmmAln.fasta
conda deactivate
python2 ${RESOURCEFOLDER}/NifH_Clusters.py ${WORKFOLDER}/nifH_corr_prot_hmmAln 1
grep '>' ${WORKFOLDER}/nifH_corr_prot_hmmAln_Clusters.fasta |  sed 's/>//' | awk '{print $1, $5, $8}' > ${RESULTSFOLDER}/nifH_corr_prot.zehr.classification

# 9. Classify sequences using BLASTP
wget http://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz ./
tar -zxvf taxdb.tar.gz

blastp -query ${RESULTSFOLDER}/nifH_corr_prot.fasta -db /proj/Resources/BLAST/NCBI/refseq_protein -num_threads 80 -max_target_seqs 1 -out ${RESULTSFOLDER}/nifH_blastp.out -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sblastnames sskingdoms stitle"

# 8. Return results and clean up
#return results
#cp ${RESULTSFOLDER} $HOMEFOLDER/

#cleanup
rm -r ${RESULTSFOLDER}/tmp.log ${RESULTSFOLDER}/unique.fa ${RESULTSFOLDER}/unoise3.txt
rm -r ${WORKFOLDER}
rm taxdb*
