#!/bin/bash -
#title          :NifMAP.sh
#description    :NifMAP - a bioinformatics pipline for analysing nifH amplicon data
#author         :Roey Angel
#date           :20180122
#version        :1.0
#usage          :./nifmap.sh $inputReads.fa
#notes          :This pipeline is designed to work in concert with the pipeline described by Herbold et al. 2015. See Angel et al. 2017 for details.
#dependencies   : 1. HMMER - http://hmmer.org/
#                 2. FrameBot - https://github.com/rdpstaff/Framebot
#                 3. seqmagick - https://fhcrc.github.io/seqmagick/
#                 4. CART - https://wwwzehr.pmc.ucsc.edu/CART_model_public/
#                 5. R - https://www.r-project.org/
#bash_version   :4.3.48(1)-release
#============================================================================

# 1. Merge MiSeq aplicon reads

# 2. Quality filter merged reads

# 3. Filter sequences using nifH HMM
# Define variables
inputReads=$1 # input fasta file containing merged MiSeq reads
HOMEFOLDER=`pwd`
FRAMEBOTPATH="/apps/rdptools/20160204"
WORKFOLDER=/tmp/$USER/NifMAP_work
RESULTSFOLDER=/tmp/$USER/NifMAP_results
RESOURCEFOLDER=$HOMEFOLDER/Resources

if [ ! -d /tmp/$USER ]
then
 mkdir /tmp/$USER
fi
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
cp ${HOMEFOLDER}/${inputReads} ${WORKFOLDER}

# Screen merged reads using HMMer
hmmsearch --domtblout ${WORKFOLDER}/hmmOut.out -o ${WORKFOLDER}/junk ${RESOURCEFOLDER}/hmm_nuc_1160_nifH.hmm ${WORKFOLDER}/${inputReads}
awk '{print $1}' ${WORKFOLDER}/hmmOut.out | grep -v "#" > ${WORKFOLDER}/acceptable_hits
grep ">" ${WORKFOLDER}/${inputReads} | grep -v -F -f ${WORKFOLDER}/acceptable_hits >${WORKFOLDER}/shitHits
totalUnique=`grep ">" ${WORKFOLDER}/${inputReads} | wc -l`
totalAccepted=`cat ${WORKFOLDER}/acceptable_hits | wc -l`
totalRemoved=`cat ${WORKFOLDER}/shitHits | wc -l`
echo "hmmscreen of nifH removed ${totalRemoved} sequences out of ${totalUnique} sequences. ${totalAccepted} unique sequences retained"> ${RESULTSFOLDER}/NifMAP_log.txt
mv ${WORKFOLDER}/${inputReads} ${WORKFOLDER}/${inputReads%.fa}_prehmm.fa
awk 'BEGIN{FS="\n";RS=">"};NR>1{print(">"$1);for(i=2;i<=NF;i++){printf($i)};print("")}' ${WORKFOLDER}/${inputReads%.fa}_prehmm.fa | grep -A 1 -F -f ${WORKFOLDER}/acceptable_hits | grep -v "^\-\-$" >${RESULTSFOLDER}/${inputReads}

# 4. Generate OTUs

# 5.Translate and correct OTU representatives using Framebot
eVal_chL=1e-50
score_chL=150
eVal_bchX=1e-50
score_bchX=150

cd ${WORKFOLDER}
java -jar ${FRAMEBOTPATH}/FrameBot.jar framebot -N -l 30 -i 0.4 -o ${inputReads%.fa} ${FRAMEBOTPATH}/Framebot/refset/nifh_prot_ref.fasta ${RESULTSFOLDER}/${inputReads}

# 6. Filter out homologous genes (bchL; chlL; bchX; parA) using HMM:
# Corrected AA sequences are in ${inputReads%.fa}_corr_prot.fasta
# Screen with hmm to identify all hits
hmmscan --domtblout ${WORKFOLDER}/hmmOut.out -o ${WORKFOLDER}/junk_NifH_ChlL_bchX ${RESOURCEFOLDER}/nifH_chlL_bchX.hmm ${WORKFOLDER}/${inputReads%.fa}_corr_prot.fasta
Rscript --vanilla < ${HOMEFOLDER}/nifH_bch_hmmEvaluation.R ${WORKFOLDER}/hmmOut.out
mv nifH_bch_hmmEvaluation.pdf ${RESULTSFOLDER}
cat hmmOut.out | awk 'NR>3{if($8>bitarray[$4]){bitarray[$4]=$8;outArray[$4]=$1"\t"$4}}END{for(entry in outArray){print outArray[entry]}}' > assignments.txt
cp hmmOut.out ${RESULTSFOLDER}/nifH_bch_hmmEvaluation.hmm.out
grep "nifH" assignments.txt | awk '{print $2"_"}' | sort > ${WORKFOLDER}/acceptable_hits
grep ">" ${RESULTSFOLDER}/${inputReads} | awk '{print $1"_"}' | grep -v -F -f ${WORKFOLDER}/acceptable_hits | sed 's/>//'>${WORKFOLDER}/shitHits
totalOTUs=`grep ">" ${RESULTSFOLDER}/${inputReads} | wc -l`
totalAccepted=`cat ${WORKFOLDER}/acceptable_hits | wc -l`
totalRemoved=`cat ${WORKFOLDER}/shitHits | wc -l`
echo "FRAMEBOT and hmmscreen of ${GENE} removed ${totalRemoved} sequences out of ${totalOTUs} OTU sequences. ${totalAccepted} otu reps retained">>${RESULTSFOLDER}/NifMAP_log.txt

mv ${RESULTSFOLDER}/${inputReads} ${WORKFOLDER}/${inputReads}
cat ${WORKFOLDER}/${inputReads} | sed 's/ //g' | awk 'BEGIN{RS=">";FS="\n"};NR>1{printf $1"_\t";for(i=2;i<=NF;i++){printf($i)}print("")}' | grep -F -f ${WORKFOLDER}/acceptable_hits | awk '{gsub("_$","",$1);print(">"$1"\n"$2)}' > ${RESULTSFOLDER}/${inputReads}
cat ${WORKFOLDER}/${inputReads%.fa}_corr_prot.fasta | sed 's/ //g' |  awk 'BEGIN{RS=">";FS="\n"};NR>1{printf $1"_\t";for(i=2;i<=NF;i++){printf($i)}print("")}' | grep -F -f ${WORKFOLDER}/acceptable_hits | awk '{gsub("_$","",$1);print(">"$1"\n"$2)}' > ${RESULTSFOLDER}/${inputReads%.fa}_corr_prot.fasta
cat ${WORKFOLDER}/${inputReads%.fa}_corr_prot.fasta | sed 's/ //g' |  awk 'BEGIN{RS=">";FS="\n"};NR>1{printf $1"_\t";for(i=2;i<=NF;i++){printf($i)}print("")}' | grep -F -f ${WORKFOLDER}/shitHits | awk '{gsub("_$","",$1);print(">"$1"\n"$2)}' > ${RESULTSFOLDER}/${inputReads%.fa}_rej_prot.fasta

# 7. Classify sequences using BLASTP, CART and RAxML-EPA
#RAxML-EPA
INPUTREFERENCEALIGNMENT=cluster.rep.nifH-chL-bchX.0.9_noDot.fasta.gz
INPUTREFERENCETREE=RAxML_bipartitions.MultipleOriginal.tree.gz
cp ${RESULTSFOLDER}/${inputReads%.fa}_corr_prot.fasta ./
cp ${RESOURCEFOLDER}/${INPUTREFERENCEALIGNMENT} ./
cp ${RESOURCEFOLDER}/${INPUTREFERENCETREE} ./

gunzip ${INPUTREFERENCEALIGNMENT}
gunzip ${INPUTREFERENCETREE}

mafft --add ${inputReads%.fa}_corr_prot.fasta --thread 16 ${INPUTREFERENCEALIGNMENT%.gz} > completeAlignment.fa

cat completeAlignment.fa | awk '{print $1}' | sed 's/|$//'> RAxml.compatible.fa
raxmlHPC-PTHREADS-SSE3 -f v -T 16 -s RAxml.compatible.fa -m PROTCATJTT -t RAxML_bipartitions.MultipleOriginal.tree -n EPAplaced -p 123
mkdir EPA
cp *EPAplaced* EPA
cp completeAlignment.fa EPA
mv EPA ${RESULTSFOLDER}

# CART
cat ${RESOURCEFOLDER}/CART_nifH/AztVine8_AA.fasta ${RESULTSFOLDER}/${inputReads%.fa}_corr_prot.fasta > ${WORKFOLDER}/${inputReads%.fa}_corr_prot_4classification.fa
hmmalign --amino ${RESOURCEFOLDER}/Zehr_2014_1812genomes_nifH_AA_noGaps.hmm ${WORKFOLDER}/${inputReads%.fa}_corr_prot_4classification.fa > ${WORKFOLDER}/${inputReads%.fa}_corr_prot_hmm.sth
seqmagick convert  ${WORKFOLDER}/${inputReads%.fa}_corr_prot_hmm.sth  ${WORKFOLDER}/${inputReads%.fa}_corr_prot_hmmAln.fasta
python ${RESOURCEFOLDER}/CART_nifH/NifH_Clusters.py ${WORKFOLDER}/${inputReads%.fa}_corr_prot_hmmAln 1
grep '>' ${WORKFOLDER}/${inputReads%.fa}_corr_prot_hmmAln_Clusters.fasta |  sed 's/>//' | awk '{print $1, $5, $8}' > ${RESULTSFOLDER}/${inputReads%.fa}_corr_prot.zehr.classification

# BLASTP
wget http://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz ./
tar -zxvf taxdb.tar.gz
blastp -query ${RESULTSFOLDER}/${inputReads%.fa}_corr_prot.fasta -db /localmirror/monthly/blast/ncbi/refseq_protein -num_threads 16 -max_target_seqs 1 -out ${RESULTSFOLDER}/${inputReads%.fa}_refseq-prot_nifH.out -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sblastnames sskingdoms stitle"

# 8. Return results and clean up
#return results
cp ${RESULTSFOLDER} $HOMEFOLDER/

#cleanup
rm -r ${WORKFOLDER}
rm -r ${RESULTSFOLDER}
