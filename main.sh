#!/bin/sh
#$ -S /bin/sh
#$ -cwd

###### write path of the required tools #######
python=python
java=java
NanoFilt=NanoFilt
canu=canu
flye=flye
minimap2=minimap2
samtools=samtools
blastn=blastn
jvarkit_dist=jvarkit/dist
php=php
multi_csar=multi-csar.php
###############################################

###### simple usage ###########################
#
# bash scALA.sh -i PooledSAGLongRead.fastq -o OutDir -g 3.1m -d 50
#
###############################################


FLG_fastq="FALSE"
FLG_outdir="FALSE"
FLG_maxdepth="FALSE"
FLG_genomesize="FALSE"
FLG_reference="FALSE"
FLG_SRcontig="FALSE"
FLG_SubDebias="FALSE"

while getopts i:o:d:g:r:s:a OPT
do
 case $OPT in
  "i" ) FLG_seqdir="TRUE" ; ifastq=$OPTARG ;;         # required. Input LongRead Fastq File (Not Compressed)
  "o" ) FLG_outdir="TRUE" ; outdir=$OPTARG ;;         # required. Output Directory
  "d" ) FLG_maxdepth="TRUE" ; maxdepth=$OPTARG ;;     # required. Threshold of Depth in Debiasing Step
  "g" ) FLG_genomesize="TRUE" ; genomesize=$OPTARG ;; # required. Estimated Genome Size of Target Bacteria
  "r" ) FLG_reference="TRUE" ; reference=$OPTARG ;;   # optional. If you have a known reference genome.
  "s" ) FLG_SRcontig="TRUE" ; SRcontig=$OPTARG ;;     # optional. If you want to use short-read contig in the scaffolding step. 
 esac
done

if [ "$FLG_seqdir" != "TRUE" ] || [ "$FLG_outdir" != "TRUE" ] || [ "$FLG_maxdepth" != "TRUE" ] || [ "$FLG_genomesize" != "TRUE" ]; then
 echo -e "no required option"
 exit
fi


mkdir $outdir -p

## QC of LR fastq
mkdir $outdir/Fastq -p
samplename=$(basename $ifastq .fastq)
fastq=$outdir/Fastq/${samplename}_filt.fastq
$NanoFilt -q 10 -l 1000 --headcrop 75 $ifastq > $fastq

## First assembly to get reference contig for debiasing
mkdir $outdir/Assembly -p
if [ $FLG_reference == "FALSE" ]; then
 refcanuout=${outdir}/Assembly/Refcanu
 $canu -d $refcanuout -p ${samplename}_ref genomeSize=$genomesize \
       -nanopore-raw $fastq useGrid=false readSamplingCoverage=$maxdepth \
       saveReads=true stopAfter=trimming
 gunzip ${refcanuout}/${samplename}_ref.trimmedReads.fasta.gz
 awk '{if( (NR-1)%2 ) print; else printf(">Indx-%d\n",cnt++)}' ${refcanuout}/${samplename}_ref.trimmedReads.fasta \
     > ${refcanuout}/${samplename}_ref.trimmedReads.rename.fasta
 gzip ${refcanuout}/${samplename}_ref.trimmedReads.fasta
 $flye --nano-raw ${refcanuout}/${samplename}_ref.trimmedReads.rename.fasta \
       --out-dir ${refcanuout}/flye --genome-size $genomesize --plasmids
 cp ${refcanuout}/flye/assembly.fasta ${outdir}/${samplename}_ref.fasta
 reference=${outdir}/${samplename}_ref.fasta
fi

## Debiasing-Assembly Loop
count=1

mkdir $outdir/Mapping_debias -p

while [ $count -lt 5 ]; do
 # Debiasing
 refname=$(basename $reference .fasta) 
 $minimap2 -ax map-ont $reference $fastq | $samtools view -b | $samtools sort > $outdir/Mapping_debias/minimap2_${refname}_sorted.bam
 $java -jar -Xmx16g -Xms10g $jvarkit_dist/sortsamrefname.jar --samoutputformat BAM $outdir/Mapping_debias/minimap2_${refname}_sorted.bam \
                                                             -o $outdir/Mapping_debias/minimap2_${refname}_namesorted.bam
 $java -jar -Xmx16g -Xms10g $jvarkit_dist/biostar154220.jar -n $maxdepth --samoutputformat BAM \
                                                            $outdir/Mapping_debias/minimap2_${refname}_namesorted.bam \ 
                            | $samtools sort > $outdir/Mapping_debias/minimap2_${refname}_namesorted_maxdepth${maxdepth}.bam
 $samtools fastq $outdir/Mapping_debias/minimap2_${refname}_namesorted_maxdepth${maxdepth}.bam > $outdir/Fastq/${samplename}_debiased_${count}.fastq
 rm $outdir/Mapping_debias/minimap2_${refname}_*bam
 
 # re-assembly
 assemblyout=${outdir}/Assembly/Debiased_$count
 $canu -d ${assemblyout} -p debiased_$count genomeSize=$genomesize -nanopore-raw \
       $outdir/Fastq/${samplename}_debiased_${count}.fastq useGrid=false readSamplingCoverage=$maxdepth saveReads=true stopAfter=trimming
 gunzip ${assemblyout}/debiased_${count}.trimmedReads.fasta.gz
 awk '{if( (NR-1)%2 ) print; else printf(">Indx-%d\n",cnt++)}' ${assemblyout}/debiased_${count}.trimmedReads.fasta \
      > ${assemblyout}/debiased_${count}.trimmedReads.rename.fasta
 gzip ${assemblyout}/debiased_${count}.trimmedReads.fasta
 $flye --nano-corr ${assemblyout}/debiased_${count}.trimmedReads.rename.fasta \
       --out-dir ${assemblyout}/flye --genome-size $genomesize --plasmids
 cp ${assemblyout}/flye/assembly.fasta ${outdir}/${samplename}_debiased_${count}.fasta

 #Remove duplicated short contig
 Seq=${outdir}/${samplename}_debiased_${count}.fasta
 flag=TRUE
 while [ $flag = "TRUE" ]; do
  $blastn -query $Seq -subject $Seq -outfmt "6 qseqid qlen sseqid slen pident length qstart qend sstart send" \
          -evalue 1e-100 | sort -k1,1 -k7,7nr > ${outdir}/tmpblast.txt
  $python script/BLAST_dedup.py -i $Seq -b ${outdir}/tmpblast.txt -o ${outdir}/tmp.fasta
  Before=$(ls -l $Seq | cut -f5 -d " ")
  After=$(ls -l ${outdir}/tmp.fasta | cut -f5 -d " ")
  if [ "$Before" != "$After" ]; then
   mv ${outdir}/tmp.fasta $Seq
  else
   flag=FALSE
  fi
 done
 rm ${outdir}/tmp.fasta ${outdir}/tmpblast.txt

 reference=${outdir}/${samplename}_debiased_${count}.fasta
 count=$(( count + 1 ))
done

### BLAST scaffolding
ls ${outdir}/*debias*fasta > ${outdir}/ifiles
if [ $FLG_SRcontig = "TRUE" ]; then
 ls $SRcontig >> ${outdir}/ifiles
fi
target=$(seqkit stats  ${outdir}/*debias*fasta | sort -k4,4n -k8,8nr | head -n 2 | tail -n1 | cut -f1 -d " ")
cp $target ${outdir}/${samplename}_BLASTscaffold.fasta
subject=${outdir}/${samplename}_BLASTscaffold.fasta
cat ${outdir}/ifiles | while read line
do
 flag=TRUE
 while [ $flag = "TRUE" ]; do
  $blastn -query $line -subject $subject -outfmt "6 qseqid qlen sseqid slen pident length qstart qend sstart send" \
          -evalue 1e-100 -culling_limit 1 | sort -k1,1 -k7,7n > ${outdir}/tmpblast.txt
  $python script/BLAST_scaffold.py -b ${outdir}/tmpblast.txt -q $line -s $subject \
          -o ${outdir}/${samplename}_tmp_scaffold.fasta
  Before=$(wc -l $subject | cut -f1 -d " ")
  After=$(wc -l ${outdir}/${samplename}_tmp_scaffold.fasta | cut -f1 -d " ")
  if [ "$Before" != "$After" ]; then
   mv ${outdir}/${samplename}_tmp_scaffold.fasta $subject
  else
   flag=FALSE
  fi
 done
done
rm ${outdir}/${samplename}_tmp_scaffold.fasta ${outdir}/tmpblast.txt

### Multicsar scaffolding
mkdir -p ${outdir}/TmpConnect
cp ${outdir}/*debias*fasta ${outdir}/TmpConnect
if [ "$FLG_SRcontig" = "TRUE" ]; then
 cp $SRcontig ${outdir}/TmpConnect
fi
csarout=${outdir}/MultiCSAR
$php $multi_csar -t $target -r ${outdir}/TmpConnect --nuc -o ${csarout}
cat $target  | sed -z "s/\n/\t/g" | sed -z "s/\t>/\n>/g" | sed -r "s/\t/\n/" | sed -r "s/\t//g" > ${outdir}/TmpConnect/tmp.fasta
grep -v ">" ${csarout}/multi-csar.nuc.out | grep "_" | cut -f1 -d " " > ${csarout}/connected.txt
grep ">" ${outdir}/TmpConnect/tmp.fasta | cut -f2 -d ">" > ${csarout}/total.txt
sort ${csarout}/total.txt ${csarout}/connected.txt | uniq -u > ${csarout}/remain.txt
grep -f ${csarout}/remain.txt -A1 ${outdir}/TmpConnect/tmp.fasta | grep -v "^-" > ${csarout}/remain.fasta
cat ${csarout}/multi-csar.nuc.out.fna ${csarout}/remain.fasta > ${outdir}/${samplename}_CSARscaffold.fasta
rm ${outdir}/TmpConnect -r
