#!/bin/sh

PSHOME=/zfssz5/BC_PS/chengyuanfang/work/CGIA_3.0_module/HLA_prediction/config
NOVOALIGN_DIR=/hwfssz4/BC_COM_FP/bc_tumor/Pipeline/Cancer_Genomics_Interpretation_Analysis/CGIA_v2.0/bin/novocraft/novoalign
NUM_THREADS=1

export PERL5LIB=/hwfssz4/BC_COM_FP/bc_tumor/Pipeline/Cancer_Genomics_Interpretation_Analysis/CGIA_v2.0/module/NeoantigenPrediction/NeoantigenPrediction_netMHCpan/lib/perl_lib:$PERL5LIB
export PSHOME
export NOVOALIGN_DIR
export NUM_THREADS

#### check if an appropriate number of arguments were passed ####

if [ $# -ne 5 ]; then
        echo 1>&2 Usage: $0 bam/cram race build format outDir
	echo "	-bam: path to the BAM/CRAM file to be used for HLA typing"
	echo "	-race: ethnicity of the individual (Caucasian, Black, Asian or Unknown)"
	echo "	-build: hg38_ucsc or hg38_gencode or hg19"
	echo "	-format: fastq format (STDFQ, ILMFQ, ILM1.8 or SLXFQ; see Novoalign documentation)"
	echo "	-outDir: output directory"
        exit 127
fi

echo "Environment variables"
echo "  -PSHOME: POLYSOLVER home directory = $PSHOME"
                                                                   
bam=$1
race=$2
build=$3
format=$4
outDir=$5
iFile=0 #dont cal insertsize

ids=$PSHOME/files/ids
tag_file=$PSHOME/files/abc_v14.uniq
bam2fastq=$PSHOME/bin/bam2fastq
samtools=$PSHOME/bin/samtools
gatk=$PSHOME/bin/gatk
mkdir $outDir
mkdir $outDir/tmp
TMP_DIR=$outDir/tmp
var=`cat $ids`

# getting matching tag sequences

echo -n -e "getting matching tags\n"

if [[ $bam == *bam ]] && [ ! -f $bam".bai" ]
then
	$samtools index $bam
fi

if [[ $bam == *cram ]] && [ ! -f $bam".crai" ]
then
        $samtools index $bam
fi


$samtools view -H $bam > $outDir/tag.sam
$samtools view $bam | grep -F -f $tag_file >> $outDir/tag.sam
$samtools view -bS -o $outDir/tag.bam $outDir/tag.sam

$bam2fastq -o $outDir/tag#.fastq --aligned -q $outDir/tag.bam 
mv $outDir/tag_1.fastq $outDir/tag.1.fastq  
mv $outDir/tag_2.fastq $outDir/tag.2.fastq 

$PSHOME/bin/clean_unpaired_fastq.pl $outDir/tag.1.fastq
$PSHOME/bin/clean_unpaired_fastq.pl $outDir/tag.2.fastq
gzip $outDir/tag.1.fastq
gzip $outDir/tag.2.fastq

#getting chr6 region

echo -n -e "getting chr6 region\n"

$samtools view -H $bam > $outDir/chr6region.sam

if [ $build == "hg19" ]; then
        echo -n -e "build=hg19\n"
        $samtools view $bam 6:29909037-29913661 >> $outDir/chr6region.sam
	$samtools view $bam chr6:29909037-29913661 >> $outDir/chr6region.sam
        $samtools view $bam 6:31321649-31324964 >> $outDir/chr6region.sam
	$samtools view $bam chr6:31321649-31324964 >> $outDir/chr6region.sam
        $samtools view $bam 6:31236526-31239869 >> $outDir/chr6region.sam
	$samtools view $bam chr6:31236526-31239869 >> $outDir/chr6region.sam
###HLA-A_chr6_hap####
	$samtools view $bam chr6_apd_hap1:1148926-1153548 >> $outDir/chr6region.sam
        $samtools view $bam chr6_cox_hap2:1420791-1425413 >> $outDir/chr6region.sam
        $samtools view $bam chr6_dbb_hap3:1202397-1207037 >> $outDir/chr6region.sam
        $samtools view $bam chr6_mann_hap4:1201421-1206062 >> $outDir/chr6region.sam
        $samtools view $bam chr6_mcf_hap5:1150160-1295548 >> $outDir/chr6region.sam
        $samtools view $bam chr6_qbl_hap6:1201814-1206454 >> $outDir/chr6region.sam
        $samtools view $bam chr6_ssto_hap7:1238371-1243012 >> $outDir/chr6region.sam
###HLA-B_chr6_hap####
        $samtools view $bam chr6_cox_hap2:2834338-2837710 >> $outDir/chr6region.sam
        $samtools view $bam chr6_mann_hap4:2668104-2671477 >> $outDir/chr6region.sam
        $samtools view $bam chr6_mcf_hap5:2701429-2704815 >> $outDir/chr6region.sam
        $samtools view $bam chr6_qbl_hap6:2615165-2618494 >> $outDir/chr6region.sam
        $samtools view $bam chr6_ssto_hap7:2655408-2658765 >> $outDir/chr6region.sam
###HLA-C_chr6_hap####
        $samtools view $bam chr6_cox_hap2:2749781-2753162 >> $outDir/chr6region.sam
        $samtools view $bam chr6_dbb_hap3:2532134-2535505 >> $outDir/chr6region.sam
        $samtools view $bam chr6_mann_hap4:2583421-2586792 >> $outDir/chr6region.sam
        $samtools view $bam chr6_mcf_hap5:2617063-2620434 >> $outDir/chr6region.sam
        $samtools view $bam chr6_qbl_hap6:2529777-2533148 >> $outDir/chr6region.sam
        $samtools view $bam chr6_ssto_hap7:2570005-2573397 >> $outDir/chr6region.sam
elif [ $build == "hg38_gencode" ]; then
        echo -n -e "build=hg38_gencode\n"
        $samtools view $bam chr6:29941260-29945884 >> $outDir/chr6region.sam
	$samtools view $bam 6:29941260-29945884 >> $outDir/chr6region.sam
        $samtools view $bam chr6:31353872-31357187 >> $outDir/chr6region.sam
	$samtools view $bam 6:31353872-31357187 >> $outDir/chr6region.sam
        $samtools view $bam chr6:31268749-31272092 >> $outDir/chr6region.sam
	$samtools view $bam 6:31268749-31272092 >> $outDir/chr6region.sam
	$samtools view $bam chr6_GL000250v2_alt:1199010-1203632 >> $outDir/chr6region.sam
	$samtools view $bam chr6_GL000251v2_alt:1420685-1425307 >> $outDir/chr6region.sam
	$samtools view $bam chr6_GL000252v2_alt:1196812-1201452 >> $outDir/chr6region.sam
	$samtools view $bam chr6_GL000253v2_alt:1195801-1200442 >> $outDir/chr6region.sam
	$samtools view $bam chr6_GL000255v2_alt:1196218-1200858 >> $outDir/chr6region.sam
	$samtools view $bam chr6_GL000256v2_alt:1239073-1243714 >> $outDir/chr6region.sam
	$samtools view $bam chr6_GL000251v2_alt:2834232-2837604 >> $outDir/chr6region.sam
	$samtools view $bam chr6_GL000253v2_alt:2662484-2665857 >> $outDir/chr6region.sam
	$samtools view $bam chr6_GL000254v2_alt:2695844-2699230 >> $outDir/chr6region.sam
	$samtools view $bam chr6_GL000255v2_alt:2609569-2612898 >> $outDir/chr6region.sam
	$samtools view $bam chr6_GL000256v2_alt:2656110-2659467 >> $outDir/chr6region.sam
	$samtools view $bam chr6_GL000251v2_alt:2749675-2753056 >> $outDir/chr6region.sam
	$samtools view $bam chr6_GL000252v2_alt:2526549-2529920 >> $outDir/chr6region.sam
	$samtools view $bam chr6_GL000253v2_alt:2577801-2581172 >> $outDir/chr6region.sam
	$samtools view $bam chr6_GL000254v2_alt:2611478-2614849 >> $outDir/chr6region.sam
	$samtools view $bam chr6_GL000255v2_alt:2524181-2527552 >> $outDir/chr6region.sam
	$samtools view $bam chr6_GL000256v2_alt:2570707-2574099 >> $outDir/chr6region.sam
elif [ $build == "hg38_ucsc" ]; then
	echo -n -e "build=hg38_ucsc\n"
	$samtools view $bam chr6:29941260-29945884 >> $outDir/chr6region.sam
        $samtools view $bam 6:29941260-29945884 >> $outDir/chr6region.sam
        $samtools view $bam chr6:31353872-31357187 >> $outDir/chr6region.sam
        $samtools view $bam 6:31353872-31357187 >> $outDir/chr6region.sam
        $samtools view $bam chr6:31268749-31272092 >> $outDir/chr6region.sam
        $samtools view $bam 6:31268749-31272092 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*01:01:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*01:01:01:02N >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*01:01:38L >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*01:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*01:04N >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*01:09 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*01:11N >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*01:14 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*01:16N >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*01:20 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:01:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:01:01:02L >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:01:01:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:01:01:04 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:03:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:03:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:05:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:06:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:07:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:10 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:251 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:259 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:264 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:265 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:266 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:269 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:279 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:32N >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:376 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:43N >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:455 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:48 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:51 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:533 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:53N >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:57 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:60:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:65 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:68 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:77 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:81 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:89 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*02:95 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*03:01:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*03:01:01:02N >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*03:01:01:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*03:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*03:11N >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*03:21N >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*03:36N >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*11:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*11:01:18 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*11:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*11:05 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*11:110 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*11:25 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*11:50Q >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*11:60 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*11:69N >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*11:74 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*11:75 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*11:77 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*23:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*23:09 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*23:38N >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*24:02:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*24:02:01:02L >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*24:02:01:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*24:02:03Q >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*24:02:10 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*24:03:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*24:07:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*24:08 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*24:09N >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*24:10:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*24:11N >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*24:152 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*24:20 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*24:215 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*24:61 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*24:86N >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*25:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*26:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*26:11N >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*26:15 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*26:50 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*29:01:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*29:01:01:02N >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*29:02:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*29:02:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*29:46 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*30:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*30:02:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*30:02:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*30:04:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*30:89 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*31:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*31:01:23 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*31:04 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*31:14N >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*31:46 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*32:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*32:06 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*33:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*33:03:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*33:07 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*34:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*34:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*36:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*43:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*66:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*66:17 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*68:01:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*68:01:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*68:01:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*68:01:02:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*68:02:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*68:02:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*68:02:01:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*68:02:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*68:03:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*68:08:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*68:113 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*68:17 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*68:18N >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*68:22 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*68:71 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*69:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*74:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*74:02:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*74:02:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*80:01:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-A*80:01:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*07:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*07:05:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*07:06 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*07:156 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*07:33:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*07:41 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*07:44 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*07:50 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*08:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*08:08N >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*08:132 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*08:134 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*08:19N >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*08:20 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*08:33 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*08:79 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*13:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*13:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*13:02:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*13:02:09 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*13:08 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*13:15 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*13:25 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*14:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*14:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*14:07N >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*15:01:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*15:01:01:02N >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*15:01:01:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*15:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*15:03:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*15:04:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*15:07:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*15:108 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*15:10:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*15:11:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*15:13:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*15:16:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*15:17:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*15:17:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*15:18:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*15:220 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*15:25:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*15:27:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*15:32:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*15:42 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*15:58 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*15:66 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*15:77 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*15:83 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*18:01:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*18:01:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*18:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*18:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*18:17N >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*18:26 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*18:94N >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*27:04:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*27:05:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*27:05:18 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*27:06 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*27:07:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*27:131 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*27:24 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*27:25 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*27:32 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*35:01:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*35:01:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*35:01:22 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*35:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*35:03:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*35:05:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*35:08:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*35:14:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*35:241 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*35:41 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*37:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*37:01:05 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*38:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*38:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*38:14 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*39:01:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*39:01:01:02L >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*39:01:01:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*39:01:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*39:01:16 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*39:01:21 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*39:05:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*39:06:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*39:10:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*39:13:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*39:14 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*39:34 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*39:38Q >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*40:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*40:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*40:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*40:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*40:06:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*40:06:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*40:10:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*40:150 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*40:40 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*40:72:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*40:79 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*41:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*41:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*42:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*42:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*42:08 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*44:02:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*44:02:01:02S >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*44:02:01:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*44:02:17 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*44:02:27 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*44:03:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*44:03:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*44:04 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*44:09 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*44:138Q >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*44:150 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*44:23N >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*44:26 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*44:46 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*44:49 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*44:56N >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*45:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*45:04 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*46:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*46:01:05 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*47:01:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*47:01:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*48:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*48:03:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*48:04 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*48:08 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*49:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*49:32 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*50:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*51:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*51:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*51:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*51:07:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*51:42 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*52:01:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*52:01:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*52:01:01:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*52:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*53:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*53:11 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*54:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*54:18 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*55:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*55:01:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*55:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*55:12 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*55:24 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*55:48 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*56:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*56:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*56:04 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*57:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*57:03:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*57:06 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*57:11 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*57:29 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*58:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*58:31N >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*59:01:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*59:01:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*67:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*67:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*67:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*73:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*78:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*81:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-B*82:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*01:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*01:02:11 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*01:02:29 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*01:02:30 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*01:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*01:06 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*01:08 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*01:14 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*01:21 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*01:30 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*01:40 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*02:02:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*02:02:02:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*02:10 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*02:11 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*02:16:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*02:69 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*02:85 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*02:86 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*02:87 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*03:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*03:02:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*03:02:02:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*03:02:02:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*03:03:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*03:04:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*03:04:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*03:04:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*03:04:04 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*03:05 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*03:06 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*03:100 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*03:13:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*03:20N >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*03:219 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*03:261 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*03:40:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*03:41:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*03:46 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*03:61 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*04:01:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*04:01:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*04:01:01:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*04:01:01:04 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*04:01:01:05 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*04:01:62 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*04:03:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*04:06 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*04:09N >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*04:128 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*04:161 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*04:177 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*04:70 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*04:71 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*05:01:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*05:01:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*05:08 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*05:09:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*05:93 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*06:02:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*06:02:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*06:02:01:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*06:23 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*06:24 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*06:46N >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:01:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:01:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:01:19 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:01:27 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:01:45 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:02:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:02:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:02:01:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:02:01:04 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:02:01:05 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:02:05 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:02:06 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:02:64 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:04:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:04:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:06 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:149 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:18 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:19 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:26 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:30 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:32N >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:384 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:385 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:386 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:391 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:392 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:49 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:56:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:66 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*07:67 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*08:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*08:01:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*08:02:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*08:02:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*08:03:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*08:04:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*08:112 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*08:20 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*08:21 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*08:22 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*08:24 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*08:27 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*08:36N >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*08:40 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*08:41 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*08:62 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*12:02:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*12:03:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*12:03:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*12:08 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*12:13 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*12:19 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*12:22 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*12:99 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*14:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*14:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*14:21N >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*14:23 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*15:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*15:05:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*15:05:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*15:13 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*15:16 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*15:17 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*15:96Q >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*16:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*16:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*16:04:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*17:01:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*17:01:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*17:01:01:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*17:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-C*18:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQA1*01:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQA1*01:02:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQA1*01:02:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQA1*01:02:01:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQA1*01:02:01:04 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQA1*01:03:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQA1*01:03:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQA1*01:04:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQA1*01:04:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQA1*01:05:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQA1*01:07 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQA1*01:10 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQA1*01:11 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQA1*02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQA1*03:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQA1*03:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQA1*03:03:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQA1*04:01:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQA1*04:01:02:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQA1*04:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQA1*05:01:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQA1*05:01:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQA1*05:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQA1*05:05:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQA1*05:05:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQA1*05:05:01:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQA1*05:11 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQA1*06:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQB1*02:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQB1*02:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQB1*03:01:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQB1*03:01:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQB1*03:01:01:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQB1*03:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQB1*03:03:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQB1*03:03:02:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQB1*03:03:02:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQB1*03:05:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQB1*05:01:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQB1*05:01:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQB1*05:03:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQB1*05:03:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQB1*06:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQB1*06:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQB1*06:03:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DQB1*06:09:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DRB1*01:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DRB1*01:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DRB1*03:01:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DRB1*03:01:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DRB1*04:03:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DRB1*07:01:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DRB1*07:01:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DRB1*08:03:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DRB1*09:21 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DRB1*10:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DRB1*11:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DRB1*11:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DRB1*11:04:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DRB1*12:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DRB1*12:17 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DRB1*13:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DRB1*13:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DRB1*14:05:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DRB1*14:54:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DRB1*15:01:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DRB1*15:01:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DRB1*15:01:01:03 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DRB1*15:01:01:04 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DRB1*15:02:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DRB1*15:03:01:01 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DRB1*15:03:01:02 >> $outDir/chr6region.sam
	$samtools view $bam HLA-DRB1*16:02:01 >> $outDir/chr6region.sam
else
	echo -n -e "undefind refenerce\n"
	exit 2
fi

$samtools view -bS -o $outDir/chr6region.bam $outDir/chr6region.sam

$bam2fastq -o $outDir/chr6region#.fastq --aligned -q $outDir/chr6region.bam  # LWH 2017-06-28
mv $outDir/chr6region_1.fastq $outDir/chr6region.1.fastq  # LWH 2017-06-28
mv $outDir/chr6region_2.fastq $outDir/chr6region.2.fastq  # LWH 2017-06-28

$PSHOME/bin/clean_unpaired_fastq.pl $outDir/chr6region.1.fastq
$PSHOME/bin/clean_unpaired_fastq.pl $outDir/chr6region.2.fastq
gzip $outDir/chr6region.1.fastq
gzip $outDir/chr6region.2.fastq

# merge the two sets of fastqs

echo -n -e "merge the two sets of fastqs\n"

$PSHOME/bin/perl $PSHOME/bin/merge_fastq.pl $outDir/tag $outDir/chr6region $outDir/merged

rm -f $outDir/*sam

# alignment

echo -n -e "aligning to HLA library\n"

cat $PSHOME/files/novoalign_complete_header.sam > $outDir/nv.complete.chr6region.R0k6.sam

$PSHOME/bin/perl $PSHOME/bin/align_fork_fh.pl $outDir/merged.1.fastq $outDir/merged.2.fastq $NUM_THREADS $format $PSHOME/files/abc_complete.nix $outDir/nv.complete.chr6region.R0k6.sam 0 $NOVOALIGN_DIR

$samtools view -bS -o $outDir/nv.complete.chr6region.R0k6.bam $outDir/nv.complete.chr6region.R0k6.sam

rm -f $outDir/nv.complete.chr6region.R0k6.sam

echo -n -e "sorting\n"

$gatk --java-options "-Xmx5G -XX:+PrintCommandLineFlags -XX:ParallelGCThreads=1" SortSam \
--TMP_DIR $TMP_DIR \
--INPUT $outDir/nv.complete.chr6region.R0k6.bam \
--OUTPUT $outDir/nv.complete.chr6region.R0k6.csorted.bam \
--VALIDATION_STRINGENCY SILENT \
--SORT_ORDER coordinate \
--CREATE_INDEX true

t=`ls -lh $outDir/nv.complete.chr6region.R0k6.csorted.bam`
echo -n -e "size of bam = $t\n"

rm -f $outDir/nv.complete.chr6region.R0k6.bam

# first allele calculations

echo -n -e "calculating lik1\n"
date

$PSHOME/bin/perl $PSHOME/bin/first_allele_calculations_fork.pl $ids $PSHOME $samtools $NUM_THREADS $race $iFile $outDir $outDir/nv.complete.chr6region.R0k6.csorted.bam $PSHOME/bin/perl

date

echo -n -e "get first winners\n"
date

rm $outDir/counts1.R0k6

for i in $var; do
        a=`tail -1 $outDir/$i.lik1 | cut -f2`
        echo -n -e "$i\t$a\n" >> $outDir/counts1.R0k6
done

winner1_a=`sort -k2,2rn $outDir/counts1.R0k6 | grep hla_a | head -1 | cut -f1`
winner1_b=`sort -k2,2rn $outDir/counts1.R0k6 | grep hla_b | head -1 | cut -f1`
winner1_c=`sort -k2,2rn $outDir/counts1.R0k6 | grep hla_c | head -1 | cut -f1`

date


# second allele calculations

echo -n -e "calculating lik2\n"

$PSHOME/bin/second_allele_calculations.pl $race $outDir/counts1.R0k6 $ids 1 $PSHOME $outDir

date

echo -n -e "get second winners\n"

rm $outDir/counts2.R0k6

for i in $var; do
        a=`tail -1 $outDir/$i.lik2 | cut -f2`
        echo -n -e "$i\t$a\n" >> $outDir/counts2.R0k6
done

winner2_a=`sort -k2,2rn $outDir/counts2.R0k6 | grep hla_a | head -1 | cut -f1`
winner2_b=`sort -k2,2rn $outDir/counts2.R0k6 | grep hla_b | head -1 | cut -f1`
winner2_c=`sort -k2,2rn $outDir/counts2.R0k6 | grep hla_c | head -1 | cut -f1`

echo -n -e "winners1\t$winner1_a\t$winner1_b\t$winner1_c\n"
echo -n -e "winners2\t$winner2_a\t$winner2_b\t$winner2_c\n"


echo -n -e "HLA-A\t$winner1_a\t$winner2_a\nHLA-B\t$winner1_b\t$winner2_b\nHLA-C\t$winner1_c\t$winner2_c\n" > $outDir/winners.hla.txt

$PSHOME/bin/perl $PSHOME/bin/change_hla_format.pl $outDir/winners.hla.txt $outDir/winners.hla.list
# cleanup

echo -n -e "cleanup\n"

rm -f $outDir/temp* $outDir/*sam $outDir/tag*fastq0* $outDir/ids_*
rm -rf $outDir/*REF* $outDir/tmp
#rm -rf $outDir/nv*
rm -rf $outDir/*fastq
rm -rf $outDir/tag*
rm -rf $outDir/chr6region*
rm -rf $outDir/merged*
#rm -rf $outDir/counts*
rm -f $outDir/*lik1 $outDir/*lik2

