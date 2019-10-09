HLA分型预测是通过病人的测序数据（bam/cram文件为输入）计算其HLA分型的模块，HLA分型是新抗原预测模块的输入文件；
该模块的使用示例为(work.sh是一个使用实例)：
/zfssz5/BC_PS/chengyuanfang/work/CGIA_3.0_module/HLA_prediction/call_hla_type.sh SAMPLE.bam Asian hg19 STDFQ OUTDIR
参数说明:
SAMPLE.bam——样本比对后的bam文件,或者cram文件;
Asian——人种，可选Caucasian、Black、Asian和Unknown;
hg19——参考基因组版本，可选hg19/hg38_ucsc/hg38_gencode;
STDFQ——novoalign的参数，设置为STDFQ即可;
OUTDIR——输出目录，不必提前建立,输出结果为OUTDIR/winners.hla.list;

说明：
1.OUTDIR不用提前建立，在运行时候流程会新建这个文件夹；
2.推荐使用Asian hg19 STDFQ这三个参数；
3.对于肿瘤配对样本，可以同时分析一下normal.bam和tumor.bam,检查一下是否发生了HLA-LOH；
4.该流程以bam文件为输入，因此monitor描述依赖关系时可以放在bam文件之后，和变异检测那些模块并行，新抗原模块串行在其后，节省时间；

