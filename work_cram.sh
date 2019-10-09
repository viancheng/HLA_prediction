#!/bin/bash
set -e
set -o pipefail
echo hostname: `hostname`
echo ==========start at : `date` ==========

echo `date` > /zfssz5/BC_PS/chengyuanfang/work/CGIA_3.0_module/HLA_prediction/work_cram.sh.sign && \
/zfssz5/BC_PS/chengyuanfang/work/CGIA_3.0_module/HLA_prediction/call_hla_type.sh /zfssz5/BC_PS/chengyuanfang/work/CGIA_3.0_module/HLA_prediction/data/BGI-0107-ESCC-014T.cram Asian hg19 STDFQ /zfssz5/BC_PS/chengyuanfang/work/CGIA_3.0_module/HLA_prediction/output_cram && \
echo `date` >> /zfssz5/BC_PS/chengyuanfang/work/CGIA_3.0_module/HLA_prediction/work_cram.sh.sign && \

if [ -e /zfssz5/BC_PS/chengyuanfang/work/CGIA_3.0_module/HLA_prediction/work_cram.sh.e* ]; then ls /zfssz5/BC_PS/chengyuanfang/work/CGIA_3.0_module/HLA_prediction/work_cram.sh.e* | while read dd ; do jobid=`echo $dd | sed -r 's/^.*\.e([0-9]*)$/\1/g'`; echo $jobid; qstat -j $jobid | grep 'usage'; done ; fi && \
echo ==========end at : `date` ==========
