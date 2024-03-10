#========================================================
# obtain clearn data from raw fastq files.
# 
conda activate sunbeam3
sunbeam init /home/dell/dataanalysis/projects/longevit/sunbeam_results \
--data_fp /home/dell/dataanalysis/projects/longevit/fastq_datasets  --force
# edit yml file with add path of reference database
cd sunbeam_results
vi sunbeam_config.yml
# add path for each database
host_fp: '/home/dell/dataanalysis/pipelines/databases/sunbeam/refernce_db/sunbeam.host.genome'
# qc _decom

sunbeam run --configfile ~/dataanalysis/projects/longevit/sunbeam_results/sunbeam_config.yml all_assembly
conda deactivate 
#========================================================
# Classification 
conda activate kraken_pipeline
#1. For reads 
for file in /home/dell/dataanalysis/projects/longevit/sunbeam_results/sunbeam_output/qc/decontam/*_1.fastq.gz
do
	file1=$file
	file2=`echo $file1 | awk -F "_1" '{print $1 "_2" $2}'`
    output=`echo $file | awk -F "_1.fastq.gz" '{print $1 "_k2output.txt"}'`
    report=`echo $file | awk -F "_1.fastq.gz" '{print $1 "_k2report.txt"}'`
    kraken2 --db /home/dell/dataanalysis/pipelines/databases/karken_database_20220607 --threads 12 --use-names \
--minimum-hit-groups 3 --report-minimizer-data --paired $file1 $file2 \
--output $output \
--report $report
    bracken_results=`echo $file | awk -F "_1.fastq.gz" '{print $1 "_d_k2_brackentaxa.tsv"}'`
    bracken -d /home/dell/dataanalysis/pipelines/databases/karken_database_20220607 \
    -i $report \
    -o $bracken_results -l S

done

conda deactivate 
