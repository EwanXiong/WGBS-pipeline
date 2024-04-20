#######===============================================================dir_set===========================================================######
###另外注意更改相关性计算时的样本数，即1的个数
file_name=name.txt
group=group.txt
input=rawdata
output=output
genesome=hg19/bowtie2/
ref_gene=hg19
mlyth_ref=hg19_CpGi_bed
anno_ref=hg19_RefSeq_bed12
dbsnp=BisSNP/dbsnp_135.hg19.sort.vcf
snp=BisSNP/1000G_phase1.snps.high_confidence.hg19.sites.vcf
Mills=BisSNP/Mills_and_1000G_gold_standard.indels.hg19.sites.sort.vcf
CNVkit_ref1=CNVkit/refFlat_hg19.txt
CNVkit_ref2=CNVkit/access-5kb-mappable.hg19.bed

echo "[`date`] ======================================================WGBS START================================================================"
name=($(cat $file_name))
group=$(cat $group)
if [ ! -d "$output/1_qc" ]; then
		echo
		echo 
		echo "[`date`] qc"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/1_qc
		myvar=0
		for i in ${name[@]}
		do
			time /opt/tsinghua/software/fastqc/fastqc/FastQC/fastqc \
			--outdir $output/1_qc  \
			--threads 20  $input/"$i"* &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "20" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		multiqc $output/1_qc -o $output/1_qc/
		echo '[`date`] Run complete'
		echo '------------------------------------------------------------------------------------------------------------------------'
fi

if [ ! -d "$output/2_low_complexity_filter" ]; then
		echo
		echo 
		echo "[`date`] low_complexity_filter"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/2_low_complexity_filter
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/2_low_complexity_filter/"$i"
			time fastp \
			--thread 30 \
			--correction \
			--detect_adapter_for_pe \
			--json $output/2_low_complexity_filter/"$i"/"$i".json \
			--html $output/2_low_complexity_filter/"$i"/"$i".html \
			--report_title "$i" \
			-i $input/"$i"_1.fq.gz \
			-o $output/2_low_complexity_filter/"$i"/"$i".pair1.truncated.gz \
			-I $input/"$i"_2.fq.gz \
			-O $output/2_low_complexity_filter/"$i"/"$i".pair2.truncated.gz &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "20" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo '[`date`] Run complete'
		echo '------------------------------------------------------------------------------------------------------------------------'
fi

if [ ! -d "$output/3_bismark" ]; then
		echo
		echo 
		echo "[`date`] bismark"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/3_bismark
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/3_bismark/"$i"
			time bismark \
			-q \
			--phred33-quals \
			--parallel 20 \
			--bowtie2 \
			--un \
			--multicore 4 \
			--output_dir $output/3_bismark/"$i" \
			--temp_dir /opt/intermediata \
			--genome_folder $genesome \
			-1 $output/2_low_complexity_filter/"$i"/"$i".pair1.truncated.gz \
			-2 $output/2_low_complexity_filter/"$i"/"$i".pair2.truncated.gz \
			> $output/3_bismark/"$i"/step3.log 2>&1 &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "8" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo '[`date`] Run complete'
		echo '------------------------------------------------------------------------------------------------------------------------'
fi

if [ ! -d "$output/4_deduplicate_bismark" ]; then
		echo
		echo 
		echo "[`date`] deduplicate_bismark"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/4_deduplicate_bismark
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/4_deduplicate_bismark/"$i"
			time deduplicate_bismark \
			--output_dir $output/4_deduplicate_bismark/"$i" \
			--paired \
			--bam $output/3_bismark/"$i"/"$i".pair1.truncated.gz_bismark_bt2_pe.bam &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "30" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo '[`date`] Run complete'
		echo '------------------------------------------------------------------------------------------------------------------------'
fi

####此步骤提取甲基化信息

if [ ! -d "$output/5_bismark_methylation_extractor" ]; then
		echo
		echo 
		echo "[`date`] bismark_methylation_extractor"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/5_bismark_methylation_extractor
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/5_bismark_methylation_extractor/"$i"
			time bismark_methylation_extractor \
			--gzip \
			--bedGraph \
			--buffer_size 10G \
			-p \
			--parallel 30 \
			--output $output/5_bismark_methylation_extractor/"$i" \
			$output/4_deduplicate_bismark/"$i"/"$i".pair1.truncated.gz_bismark_bt2_pe.deduplicated.bam &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "10" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo '[`date`] Run complete'
		echo '------------------------------------------------------------------------------------------------------------------------'
fi

if [ ! -d "$output/6_bismark_methylation_report" ]; then
		echo
		echo 
		echo "[`date`] bismark_methylation_report"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/6_bismark_methylation_report
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/6_bismark_methylation_report/"$i"
			time bismark2report \
			--alignment_report $output/3_bismark/"$i"/"$i".pair1.truncated.gz_bismark_bt2_PE_report.txt \
			--dir $output/6_bismark_methylation_report/"$i"  &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "20" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo '[`date`] Run complete'
		echo '------------------------------------------------------------------------------------------------------------------------'
fi

if [ ! -d "$output/summary" ]; then
		echo
		echo 
		echo "[`date`] summary"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/summary
		mkdir -p $output/summary/bismark
		mkdir -p $output/summary/deduplicate_bismark
		myvar=0
		for i in ${name[@]}
		do
			time cp $output/3_bismark/"$i"/"$i".pair1.truncated.gz_bismark_bt2_PE_report.txt $output/summary/bismark &
			time cp $output/4_deduplicate_bismark/"$i"/"$i".pair1.truncated.gz_bismark_bt2_pe.deduplication_report.txt $output/summary/deduplicate_bismark &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "70" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		multiqc $output/summary/bismark/*report.txt -o  $output/summary/bismark ;time rm -rf $output/summary/bismark/*.txt &
		multiqc $output/summary/deduplicate_bismark/*report.txt -o  $output/summary/deduplicate_bismark ;time rm -rf $output/summary/deduplicate_bismark/*.txt&
		echo '[`date`] Run complete'
		echo '------------------------------------------------------------------------------------------------------------------------'
fi

if [ ! -d "$output/7_bamsort" ]; then
		echo
		echo 
		echo "[`date`] bamsort"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/7_bamsort
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/7_bamsort/"$i"
			time samtools sort \
			-@ 18 \
			-o $output/7_bamsort/"$i"/"$i".pair1_sorted.bam \
			$output/4_deduplicate_bismark/"$i"/"$i".pair1.truncated.gz_bismark_bt2_pe.deduplicated.bam &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "30" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		myvar=0
		for i in ${name[@]}
		do
			time samtools index -@ 18 $output/7_bamsort/"$i"/"$i".pair1_sorted.bam &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "30" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo '[`date`] Run complete'
		echo '------------------------------------------------------------------------------------------------------------------------'
fi

if [ ! -d "$output/8_qualimap" ]; then
		echo
		echo 
		echo "[`date`] qualimap"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/8_qualimap
		for i in ${name[@]}
		do
			echo $output/7_bamsort/"$i"/"$i".pair1_sorted.bam >>$output/8_qualimap/site.txt
		done
		wait
		paste /opt/tsinghua/NuoHe/WGBS_22.10.8/name.txt $output/8_qualimap/site.txt >$output/8_qualimap/ref.txt
		wait
		time qualimap multi-bamqc \
		-r \
		-d $output/8_qualimap/ref.txt -outdir $output/8_qualimap/ \
		-outformat PDF:HTML \
		--java-mem-size=20G 
		echo '[`date`] Run complete'
		echo '------------------------------------------------------------------------------------------------------------------------'
fi

if [ ! -d "$output/9_cpg_collect" ]; then
		echo
		echo 
		echo "[`date`] cpg_collect"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/9_cpg_collect
		myvar=0
		for i in ${name[@]}
		do
			time Rscript /root/anaconda3/envs/cfDNApipe2.0/lib/python3.6/site-packages/cfDNApipe/data/Rscript/methylKit_processBam.r \
			--bam $output/7_bamsort/"$i"/"$i".pair1_sorted.bam --id "$i" \
			--outputdir $output/9_cpg_collect \
			--assembly $ref_gene \
			--cx CpG CHG CHH \
			--minqual 20 \
			--mincov 0 &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "30" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		mkdir -p $output/summary/Bis_conver
		for i in ${name[@]}
		do
			cat $output/3_bismark/"$i"/"$i".pair1.truncated.gz_bismark_bt2_PE_report.txt|grep 'Total methylated'|awk -F '\t' '{print$2}'|sed '4d'|awk '{sum+=$1} END {print sum}' >>$output/summary/Bis_conver/mCG.txt ;\
			cat $output/3_bismark/"$i"/"$i".pair1.truncated.gz_bismark_bt2_PE_report.txt|grep 'Total methylated'|awk -F '\t' '{print$2}'|awk '{sum+=$1} END {print sum}' >>$output/summary/Bis_conver/mC.txt
		done
		wait
		paste -d '\t' $file_name $output/summary/Bis_conver/mCG.txt  $output/summary/Bis_conver/mC.txt >> $output/summary/Bis_conver/bis.txt
		echo '[`date`] Run complete'
		echo '------------------------------------------------------------------------------------------------------------------------'
fi

# treatment样本数需要改，几个样本几个1
cpg_file=($(ls $output/9_cpg_collect/*CpG.txt))

if [ ! -d "$output/10_correlation" ]; then
		echo
		echo 
		echo "[`date`] correlation"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/10_correlation
		time Rscript /root/anaconda3/envs/cfDNApipe2.0/lib/python3.6/site-packages/cfDNApipe/data/Rscript/methylKit_summary_diff.r \
		--files "${cpg_file[@]}" \
		--ids "${name[@]}" \
		--treatment 1 1  \
		--assembly $ref_gene \
		--context CpG \
		--cores 18 \
		-o $output/10_correlation &
		wait
		echo '[`date`] Run complete'
		echo '------------------------------------------------------------------------------------------------------------------------'
fi

if [ ! -d "$output/11_methylLevel" ]; then
		echo
		echo 
		echo "[`date`] methylLevel"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/11_methylLevel
		
		time Rscript /root/anaconda3/envs/cfDNApipe2.0/lib/python3.6/site-packages/cfDNApipe/data/Rscript/meth_levels.r \
		--files "${cpg_file[@]}" \
		--ids "${name[@]}"  \
		--treatment 1 1 \
		--assembly $ref_gene \
		--context CpG \
		--func regions_bins regions_all bed \
		-o $output/11_methylLevel \
		--mincov 10 \
		--num 20 \
		--len 200 \
		--window 10000 \
		--bed_CpGi $mlyth_ref &
		wait
		echo '[`date`] Run complete'
		echo '------------------------------------------------------------------------------------------------------------------------'
fi

if [ ! -d "$output/12_methylstat" ]; then
		echo
		echo 
		echo "[`date`] methylstat"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/12_methylstat
		time  Rscript /root/anaconda3/envs/cfDNApipe2.0/lib/python3.6/site-packages/cfDNApipe/data/Rscript/methylKit_methylstat.r \
		--files "${cpg_file[@]}" \
		--ids "${name[@]}" \
		--treatment 1 1 \
		--assembly $ref_gene \
		--context CpG \
		-o $output/12_methylstat &
		wait
		echo '[`date`] Run complete'
		echo '------------------------------------------------------------------------------------------------------------------------'
fi

# treatment 分组重复：1 0 == 1样本1对照； 1 1 0 0 == 2样本2对照

if [ ! -d "$output/13_DMR" ]; then
		echo
		echo 
		echo "[`date`] DMR"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/13_DMR
		myvar=0
		for i in ${group[@]}
		do
			mkdir -p $output/13_DMR/"$i"
			time Rscript /root/anaconda3/envs/cfDNApipe2.0/lib/python3.6/site-packages/cfDNApipe/data/Rscript/methylKit_diffMethyl.r \
			--files $output/9_cpg_collect/"${i%*_*}"_CpG.txt $output/9_cpg_collect/"${i#*_*}"_CpG.txt \
			--ids "${i%*_*}" "${i#*_*}" \
			--treatment 1 0 \
			--assembly $ref_gene \
			--context CpG \
			--cores 18 \
			-o $output/13_DMR/"$i" &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "30" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		for i in ${group[@]}
		do
			time Rscript /root/anaconda3/envs/cfDNApipe2.0/lib/python3.6/site-packages/cfDNApipe/data/Rscript/wgbs_circos.r \
			--input $output/13_DMR/"$i"/DMR_all.txt \
			--output $output/13_DMR/"$i"/circlize.pdf &
		done
		wait
		echo '[`date`] Run complete'
		echo '------------------------------------------------------------------------------------------------------------------------'
fi

if [ ! -d "$output/14_DMR_annotation" ]; then
		echo
		echo 
		echo "[`date`] DMR_annotation"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/14_DMR_annotation
				myvar=0
		for i in ${group[@]}
		do
			mkdir -p $output/14_DMR_annotation/"$i"
			time Rscript /root/anaconda3/envs/cfDNApipe2.0/lib/python3.6/site-packages/cfDNApipe/data/Rscript/methyl_annotation.r \
			-file $output/13_DMR/"$i"/DMR_all.txt \
			-bed12 $anno_ref \
			-CpGi $mlyth_ref \
			-o $output/14_DMR_annotation/"$i"/ \
			--name DMR_all \
			--upFlank 1000 \
			--downFlank 1000 \
			--flankCpGi 2000 \
			--sep tab &
			time Rscript /root/anaconda3/envs/cfDNApipe2.0/lib/python3.6/site-packages/cfDNApipe/data/Rscript/methyl_annotation.r \
			-file $output/13_DMR/"$i"/DMC_all.txt \
			-bed12 $anno_ref \
			-CpGi $mlyth_ref \
			-o $output/14_DMR_annotation/"$i"/ \
			--name DMC_all \
			--upFlank 1000 \
			--downFlank 1000 \
			--flankCpGi 2000 \
			--sep tab &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "30" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		for i in ${group[@]}
		do
			awk 'BEGIN{FS=","} {print $10}' $output/14_DMR_annotation/"$i"/DMR_all_annotation.csv  |  sort -u > $output/14_DMR_annotation/"$i"/DMR_gene.txt ;\
			awk 'BEGIN{FS=","} {print $10}' $output/14_DMR_annotation/"$i"/DMC_all_annotation.csv  |  sort -u > $output/14_DMR_annotation/"$i"/DMC_gene.txt ;\
			time Rscript /root/anaconda3/envs/cfDNApipe2.0/lib/python3.6/site-packages/cfDNApipe/data/Rscript/NC_ID_convert.r \
			--input $output/14_DMR_annotation/"$i"/DMC_gene.txt \
			--output $output/14_DMR_annotation/"$i" &
		done
		echo '[`date`] Run complete'
		echo '------------------------------------------------------------------------------------------------------------------------'
fi

echo "[`date`] ======================================================WGBS SNP START================================================================"

if [ ! -d "$output/16_Add_RG" ]; then
		echo
		echo 
		echo "[`date`] Add_RG"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/16_Add_RG
		myvar=0
		for i in ${name[@]}
		do
			time gatk --java-options ""-Xmx6G"" \
			AddOrReplaceReadGroups \
			--INPUT $output/7_bamsort/"$i"/"$i".pair1_sorted.bam \
			--OUTPUT $output/16_Add_RG/"$i".clean-RG.bam \
			--SORT_ORDER coordinate \
			--RGID "$i" \
			--RGLB "$i" \
			--RGPL ILLUMINA \
			--RGSM "$i" \
			--RGPU Null \
			--CREATE_INDEX &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "20" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo '[`date`] Run complete'
		echo '------------------------------------------------------------------------------------------------------------------------'
fi

if [ ! -d "$output/17_BisulfiteCountCovariates" ]; then
		echo
		echo 
		echo "[`date`] BisulfiteCountCovariates"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/17_BisulfiteCountCovariates
		myvar=0
		for i in ${name[@]}
		do
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-R $genesome/"$ref_gene".fa \
			-T BisulfiteCountCovariates \
			-I $output/16_Add_RG/"$i".clean-RG.bam \
			--knownSites $dbsnp \
			--knownSites $snp \
			--knownSites $Mills \
			-cov ReadGroupCovariate \
			-cov QualityScoreCovariate \
			-cov CycleCovariate \
			-recalFile $output/17_BisulfiteCountCovariates/"$i".recalFile_before.csv \
			-nt 1 &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "50" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo '[`date`] Run complete'
		echo '------------------------------------------------------------------------------------------------------------------------'
fi

if [ ! -d "$output/18_BisulfiteTableRecalibration" ]; then
		echo
		echo 
		echo "[`date`] BisulfiteTableRecalibration"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/18_BisulfiteTableRecalibration
		myvar=0
		for i in ${name[@]}
		do
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T BisulfiteTableRecalibration \
			-R $genesome/"$ref_gene".fa \
			-I $output/16_Add_RG/"$i".clean-RG.bam \
			-recalFile $output/17_BisulfiteCountCovariates/"$i".recalFile_before.csv \
			-o $output/18_BisulfiteTableRecalibration/"$i"-BQSR.bam &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "30" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo '[`date`] Run complete'
		echo '------------------------------------------------------------------------------------------------------------------------'
fi 

if [ ! -d "$output/19_BisulfiteGenotyper" ]; then
		echo
		echo 
		echo "[`date`] BisulfiteGenotyper"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/19_BisulfiteGenotyper
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/19_BisulfiteGenotyper/"$i"
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T BisulfiteGenotyper \
			-R $genesome/"$ref_gene".fa \
			-I $output/18_BisulfiteTableRecalibration/"$i"-BQSR.bam \
			-D $dbsnp \
			-L chr1 \
			-vfn1 $output/19_BisulfiteGenotyper/"$i"/"$i"_chr1.snp.raw.vcf \
			-out_modes EMIT_VARIANTS_ONLY \
			-stand_call_conf 10 \
			-stand_emit_conf 0 \
			-nt 1 \
			-mmq 20 \
			-mbq 5 &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T BisulfiteGenotyper \
			-R $genesome/"$ref_gene".fa \
			-I $output/18_BisulfiteTableRecalibration/"$i"-BQSR.bam \
			-D $dbsnp \
			-L chr2 \
			-vfn1 $output/19_BisulfiteGenotyper/"$i"/"$i"_chr2.snp.raw.vcf \
			-out_modes EMIT_VARIANTS_ONLY \
			-stand_call_conf 10 \
			-stand_emit_conf 0 \
			-nt 1 \
			-mmq 20 \
			-mbq 5 &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T BisulfiteGenotyper \
			-R $genesome/"$ref_gene".fa \
			-I $output/18_BisulfiteTableRecalibration/"$i"-BQSR.bam \
			-D $dbsnp \
			-L chr3 \
			-vfn1 $output/19_BisulfiteGenotyper/"$i"/"$i"_chr3.snp.raw.vcf \
			-out_modes EMIT_VARIANTS_ONLY \
			-stand_call_conf 10 \
			-stand_emit_conf 0 \
			-nt 1 \
			-mmq 20 \
			-mbq 5 &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T BisulfiteGenotyper \
			-R $genesome/"$ref_gene".fa \
			-I $output/18_BisulfiteTableRecalibration/"$i"-BQSR.bam \
			-D $dbsnp \
			-L chr4 \
			-vfn1 $output/19_BisulfiteGenotyper/"$i"/"$i"_chr4.snp.raw.vcf \
			-out_modes EMIT_VARIANTS_ONLY \
			-stand_call_conf 10 \
			-stand_emit_conf 0 \
			-nt 1 \
			-mmq 20 \
			-mbq 5 &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T BisulfiteGenotyper \
			-R $genesome/"$ref_gene".fa \
			-I $output/18_BisulfiteTableRecalibration/"$i"-BQSR.bam \
			-D $dbsnp \
			-L chr5 \
			-vfn1 $output/19_BisulfiteGenotyper/"$i"/"$i"_chr5.snp.raw.vcf \
			-out_modes EMIT_VARIANTS_ONLY \
			-stand_call_conf 10 \
			-stand_emit_conf 0 \
			-nt 1 \
			-mmq 20 \
			-mbq 5 &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T BisulfiteGenotyper \
			-R $genesome/"$ref_gene".fa \
			-I $output/18_BisulfiteTableRecalibration/"$i"-BQSR.bam \
			-D $dbsnp \
			-L chr6 \
			-vfn1 $output/19_BisulfiteGenotyper/"$i"/"$i"_chr6.snp.raw.vcf \
			-out_modes EMIT_VARIANTS_ONLY \
			-stand_call_conf 10 \
			-stand_emit_conf 0 \
			-nt 1 \
			-mmq 20 \
			-mbq 5 &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T BisulfiteGenotyper \
			-R $genesome/"$ref_gene".fa \
			-I $output/18_BisulfiteTableRecalibration/"$i"-BQSR.bam \
			-D $dbsnp \
			-L chr7 \
			-vfn1 $output/19_BisulfiteGenotyper/"$i"/"$i"_chr7.snp.raw.vcf \
			-out_modes EMIT_VARIANTS_ONLY \
			-stand_call_conf 10 \
			-stand_emit_conf 0 \
			-nt 1 \
			-mmq 20 \
			-mbq 5 &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T BisulfiteGenotyper \
			-R $genesome/"$ref_gene".fa \
			-I $output/18_BisulfiteTableRecalibration/"$i"-BQSR.bam \
			-D $dbsnp \
			-L chr8 \
			-vfn1 $output/19_BisulfiteGenotyper/"$i"/"$i"_chr8.snp.raw.vcf \
			-out_modes EMIT_VARIANTS_ONLY \
			-stand_call_conf 10 \
			-stand_emit_conf 0 \
			-nt 1 \
			-mmq 20 \
			-mbq 5 &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T BisulfiteGenotyper \
			-R $genesome/"$ref_gene".fa \
			-I $output/18_BisulfiteTableRecalibration/"$i"-BQSR.bam \
			-D $dbsnp \
			-L chr9 \
			-vfn1 $output/19_BisulfiteGenotyper/"$i"/"$i"_chr9.snp.raw.vcf \
			-out_modes EMIT_VARIANTS_ONLY \
			-stand_call_conf 10 \
			-stand_emit_conf 0 \
			-nt 1 \
			-mmq 20 \
			-mbq 5 &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T BisulfiteGenotyper \
			-R $genesome/"$ref_gene".fa \
			-I $output/18_BisulfiteTableRecalibration/"$i"-BQSR.bam \
			-D $dbsnp \
			-L chr10 \
			-vfn1 $output/19_BisulfiteGenotyper/"$i"/"$i"_chr10.snp.raw.vcf \
			-out_modes EMIT_VARIANTS_ONLY \
			-stand_call_conf 10 \
			-stand_emit_conf 0 \
			-nt 1 \
			-mmq 20 \
			-mbq 5 &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T BisulfiteGenotyper \
			-R $genesome/"$ref_gene".fa \
			-I $output/18_BisulfiteTableRecalibration/"$i"-BQSR.bam \
			-D $dbsnp \
			-L chr11 \
			-vfn1 $output/19_BisulfiteGenotyper/"$i"/"$i"_chr11.snp.raw.vcf \
			-out_modes EMIT_VARIANTS_ONLY \
			-stand_call_conf 10 \
			-stand_emit_conf 0 \
			-nt 1 \
			-mmq 20 \
			-mbq 5 &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T BisulfiteGenotyper \
			-R $genesome/"$ref_gene".fa \
			-I $output/18_BisulfiteTableRecalibration/"$i"-BQSR.bam \
			-D $dbsnp \
			-L chr12 \
			-vfn1 $output/19_BisulfiteGenotyper/"$i"/"$i"_chr12.snp.raw.vcf \
			-out_modes EMIT_VARIANTS_ONLY \
			-stand_call_conf 10 \
			-stand_emit_conf 0 \
			-nt 1 \
			-mmq 20 \
			-mbq 5 &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T BisulfiteGenotyper \
			-R $genesome/"$ref_gene".fa \
			-I $output/18_BisulfiteTableRecalibration/"$i"-BQSR.bam \
			-D $dbsnp \
			-L chr13 \
			-vfn1 $output/19_BisulfiteGenotyper/"$i"/"$i"_chr13.snp.raw.vcf \
			-out_modes EMIT_VARIANTS_ONLY \
			-stand_call_conf 10 \
			-stand_emit_conf 0 \
			-nt 1 \
			-mmq 20 \
			-mbq 5 &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T BisulfiteGenotyper \
			-R $genesome/"$ref_gene".fa \
			-I $output/18_BisulfiteTableRecalibration/"$i"-BQSR.bam \
			-D $dbsnp \
			-L chr14 \
			-vfn1 $output/19_BisulfiteGenotyper/"$i"/"$i"_chr14.snp.raw.vcf \
			-out_modes EMIT_VARIANTS_ONLY \
			-stand_call_conf 10 \
			-stand_emit_conf 0 \
			-nt 1 \
			-mmq 20 \
			-mbq 5 &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T BisulfiteGenotyper \
			-R $genesome/"$ref_gene".fa \
			-I $output/18_BisulfiteTableRecalibration/"$i"-BQSR.bam \
			-D $dbsnp \
			-L chr15 \
			-vfn1 $output/19_BisulfiteGenotyper/"$i"/"$i"_chr15.snp.raw.vcf \
			-out_modes EMIT_VARIANTS_ONLY \
			-stand_call_conf 10 \
			-stand_emit_conf 0 \
			-nt 1 \
			-mmq 20 \
			-mbq 5 &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T BisulfiteGenotyper \
			-R $genesome/"$ref_gene".fa \
			-I $output/18_BisulfiteTableRecalibration/"$i"-BQSR.bam \
			-D $dbsnp \
			-L chr16 \
			-vfn1 $output/19_BisulfiteGenotyper/"$i"/"$i"_chr16.snp.raw.vcf \
			-out_modes EMIT_VARIANTS_ONLY \
			-stand_call_conf 10 \
			-stand_emit_conf 0 \
			-nt 1 \
			-mmq 20 \
			-mbq 5 &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T BisulfiteGenotyper \
			-R $genesome/"$ref_gene".fa \
			-I $output/18_BisulfiteTableRecalibration/"$i"-BQSR.bam \
			-D $dbsnp \
			-L chr17 \
			-vfn1 $output/19_BisulfiteGenotyper/"$i"/"$i"_chr17.snp.raw.vcf \
			-out_modes EMIT_VARIANTS_ONLY \
			-stand_call_conf 10 \
			-stand_emit_conf 0 \
			-nt 1 \
			-mmq 20 \
			-mbq 5 &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T BisulfiteGenotyper \
			-R $genesome/"$ref_gene".fa \
			-I $output/18_BisulfiteTableRecalibration/"$i"-BQSR.bam \
			-D $dbsnp \
			-L chr18 \
			-vfn1 $output/19_BisulfiteGenotyper/"$i"/"$i"_chr18.snp.raw.vcf \
			-out_modes EMIT_VARIANTS_ONLY \
			-stand_call_conf 10 \
			-stand_emit_conf 0 -nt 1 -mmq 20 -mbq 5 &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T BisulfiteGenotyper \
			-R $genesome/"$ref_gene".fa \
			-I $output/18_BisulfiteTableRecalibration/"$i"-BQSR.bam \
			-D $dbsnp \
			-L chr19 \
			-vfn1 $output/19_BisulfiteGenotyper/"$i"/"$i"_chr19.snp.raw.vcf \
			-out_modes EMIT_VARIANTS_ONLY \
			-stand_call_conf 10 \
			-stand_emit_conf 0 \
			-nt 1 \
			-mmq 20 \
			-mbq 5 &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T BisulfiteGenotyper \
			-R $genesome/"$ref_gene".fa \
			-I $output/18_BisulfiteTableRecalibration/"$i"-BQSR.bam \
			-D $dbsnp \
			-L chr20 \
			-vfn1 $output/19_BisulfiteGenotyper/"$i"/"$i"_chr20.snp.raw.vcf \
			-out_modes EMIT_VARIANTS_ONLY \
			-stand_call_conf 10 \
			-stand_emit_conf 0 \
			-nt 1 \
			-mmq 20 \
			-mbq 5 &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T BisulfiteGenotyper \
			-R $genesome/"$ref_gene".fa \
			-I $output/18_BisulfiteTableRecalibration/"$i"-BQSR.bam \
			-D $dbsnp \
			-L chr21 \
			-vfn1 $output/19_BisulfiteGenotyper/"$i"/"$i"_chr21.snp.raw.vcf \
			-out_modes EMIT_VARIANTS_ONLY \
			-stand_call_conf 10 \
			-stand_emit_conf 0 \
			-nt 1 \
			-mmq 20 \
			-mbq 5 &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T BisulfiteGenotyper \
			-R $genesome/"$ref_gene".fa \
			-I $output/18_BisulfiteTableRecalibration/"$i"-BQSR.bam \
			-D $dbsnp \
			-L chr22 \
			-vfn1 $output/19_BisulfiteGenotyper/"$i"/"$i"_chr22.snp.raw.vcf \
			-out_modes EMIT_VARIANTS_ONLY \
			-stand_call_conf 10 \
			-stand_emit_conf 0 \
			-nt 1 \
			-mmq 20 \
			-mbq 5 &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T BisulfiteGenotyper \
			-R $genesome/"$ref_gene".fa \
			-I $output/18_BisulfiteTableRecalibration/"$i"-BQSR.bam \
			-D $dbsnp \
			-L chrX \
			-vfn1 $output/19_BisulfiteGenotyper/"$i"/"$i"_chrX.snp.raw.vcf \
			-out_modes EMIT_VARIANTS_ONLY \
			-stand_call_conf 10 \
			-stand_emit_conf 0 \
			-nt 1 \
			-mmq 20 \
			-mbq 5 &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T BisulfiteGenotyper \
			-R $genesome/"$ref_gene".fa \
			-I $output/18_BisulfiteTableRecalibration/"$i"-BQSR.bam \
			-D $dbsnp \
			-L chrY \
			-vfn1 $output/19_BisulfiteGenotyper/"$i"/"$i"_chrY.snp.raw.vcf \
			-out_modes EMIT_VARIANTS_ONLY \
			-stand_call_conf 10 \
			-stand_emit_conf 0 \
			-nt 1 \
			-mmq 20 \
			-mbq 5 &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T BisulfiteGenotyper \
			-R $genesome/"$ref_gene".fa \
			-I $output/18_BisulfiteTableRecalibration/"$i"-BQSR.bam \
			-D $dbsnp \
			-L chrM \
			-vfn1 $output/19_BisulfiteGenotyper/"$i"/"$i"_chrM.snp.raw.vcf \
			-out_modes EMIT_VARIANTS_ONLY \
			-stand_call_conf 10 \
			-stand_emit_conf 0 \
			-nt 1 \
			-mmq 20 \
			-mbq 5 &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "2" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo '[`date`] Run complete'
		echo '------------------------------------------------------------------------------------------------------------------------'
fi

if [ ! -d "$output/20_VCFpostprocess" ]; then
		echo
		echo 
		echo "[`date`] VCFpostprocess"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/20_VCFpostprocess
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/20_VCFpostprocess/"$i"
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T VCFpostprocess \
			-R  $genesome/"$ref_gene".fa \
			-oldVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr1.snp.raw.vcf \
			-newVcf  $output/20_VCFpostprocess/"$i"/"$i"_chr1.filter.vcf \
			-snpVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr1.snp.raw.vcf \
			-o $output/20_VCFpostprocess/"$i"/"$i"_chr1.filter.summary.txt  &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T VCFpostprocess \
			-R  $genesome/"$ref_gene".fa \
			-oldVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr2.snp.raw.vcf \
			-newVcf  $output/20_VCFpostprocess/"$i"/"$i"_chr2.filter.vcf \
			-snpVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr2.snp.raw.vcf \
			-o $output/20_VCFpostprocess/"$i"/"$i"_chr2.filter.summary.txt  &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T VCFpostprocess \
			-R  $genesome/"$ref_gene".fa \
			-oldVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr3.snp.raw.vcf \
			-newVcf  $output/20_VCFpostprocess/"$i"/"$i"_chr3.filter.vcf \
			-snpVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr3.snp.raw.vcf \
			-o $output/20_VCFpostprocess/"$i"/"$i"_chr3.filter.summary.txt  &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T VCFpostprocess \
			-R  $genesome/"$ref_gene".fa \
			-oldVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr4.snp.raw.vcf \
			-newVcf  $output/20_VCFpostprocess/"$i"/"$i"_chr4.filter.vcf \
			-snpVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr4.snp.raw.vcf \
			-o $output/20_VCFpostprocess/"$i"/"$i"_chr4.filter.summary.txt  &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T VCFpostprocess \
			-R  $genesome/"$ref_gene".fa \
			-oldVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr5.snp.raw.vcf \
			-newVcf  $output/20_VCFpostprocess/"$i"/"$i"_chr5.filter.vcf \
			-snpVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr5.snp.raw.vcf \
			-o $output/20_VCFpostprocess/"$i"/"$i"_chr5.filter.summary.txt  &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T VCFpostprocess \
			-R  $genesome/"$ref_gene".fa \
			-oldVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr6.snp.raw.vcf \
			-newVcf  $output/20_VCFpostprocess/"$i"/"$i"_chr6.filter.vcf \
			-snpVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr6.snp.raw.vcf \
			-o $output/20_VCFpostprocess/"$i"/"$i"_chr6.filter.summary.txt  &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T VCFpostprocess \
			-R  $genesome/"$ref_gene".fa \
			-oldVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr7.snp.raw.vcf \
			-newVcf  $output/20_VCFpostprocess/"$i"/"$i"_chr7.filter.vcf \
			-snpVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr7.snp.raw.vcf \
			-o $output/20_VCFpostprocess/"$i"/"$i"_chr7.filter.summary.txt  &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T VCFpostprocess \
			-R  $genesome/"$ref_gene".fa \
			-oldVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr8.snp.raw.vcf \
			-newVcf  $output/20_VCFpostprocess/"$i"/"$i"_chr8.filter.vcf \
			-snpVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr8.snp.raw.vcf \
			-o $output/20_VCFpostprocess/"$i"/"$i"_chr8.filter.summary.txt  &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T VCFpostprocess \
			-R  $genesome/"$ref_gene".fa \
			-oldVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr9.snp.raw.vcf \
			-newVcf  $output/20_VCFpostprocess/"$i"/"$i"_chr9.filter.vcf \
			-snpVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr9.snp.raw.vcf \
			-o $output/20_VCFpostprocess/"$i"/"$i"_chr9.filter.summary.txt  &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T VCFpostprocess \
			-R  $genesome/"$ref_gene".fa \
			-oldVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr10.snp.raw.vcf \
			-newVcf  $output/20_VCFpostprocess/"$i"/"$i"_chr10.filter.vcf \
			-snpVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr10.snp.raw.vcf \
			-o $output/20_VCFpostprocess/"$i"/"$i"_chr10.filter.summary.txt  &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G -Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T VCFpostprocess \
			-R  $genesome/"$ref_gene".fa \
			-oldVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr11.snp.raw.vcf \
			-newVcf  $output/20_VCFpostprocess/"$i"/"$i"_chr11.filter.vcf \
			-snpVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr11.snp.raw.vcf \
			-o $output/20_VCFpostprocess/"$i"/"$i"_chr11.filter.summary.txt  &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T VCFpostprocess \
			-R  $genesome/"$ref_gene".fa \
			-oldVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr12.snp.raw.vcf \
			-newVcf  $output/20_VCFpostprocess/"$i"/"$i"_chr12.filter.vcf \
			-snpVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr12.snp.raw.vcf \
			-o $output/20_VCFpostprocess/"$i"/"$i"_chr12.filter.summary.txt  &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T VCFpostprocess \
			-R  $genesome/"$ref_gene".fa \
			-oldVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr13.snp.raw.vcf \
			-newVcf  $output/20_VCFpostprocess/"$i"/"$i"_chr13.filter.vcf \
			-snpVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr13.snp.raw.vcf \
			-o $output/20_VCFpostprocess/"$i"/"$i"_chr13.filter.summary.txt  &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T VCFpostprocess \
			-R  $genesome/"$ref_gene".fa \
			-oldVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr14.snp.raw.vcf \
			-newVcf  $output/20_VCFpostprocess/"$i"/"$i"_chr14.filter.vcf \
			-snpVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr14.snp.raw.vcf \
			-o $output/20_VCFpostprocess/"$i"/"$i"_chr14.filter.summary.txt  &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T VCFpostprocess \
			-R  $genesome/"$ref_gene".fa \
			-oldVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr15.snp.raw.vcf \
			-newVcf  $output/20_VCFpostprocess/"$i"/"$i"_chr15.filter.vcf \
			-snpVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr15.snp.raw.vcf \
			-o $output/20_VCFpostprocess/"$i"/"$i"_chr15.filter.summary.txt  &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T VCFpostprocess \
			-R  $genesome/"$ref_gene".fa \
			-oldVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr16.snp.raw.vcf \
			-newVcf  $output/20_VCFpostprocess/"$i"/"$i"_chr16.filter.vcf \
			-snpVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr16.snp.raw.vcf \
			-o $output/20_VCFpostprocess/"$i"/"$i"_chr16.filter.summary.txt  &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T VCFpostprocess \
			-R  $genesome/"$ref_gene".fa \
			-oldVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr17.snp.raw.vcf \
			-newVcf  $output/20_VCFpostprocess/"$i"/"$i"_chr17.filter.vcf \
			-snpVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr17.snp.raw.vcf \
			-o $output/20_VCFpostprocess/"$i"/"$i"_chr17.filter.summary.txt  &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T VCFpostprocess \
			-R  $genesome/"$ref_gene".fa \
			-oldVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr18.snp.raw.vcf \
			-newVcf  $output/20_VCFpostprocess/"$i"/"$i"_chr18.filter.vcf \
			-snpVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr18.snp.raw.vcf \
			-o $output/20_VCFpostprocess/"$i"/"$i"_chr18.filter.summary.txt  &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T VCFpostprocess \
			-R  $genesome/"$ref_gene".fa \
			-oldVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr19.snp.raw.vcf \
			-newVcf  $output/20_VCFpostprocess/"$i"/"$i"_chr19.filter.vcf \
			-snpVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr19.snp.raw.vcf \
			-o $output/20_VCFpostprocess/"$i"/"$i"_chr19.filter.summary.txt  &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata -jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T VCFpostprocess \
			-R  $genesome/"$ref_gene".fa \
			-oldVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr20.snp.raw.vcf \
			-newVcf  $output/20_VCFpostprocess/"$i"/"$i"_chr20.filter.vcf \
			-snpVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr20.snp.raw.vcf \
			-o $output/20_VCFpostprocess/"$i"/"$i"_chr20.filter.summary.txt  &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T VCFpostprocess \
			-R  $genesome/"$ref_gene".fa \
			-oldVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr21.snp.raw.vcf \
			-newVcf  $output/20_VCFpostprocess/"$i"/"$i"_chr21.filter.vcf \
			-snpVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr21.snp.raw.vcf \
			-o $output/20_VCFpostprocess/"$i"/"$i"_chr21.filter.summary.txt  &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T VCFpostprocess \
			-R  $genesome/"$ref_gene".fa \
			-oldVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr22.snp.raw.vcf \
			-newVcf  $output/20_VCFpostprocess/"$i"/"$i"_chr22.filter.vcf \
			-snpVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chr22.snp.raw.vcf \
			-o $output/20_VCFpostprocess/"$i"/"$i"_chr22.filter.summary.txt  &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T VCFpostprocess \
			-R  $genesome/"$ref_gene".fa \
			-oldVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chrX.snp.raw.vcf \
			-newVcf  $output/20_VCFpostprocess/"$i"/"$i"_chrX.filter.vcf \
			-snpVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chrX.snp.raw.vcf \
			-o $output/20_VCFpostprocess/"$i"/"$i"_chrX.filter.summary.txt  &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T VCFpostprocess \
			-R  $genesome/"$ref_gene".fa \
			-oldVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chrY.snp.raw.vcf \
			-newVcf  $output/20_VCFpostprocess/"$i"/"$i"_chrY.filter.vcf \
			-snpVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chrY.snp.raw.vcf \
			-o $output/20_VCFpostprocess/"$i"/"$i"_chrY.filter.summary.txt  &
			time /opt/tsinghua/cfDNApipeTest/software/zulu13.29.9-ca-jdk13.0.2-linux_x64/bin/java \
			-Xmx6G \
			-Djava.io.tmpdir=/opt/intermediata \
			-jar /opt/tsinghua/cfDNApipeTest/software/BisSNP-0.82.2.jar \
			-T VCFpostprocess \
			-R  $genesome/"$ref_gene".fa \
			-oldVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chrM.snp.raw.vcf \
			-newVcf  $output/20_VCFpostprocess/"$i"/"$i"_chrM.filter.vcf \
			-snpVcf $output/19_BisulfiteGenotyper/"$i"/"$i"_chrM.snp.raw.vcf \
			-o $output/20_VCFpostprocess/"$i"/"$i"_chrM.filter.summary.txt  &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "2" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo '[`date`] Run complete'
		echo '------------------------------------------------------------------------------------------------------------------------'
fi

if [ ! -d "$output/21_gatherVCF" ]; then
		echo
		echo 
		echo "[`date`] gatherVCF"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/21_gatherVCF
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/21_gatherVCF/"$i"
			time gatk GatherVcfs \
			-I $output/20_VCFpostprocess/"$i"/"$i"_chr1.filter.vcf \
			-I $output/20_VCFpostprocess/"$i"/"$i"_chr2.filter.vcf \
			-I $output/20_VCFpostprocess/"$i"/"$i"_chr3.filter.vcf \
			-I $output/20_VCFpostprocess/"$i"/"$i"_chr4.filter.vcf \
			-I $output/20_VCFpostprocess/"$i"/"$i"_chr5.filter.vcf \
			-I $output/20_VCFpostprocess/"$i"/"$i"_chr6.filter.vcf \
			-I $output/20_VCFpostprocess/"$i"/"$i"_chr7.filter.vcf \
			-I $output/20_VCFpostprocess/"$i"/"$i"_chr8.filter.vcf \
			-I $output/20_VCFpostprocess/"$i"/"$i"_chr9.filter.vcf \
			-I $output/20_VCFpostprocess/"$i"/"$i"_chr10.filter.vcf \
			-I $output/20_VCFpostprocess/"$i"/"$i"_chr11.filter.vcf \
			-I $output/20_VCFpostprocess/"$i"/"$i"_chr12.filter.vcf \
			-I $output/20_VCFpostprocess/"$i"/"$i"_chr13.filter.vcf \
			-I $output/20_VCFpostprocess/"$i"/"$i"_chr14.filter.vcf \
			-I $output/20_VCFpostprocess/"$i"/"$i"_chr15.filter.vcf \
			-I $output/20_VCFpostprocess/"$i"/"$i"_chr16.filter.vcf \
			-I $output/20_VCFpostprocess/"$i"/"$i"_chr17.filter.vcf \
			-I $output/20_VCFpostprocess/"$i"/"$i"_chr18.filter.vcf \
			-I $output/20_VCFpostprocess/"$i"/"$i"_chr19.filter.vcf \
			-I $output/20_VCFpostprocess/"$i"/"$i"_chr20.filter.vcf \
			-I $output/20_VCFpostprocess/"$i"/"$i"_chr21.filter.vcf \
			-I $output/20_VCFpostprocess/"$i"/"$i"_chrX.filter.vcf \
			-I $output/20_VCFpostprocess/"$i"/"$i"_chrY.filter.vcf \
			-I $output/20_VCFpostprocess/"$i"/"$i"_chrM.filter.vcf \
			-O $output/21_gatherVCF/"$i"/"$i".filter.vcf.gz &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "10" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		myvar=0
		for i in ${name[@]}
		do
			time gatk IndexFeatureFile -I $output/21_gatherVCF/"$i"/"$i".filter.vcf.gz &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "30" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo '[`date`] Run complete'
		echo '------------------------------------------------------------------------------------------------------------------------'
fi

if [ ! -d "$output/22_VariantQC" ]; then
		echo
		echo 
		echo "[`date`] VariantQC"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/22_VariantQC
		myvar=0
		for i in ${name[@]}
		do
		echo "-V:"$i" $output/21_gatherVCF/"$i"/"$i".filter.vcf.gz" >> $output/22_VariantQC/test.txt
		done
		cat $output/22_VariantQC/test.txt|tr "\n" " " > $output/22_VariantQC/ref.txt;rm -rf $output/22_VariantQC/test.txt
		ref=($(less $output/22_VariantQC/ref.txt))
		java -jar /opt/tsinghua/software/DISCVRSeq-1.3.13/DISCVRSeq-1.3.13.jar VariantQC \
		-R $genesome/"$ref_gene".fa \
		"${ref[@]}" \
		--rd $output/22_VariantQC/reports.json \
		-O $output/22_VariantQC/snp_report.html &
		wait
		echo '[`date`] Run complete'
		echo '------------------------------------------------------------------------------------------------------------------------'
fi

if [ ! -d "$output/23_annovar" ]; then
		echo
		echo 
		echo "[`date`] annovar"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/23_annovar
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/23_annovar/"$i"
			time /opt/tsinghua/cfDNApipeTest/software/annovar/annovar/table_annovar.pl \
			$output/21_gatherVCF/"$i"/"$i".filter.vcf.gz \
			/opt/tsinghua/cfDNApipeTest/software/annovar/annovar/humandb \
			--out $output/23_annovar/"$i"/"$i".snp \
			--buildver $ref_gene \
			--protocol refGene,HGMD,omim,mimTitles,morbidmap,cytoBand,wgRna,genomicSuperDups,avsnp150,cosmic70,clinvar_20210501,gwasCatalog,1000g2015aug_eas,1000g2015aug_sas,1000g2015aug_eur,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_all,esp6500siv2_all,exac03,dbnsfp35a,gnomad_genome,gnomad211_exome,gnomad211_genome,cancer_hotspots_INDEL_v2,cancer_hotspots_SNV_v2,icgc28,intervar_20180118,pmkb_simple \
			--operation g,r,r,f,f,r,r,r,f,f,f,r,f,f,f,f,f,f,f,f,f,r,r,r,f,f,r,r,f \
			--argument ',,,,,,,,,,,,,,,,,,,,,,,,,,,,' \
			--thread 30 \
			--remove \
			--nastring . \
			--vcfinput \
			--intronhgvs 300 &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "30" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo '[`date`] Run complete'
		echo '------------------------------------------------------------------------------------------------------------------------'
fi

echo "[`date`] ======================================================WGBS CNV START================================================================"

if [ ! -d "$output/24_cnvbatch" ]; then
		echo
		echo 
		echo "[`date`] cnvbatch"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/24_cnvbatch
		for i in ${name[@]}
		do
			cp $output/7_bamsort/"$i"/"$i".pair1_sorted.bam $output/7_bamsort/
		done
		wait
		bam_sort=($(ls $output/7_bamsort/*pair1_sorted.bam))
		time cnvkit.py batch \
		"${bam_sort[@]}" \
		--normal \
		--annotate $CNVkit_ref1 \
		--access $CNVkit_ref2 \
		--fasta $genesome/"$ref_gene".fa \
		-m wgs \
		-y \
		--output-reference $output/24_cnvbatch/reference.cnn \
		--output-dir $output/24_cnvbatch \
		-p 30 &
		wait
		rm -rf $output/7_bamsort/*.pair1_sorted.bam $output/7_bamsort/*.pair1_sorted.bam.bai
		echo '[`date`] Run complete'
		echo '------------------------------------------------------------------------------------------------------------------------'
fi

if [ ! -d "$output/25_cnvHeatmap" ]; then
		echo
		echo 
		echo "[`date`] cnvHeatmap"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/25_cnvHeatmap
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/25_cnvHeatmap/"$i"
			time cnvkit.py heatmap \
			$output/24_cnvbatch/"$i".pair1_sorted.cnr \
			-d \
			-y \
			-o $output/25_cnvHeatmap/"$i"/"$i"_heatmap.pdf  &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "30" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		cnr_sort=($(ls $output/24_cnvbatch/*pair1_sorted.cnr))
		time cnvkit.py heatmap "${cnr_sort[@]}" -d -y -o $output/25_cnvHeatmap/heatmap.pdf &
		echo '[`date`] Run complete'
		echo '------------------------------------------------------------------------------------------------------------------------'
fi

if [ ! -d "$output/26_cnvPlot" ]; then
		echo
		echo 
		echo "[`date`] cnvPlot"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/26_cnvPlot
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/26_cnvPlot/"$i"
			time cnvkit.py diagram \
			$output/24_cnvbatch/"$i".pair1_sorted.cnr \
			-s $output/24_cnvbatch/"$i".pair1_sorted.cns \
			-o $output/26_cnvPlot/"$i"/"$i"_diagram.pdf \
			--title "$i" \
			--threshold 0.5 \
			--min-probes 3 -y &
			time cnvkit.py scatter \
			$output/24_cnvbatch/"$i".pair1_sorted.cnr \
			-s $output/24_cnvbatch/"$i".pair1_sorted.cns \
			-o $output/26_cnvPlot/"$i"/"$i"_scatter.pdf \
			--title "$i" &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "30" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo '[`date`] Run complete'
		echo '------------------------------------------------------------------------------------------------------------------------'
fi

if [ ! -d "$output/27_cnvTable" ]; then
		echo
		echo 
		echo "[`date`] cnvTable"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/27_cnvTable
		myvar=0
		for i in ${name[@]}
		do
			mkdir -p $output/27_cnvTable/"$i"
			time cnvkit.py breaks $output/24_cnvbatch/"$i".pair1_sorted.cnr $output/24_cnvbatch/"$i".pair1_sorted.cns  --min-probes 1 \
			-o $output/27_cnvTable/"$i"/"$i"_breaks.txt &
			time cnvkit.py genemetrics $output/24_cnvbatch/"$i".pair1_sorted.cnr -s $output/24_cnvbatch/"$i".pair1_sorted.cns --threshold 0.1 --min-probes 3 -y \
			-o $output/27_cnvTable/"$i"/"$i"_genemetrics_cnrs.txt;cnvkit.py genemetrics $output/24_cnvbatch/"$i".pair1_sorted.cnr --threshold 0.1 --min-probes 3 -y \
			-o $output/27_cnvTable/"$i"/"$i"_genemetrics_cnr.txt;tail -n +2 \
			$output/27_cnvTable/"$i"/"$i"_genemetrics_cnr.txt| cut -f 1 |sort > $output/27_cnvTable/"$i"/"$i"_cnrs_gene.txt;tail -n +2 \
			$output/27_cnvTable/"$i"/"$i"_genemetrics_cnr.txt| cut -f 1 |sort > $output/27_cnvTable/"$i"/"$i"_cnr_gene.txt;comm -12 \
			$output/27_cnvTable/"$i"/"$i"_cnrs_gene.txt $output/27_cnvTable/"$i"/"$i"_cnr_gene.txt  > $output/27_cnvTable/"$i"/"$i"_genemetrics_gene.txt  &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "10" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo '[`date`] Run complete'
		echo '------------------------------------------------------------------------------------------------------------------------'
fi

echo "[`date`] ======================================================WGBS SUMMARY COUNT================================================================"

if [ ! -d "$output/summary/filter_result" ]; then
		echo
		echo 
		echo "[`date`] filter_result"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/summary/filter_result
		myvar=0
		for i in ${name[@]}
		do
			sed -n '15p' $output/2_low_complexity_filter/"$i"/"$i".json|sed 's/\t\t\t"total_reads"://'|sed 's/,//' >>$output/summary/filter_result/clean_reads.txt ;\
			sed -n '19p' $output/2_low_complexity_filter/"$i"/"$i".json|sed 's/\t\t\t"q20_rate"://'|sed 's/,//' >>$output/summary/filter_result/Q20.txt ;\
			sed -n '20p' $output/2_low_complexity_filter/"$i"/"$i".json|sed 's/\t\t\t"q30_rate"://'|sed 's/,//' >>$output/summary/filter_result/Q30.txt ;\
			sed -n '23p' $output/2_low_complexity_filter/"$i"/"$i".json|sed 's/\t\t\t"gc_content"://' >>$output/summary/filter_result/GC_content.txt ;\
			sed -n '36p' $output/2_low_complexity_filter/"$i"/"$i".json|sed 's/\t\t\"rate"://'|sed 's/[ ]//'|sed 's/,//' >>$output/summary/filter_result/duplication.txt ;\
			sed -n '41p' $output/2_low_complexity_filter/"$i"/"$i".json|sed 's/\t\t\"peak"://'|sed 's/[ ]//'|sed 's/,//' >>$output/summary/filter_result/insert_size_peak.txt 
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "20" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		paste -d '\t' $file_name $output/summary/filter_result/clean_reads.txt $output/summary/filter_result/Q20.txt $output/summary/filter_result/Q30.txt \
		$output/summary/filter_result/GC_content.txt $output/summary/filter_result/duplication.txt $output/summary/filter_result/insert_size_peak.txt|sed '1i\Sample\tclean_reads\tQ20\tQ30\tGC_content\tduplication\tinsert_size_peak' >>$output/summary/filter_result/result.txt
		echo '[`date`] Run complete'
		echo '------------------------------------------------------------------------------------------------------------------------'
fi

if [ ! -d "$output/summary/bismark_result" ]; then
		echo
		echo 
		echo "[`date`] bismark_result"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/summary/bismark_result
		myvar=0
		for i in ${name[@]}
		do
			sed -n '9p' $output/3_bismark/"$i"/"$i".pair1.truncated.gz_bismark_bt2_PE_report.txt|sed 's/Mapping efficiency:\t//' >>$output/summary/bismark_result/Mapping_effciency.txt ;\
			sed -n '26p' $output/7_bamsort/"$i"/"$i".pair1_sorted_stats/genome_results.txt|sed 's/number of mapped paired reads (both in pair) = //'|sed 's/[ ]*//' >>$output/summary/bismark_result/Paried_map.txt ;\
			sed -n '38p' $output/7_bamsort/"$i"/"$i".pair1_sorted_stats/genome_results.txt|sed 's/mean insert size = //'|sed 's/[ ]*//' >>$output/summary/bismark_result/Insert_Size.txt ;\
			sed -n '33p' $output/7_bamsort/"$i"/"$i".pair1_sorted_stats/genome_results.txt|sed 's/duplication rate = //'|sed 's/[ ]*//' >>$output/summary/bismark_result/Duplication.txt ;\
			sed -n '72p' $output/7_bamsort/"$i"/"$i".pair1_sorted_stats/genome_results.txt|sed 's/mean coverageData = //'|sed 's/[ ]*//'|sed 's/X//' >>$output/summary/bismark_result/Average_Depth.txt ;\
			sed -n '78p' $output/7_bamsort/"$i"/"$i".pair1_sorted_stats/genome_results.txt|sed 's/There is a //'|sed 's/of reference with a coverageData >= 4X//'|sed 's/[ ]*//' >>$output/summary/bismark_result/Coverage_4X.txt ;\
			sed -n '84p' $output/7_bamsort/"$i"/"$i".pair1_sorted_stats/genome_results.txt|sed 's/There is a //'|sed 's/of reference with a coverageData >= 10X//'|sed 's/[ ]*//' >>$output/summary/bismark_result/Coverage_10X.txt ;\
			sed -n '94p' $output/7_bamsort/"$i"/"$i".pair1_sorted_stats/genome_results.txt|sed 's/There is a //'|sed 's/of reference with a coverageData >= 20X//'|sed 's/[ ]*//' >>$output/summary/bismark_result/Coverage_20X.txt 
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "20" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		paste -d '\t' $file_name $output/summary/bismark_result/Mapping_effciency.txt $output/summary/bismark_result/Paried_map.txt $output/summary/bismark_result/Insert_Size.txt \
		$output/summary/bismark_result/Duplication.txt $output/summary/bismark_result/Average_Depth.txt $output/summary/bismark_result/Coverage_4X.txt \
		$output/summary/bismark_result/Coverage_10X.txt $output/summary/bismark_result/Coverage_20X.txt|sed '1i\Sample\tMapping_effciency\tParied_map\tInsert_Size\tDuplication\tAverage_Depth\tCoverage>4X\tCoverage>10X\tCoverage>20X' >>$output/summary/bismark_result/result.txt
		echo '[`date`] Run complete'
		echo '------------------------------------------------------------------------------------------------------------------------'
fi

if [ ! -d "$output/summary/cpg_result" ]; then
		echo
		echo 
		echo "[`date`] cpg_result"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/summary/cpg_result
		myvar=0
		for i in ${name[@]}
		do
			sed -n '26p' $output/3_bismark/"$i"/"$i".pair1.truncated.gz_bismark_bt2_PE_report.txt|sed "s/Total methylated C's in CpG context:\t//" >>$output/summary/cpg_result/Total_CpG.txt ;\
			sed -n '27p' $output/3_bismark/"$i"/"$i".pair1.truncated.gz_bismark_bt2_PE_report.txt|sed "s/Total methylated C's in CHG context:\t//" >>$output/summary/cpg_result/Total_CHG.txt ;\
			sed -n '28p' $output/3_bismark/"$i"/"$i".pair1.truncated.gz_bismark_bt2_PE_report.txt|sed "s/Total methylated C's in CHH context:\t//" >>$output/summary/cpg_result/Total_CHH.txt ;\
			sed -n '36p' $output/3_bismark/"$i"/"$i".pair1.truncated.gz_bismark_bt2_PE_report.txt|sed 's/C methylated in CpG context:\t//' >>$output/summary/cpg_result/mCpG.txt ;\
			sed -n '37p' $output/3_bismark/"$i"/"$i".pair1.truncated.gz_bismark_bt2_PE_report.txt|sed 's/C methylated in CHG context:\t//' >>$output/summary/cpg_result/mCHG.txt ;\
			sed -n '38p' $output/3_bismark/"$i"/"$i".pair1.truncated.gz_bismark_bt2_PE_report.txt|sed 's/C methylated in CHH context:\t//' >>$output/summary/cpg_result/mCHH.txt
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "20" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		paste -d '\t' $file_name $output/summary/cpg_result/Total_CpG.txt $output/summary/cpg_result/Total_CHG.txt $output/summary/cpg_result/Total_CHH.txt \
		$output/summary/cpg_result/mCpG.txt $output/summary/cpg_result/mCHG.txt $output/summary/cpg_result/mCHH.txt|sed '1i\Sample\tTotal_CpG\tTotal_CHG\tTotal_CHH\tmCpG\tmCHG\tmCHH' >>$output/summary/cpg_result/result.txt
		echo '------------------------------------------------------------------------------------------------------------------------'
fi

if [ ! -d "$output/15_GO_KEGG" ]; then
		echo
		echo 
		echo "[`date`] GO_KEGG"
		echo '------------------------------------------------------------------------------------------------------------------------'
		mkdir -p $output/15_GO_KEGG
		myvar=0
		source /root/anaconda3/bin/activate r_updata
		for i in ${group[@]}
		do
			mkdir -p $output/15_GO_KEGG/"$i"
			time Rscript /root/anaconda3/envs/cfDNApipe2.0/lib/python3.6/site-packages/cfDNApipe/data/Rscript/GO_kegg.r \
			--input $output/14_DMR_annotation/"$i"/genelist.txt \
			--type SYMBOL \
			--GO \
			--KEGG \
			--output_dir $output/15_GO_KEGG/"$i" \
			--pAdjustMethod BH \
			--species Hs \
			--showCategory 5 &
			myvar=$(($myvar + 1 ))
			if [ "$myvar" = "10" ]
			then
					myvar=0
					wait
			fi
		done
		wait
		echo '[`date`] Run complete'
		echo '------------------------------------------------------------------------------------------------------------------------'
fi
echo "[`date`] ======================================================WGBS FINISH================================================================"
exit 

echo "[`date`] ======================================================WGBS FINISH================================================================"
exit 
