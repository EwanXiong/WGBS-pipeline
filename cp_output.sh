output=/opt/tsinghua/NuoHe/WGBS_23.3.27/phy/output
file_name=/opt/tsinghua/NuoHe/WGBS_23.3.27/phy/name.txt
name=($(cat $file_name))
echo "==========Final_Result collect========="
if [ ! -d "$output/Final_Result" ]; then
		echo
		echo 
		echo "[`date`] Final_Result"
		echo '-------------------------'
		mkdir -p $output/Final_Result
		mkdir -p $output/Final_Result/qc
		mkdir -p $output/Final_Result/low_complexity_filter
		mkdir -p $output/Final_Result/bismark_methylation_report
		mkdir -p $output/Final_Result/Qualimap
		mkdir -p $output/Final_Result/correlation
		mkdir -p $output/Final_Result/methylLevel
		mkdir -p $output/Final_Result/methylstat

		mkdir -p $output/Final_Result/DMR
		mkdir -p $output/Final_Result/DMR_annoation
		mkdir -p $output/Final_Result/SNV
		mkdir -p $output/Final_Result/SNV/VCF
		mkdir -p $output/Final_Result/SNV/ANNOVAR
		mkdir -p $output/Final_Result/SNV/VariantQC_report

		mkdir -p $output/Final_Result/CNV
		mkdir -p $output/Final_Result/CNV/Heatmap
		mkdir -p $output/Final_Result/CNV/CNVplot
		mkdir -p $output/Final_Result/CNV/Result_txt
		
		mkdir -p $output/Final_Result/summary
		mkdir -p $output/Final_Result/GO_KEGG
		wait
		cp -r $output/1_qc/* $output/Final_Result/qc &
		cp -r $output/6_bismark_methylation_report/* $output/Final_Result/bismark_methylation_report &
		cp $output/10_correlation/clusterSamples.pdf $output/Final_Result/correlation &
		cp $output/10_correlation/Correlation.pdf $output/Final_Result/correlation &
		cp $output/10_correlation/PCA.pdf $output/Final_Result/correlation &
		cp -r $output/11_methylLevel/* $output/Final_Result/methylLevel &
		cp -r $output/12_methylstat/* $output/Final_Result/methylstat &
		
		cp -r $output/13_DMR/* $output/Final_Result/DMR &
		cp -r $output/14_DMR_annotation/* $output/Final_Result/DMR_annoation &
		cp -r $output/21_gatherVCF/* $output/Final_Result/SNV/VCF &
		cp -r $output/22_VariantQC/* $output/Final_Result/SNV/VariantQC_report &
		cp -r $output/23_annovar/* $output/Final_Result/SNV/ANNOVAR &
		
		cp -r $output/25_cnvHeatmap/* $output/Final_Result/CNV/Heatmap &
		cp -r $output/26_cnvPlot/* $output/Final_Result/CNV/CNVplot &
		cp -r $output/summary/* $output/Final_Result/summary &
		cp -r $output/15_GO_KEGG/* $output/Final_Result/GO_KEGG &

		
		for i in ${name[@]}
		do
			cp $output/2_low_complexity_filter/"$i"/"$i".html $output/Final_Result/low_complexity_filter &
			cp $output/2_low_complexity_filter/"$i"/"$i".json $output/Final_Result/low_complexity_filter &
			cp $output/27_cnvTable/"$i"/"$i"_genemetrics_cnrs.txt $output/Final_Result/CNV/Result_txt &
			cp -r $output/7_bamsort/"$i"/"$i".pair1_sorted_stats $output/Final_Result/Qualimap &
		done
		wait
		echo '[`date`] Run complete'
		echo '------------------------'
fi
