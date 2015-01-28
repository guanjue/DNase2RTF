# if statement () vs []


########################## 
############
######
# bash format_scripts.sh allpks_bed_names.txt all_pks.bed
######
allpks2onepk_bed(){
    allinputfilename=$1
    alloutputfilename=$2
    rm $alloutputfilename
    cat $allinputfilename | while read pksfiles
    do
        #put the pks of target cell and source cell into one bed file
        cat $pksfiles >> $alloutputfilename
    done
}
############
##########################
#allpks2onepk_bed allpks_bed_names.txt all_pks.bed 

##########################
############
######
# bash format_scripts.sh bed2mid all_pks_sorted_merged_noTxexpan_normsk.bed
######
bed2mid(){
    bedfile=$1
    cat $bedfile | awk -F '\t' -v OFS='\t' '($2+$3)%2==1 {print $1":"($2+$3+1)/2}' > $bedfile'_remaider1.mid'
    cat $bedfile | awk -F '\t' -v OFS='\t' '($2+$3)%2==0 {print $1":"($2+$3)/2}' > $bedfile'_remaider0.mid'
    ### remove previous midpoints file
    rm $bedfile'.mid'
    ### put all midpoints into one file
    cat $bedfile'_remaider1.mid' >> $bedfile'.mid'
    cat $bedfile'_remaider0.mid' >> $bedfile'.mid'
    rm $bedfile'_remaider1.mid'
    rm $bedfile'_remaider0.mid'
}
############
##########################


##########################
############
######
# bash format_scripts.sh get_top_pks T_reads_signal S_reads_signal 3000 1500 0
######
get_top_pks(){
    target_sig=$1
    source_sig=$2
    top_num_sig=$3 ### number of pks regard as the top signal pks
    top_num_fold=$4 ### number of pks regard as the top fold change pks
    remove_interm=$5
    # put target signal and source signal into one file
    paste $target_sig $source_sig > TS_sig_OD
	tail -n +2 TS_sig_OD > TS_sig_OD_noheader
	cat TS_sig_OD_noheader | awk -F '\t' -v OFS='\t' '{print $1,$2+1,$4+1}' > TS_sig_noheader
	sort TS_sig_noheader | uniq > TS_sig_noheader_uniq
	echo -e "pks\tT_sig\tS_sig" | cat - TS_sig_noheader_uniq > temp && mv temp TS_sig_EdgeR
    cat TS_sig_OD_noheader | awk -F '\t' -v OFS='\t' '{print $1,$2,$4,($2+1)/($4+1)}' > TS_sig.txt
    # pull out the pks with highest signal in target cell
    sort -k2,1nr TS_sig.txt > T_sig_sorted.txt
    head -$top_num_sig T_sig_sorted.txt > T_sig_sorted_top.txt
    # pull out the pks with highest signal in source cell
    sort -k3,1nr TS_sig.txt > S_sig_sorted.txt
    head -$top_num_sig S_sig_sorted.txt > S_sig_sorted_top.txt
    # put pks with top target signal and top source signal into one file
    rm TS_sig_sorted_top.txt
    cat T_sig_sorted_top.txt >> TS_sig_sorted_top.txt
    cat S_sig_sorted_top.txt >> TS_sig_sorted_top.txt
    # remove the repeats
    sort TS_sig_sorted_top.txt | uniq > TS_sig_sorted_top_uniq.txt
    # get the top fold pks of target & source cell
    sort -k4,1nr TS_sig_sorted_top_uniq.txt |head -$top_num_fold | awk -F '\t' -v OFS='\t' '{print $1}'> target_high.pks
    sort -k4,1nr TS_sig_sorted_top_uniq.txt |tail -$top_num_fold | awk -F '\t' -v OFS='\t' '{print $1}'> source_high.pks
    # remove intermediate data
	echo 'remove intermediate data'
	rm TS_sig_OD_noheader
	rm TS_sig_OD_noheader_uniq
	rm TS_sig_noheader
    rm TS_sig.txt
    rm T_sig_sorted.txt
    rm S_sig_sorted.txt
    rm T_sig_sorted_top.txt
    rm S_sig_sorted_top.txt
    rm TS_sig_sorted_top.txt
}
############
##########################


##########################
############
######
# bash format_scripts.sh edgeR2bed target_high.pks ~/group/genomes/mm9/mm9.fa
######
edgeR2bed(){
    edgeRfile=$1
    genome_fasta=$2
    cat $edgeRfile | awk -F ':' -v OFS='\t' '{print $1,$2}' > $edgeRfile'_chrsep'
    cat $edgeRfile'_chrsep' | awk -F '-' -v OFS='\t' '{print $1,$2}' > $edgeRfile'.bed'
    bedtools getfasta -fi $genome_fasta -bed $edgeRfile'.bed' -fo $edgeRfile'.bed.fasta'
    rm $edgeRfile'_chrsep'
}
############
##########################


##########################
############
# bash format_scripts.sh bed2readscount "Mus musculus;mm9" "Zhou Islet ATACseq C57BL6-129;bowtie_unique" all_pks_sorted_merged_noTxexpan_normsk.bed.mid 150 T_reads_signal
######
bed2readscount(){
    species_genome=$1
    experiments=$2
    midfile=$3
    window=$4
    outputname=$5
    if ( "$experiments" == 'Zhou Islet ATACseq C57BL6-129;bowtie_unique' )
        then
            java -cp ~/group/code/seqcode.mahonylab.jar  edu.psu.compbio.seqcode.projects.shaun.PeaksAnalysis --species "$species_genome" --peaks $midfile --out $outputname --rdbexpt "$experiments" --win $window --counts --fixedpb 1
    else
            java -cp ~/group/code/seqcode.mahonylab.jar  edu.psu.compbio.seqcode.projects.shaun.PeaksAnalysis --species "$species_genome" --peaks $midfile --out $outputname --rdbexpt "$experiments" --win $window --counts
    fi

}
############
##########################




