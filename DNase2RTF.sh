#!/bin/sh
#  DNase2RTF.sh
#  
#
#  Created by Guanjue Xiang on 1/26/15.

DNase2RTF(){
    arg_folder=$1
    Tcell_type=$2
    Scell_type=$3
    refseq_expan=$4
    rmsk=$5
    species_genome=$6
    Targetcell_experiment=$7
    Sourcecell_experiment=$8
    window=$9


    cd $arg_folder

    ### Merge the pks of target cell and source cell
    ##########################
    echo 'Merge the pks of target cell and source cell'
    if ( bash format_scripts_bash.sh allpks2onepk_bed allpks_bed_names.txt all_pks.bed )
        then
            echo 'successfully putting all bedfiles into one bed file'
    else
            echo ERROR: 'Please recheck bed files'
            exit 1

    echo 'Sorting pks'
    if ( sortBed -i all_pks.bed > all_pks_sorted.bed )
        then
            echo 'successfully sorted all_pks_sorted.bed'
    else
            echo ERROR: 'sorting all_pks.bed, Please recheck pks files'
            exit 1
    fi

    echo 'Merging pks'
    if ( bedtools merge -i all_pks_sorted.bed > all_pks_sorted_merged.bed )
        then
            echo 'successfully merged all_pks_sorted_merged.bed'
    else
            echo ERROR: 'merging all_pks_sorted.bed'
            exit 1
    fi
    echo 'Merge the pks Successfully DONE'
    ##########################


    ### Subtract Ref Seq Tx expand regions & subtract Repeat Maskers
    ##########################
    echo 'Subtracting refseq expand transcription regions'
    if ( bedtools subtract -a all_pks_sorted_merged.bed -b $refseq_expan > all_pks_sorted_merged_noTxexpan.bed )
        then 
            echo 'successfully subtracted RefSeq Tx expanded regions'
    else
            echo ERROR: 'subtracting Ref Seq Tx expanded regions'
            exit 1
    fi
 
    echo 'Subtracting Repeat Maskers regions'
    if ( bedtools subtract -a all_pks_sorted_merged_noTxexpan.bed -b $rmsk > all_pks_sorted_merged_noTxexpan_normsk.bed )
        then 
            echo 'successfully subtracted Repeat Maskers regions'
    else
            echo ERROR: 'subtracting Repeat Maskers regions'
            exit 1
    fi
    echo 'Subtracting refseq expand transcription regions Successfully DONE'
    ##########################


    ### Get midpoints of the all_pks_sorted_merged_noTxexpan_normsk
    echo 'Getting midpoints'
    if ( bash format_scripts_bash.sh bed2mid all_pks_sorted_merged_noTxexpan_normsk.bed )
        then 
            echo 'Successfully Got midpoints'
    else
            echo ERROR: 'Getting midpoints'
            exit 1

    echo 'Getting midpoints Successfully DONE' 
    ##########################


    ### pull out reads counts of the pks
    if ( bash format_scripts.sh bed2readscount $species_genome $Targetcell_experiment all_pks_sorted_merged_noTxexpan_normsk.bed.mid $window T_reads_signal)
        then
            echo 'Successfully pull out reads counts of target cell experiments'
    else
            echo ERROR: 'in calculating target cell reads counts'
            exit 1
    fi

    if ( bash format_scripts.sh bed2readscount $species_genome $Sourcecell_experiment all_pks_sorted_merged_noTxexpan_normsk.bed.mid $window S_reads_signal)
        then
            echo 'Successfully pull out reads counts of source cell experiments'
    else
            echo ERROR: 'in calculating source cell reads counts'
            exit 1
    fi
    ##########################


    echo 'extract target cell pks with highest target cell signals and source cell pks with highest source cell signals, then further extract the pks with highest fold changes'
    if ( bash format_scripts_bash.sh get_top_pks T_reads_signal S_reads_signal 3000 1500 0 )
        then
            echo 'successfully pulled out the top pks with top fold change'
    else
            echo ERROR: 'Please recheck the reads count files: T_reads_signal & S_reads_signal'
            exit 1
    fi

    if ( bash format_scripts_bash.sh edgeR2bed target_high.pks $genome_fasta) && ( bash format_scripts_bash.sh edgeR2bed source_high.pks $genome_fasta) 
        then 
            echo 'successfully pulled out the top pks with top fold change in a bed format (target/source_high.pks.bed) and fasta format (target/source_high.pks.bed.fasta)'
    else
            echo ERROR: 'Please recheck the reads count files: *.pks files'
            exit 1
    fi
    ##########################


    ### Use DREME to detect motifs enriched in target cell (source cell as background) 
    ##########################
    echo 'Use DREME to detect motifs enriched in target cell (source cell as background)'
    if ( dreme -mink 6 -maxk 8 -oc $Tcell_type'_dreme_result'  -p target_high.pks.bed.fasta -n source_high.pks.bed.fasta -e 0.001 )
        then 
            echo 'successfully run DREME to detect motifs enriched in target cell (source cell as background)'
    else
            echo ERROR: 'DREME problem'
            exit 1
    fi
    
    cd $Tcell_type'_dreme_result'
    if ( python format_scripts_python.py dreme2meme('dreme.txt') )
        then 
            echo 'successfully converted DREME result to MEME format'
    else
            echo ERROR: 'Please recheck DREME results'
            exit 1
    fi
    ##########################


    ### Use STAMP to match the DREME found motifs with familial motifs
    ##########################
    if ( perl ~/group/software/stamp/FormatMatrices.pl dreme.txt.meme dreme.txt.meme.transfac )
        then
            echo 'successfully converted MEME result to TRANSFAC format'
    else
            echo ERROR: 'Please recheck dreme.txt & dreme.txt.meme'
            exit 1
    fi
    ### Use STAMP to match the DREME motifs to the CIS familial motifs
    if ( ~/group/software/stamp/STAMP.v1.2/code/stamp -tf dreme.txt.meme.transfac -align SWU -sd ~/group/software/stamp/ScoreDists/JaspRand_PCC_SWU.scores -out target_specific_high_dreme_motif -chp -match_top 3 -match ~/group/data/motif_datasets/CIS_familial_motifs/human_mouse_CIS_family_merged_group_motif_spe_outliers.meme.meme.transfac )
        then
            echo 'successfully use STAMP to find the DREME motifs matched familial motifs' 
    else
            echo ERROR: 'Please recheck dreme.txt.meme.transfac'
            exit 1
    fi


}
DNase2RTF ~/group/projects/guanjue/PIPELINE_results/islet_acinar islet Acinar ~/group/projects/guanjue/genome_data/refGene_outer-coords_10000expan.mm9.txt ~/group/projects/guanjue/genome_data/mm9rmsk "Mus musculus;mm9" "Zhou Islet ATACseq C57BL6-129;bowtie_unique" "Zhou Acinar DnaseSeq(50-100) C57BL6-129;bowtie2_unique" 150
