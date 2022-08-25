version 1.0

struct Samples {
    File fastq_r2
    String barcode    
}

workflow Eprint {
    
    input {
        Samples samples
        File hg19_dup_tar
        File hg19_tar
    }
    
    call UmiTools {
        input:
        fastq_r2 = samples.fastq_r2,
        barcode = samples.barcode
    }
    
    call CutAdapt {
        input:
        umi_r2 = UmiTools.umi_r2,
        barcode = samples.barcode
    }

    # call FastQC {
    #     input:
    #     cut_r2 = CutAdapt.umi_r2 

    # }
}

task UmiTools {
    input {
        File fastq_r2
        String barcode  
    }

    String r2 = basename(fastq_r2,'.fq.gz')

    command <<<
    eval "$(conda shell.bash hook)"
    conda activate eprint
    umi_tools extract \
    --random-seed 1 \
    --bc-pattern NNNNNNNNNN \
    --log ~{r2}.metrics \
    --stdin ~{fastq_r2} \
    --stdout ~{r2}.fq.gz
    >>>

    runtime  {
        cpu: 6
        memory: "16 GB"
    }

    output {
        File umi_r2 = "~{r2}.fq.gz"
    }
}

task CutAdapt {
    
    input {
        File umi_r2
        String barcode
    }
    
    String r2 = basename(umi_r2,'.fq.gz')

    command <<<
    eval "$(conda shell.bash hook)" 
    conda activate eprint
    cutadapt \
    --cores=16 \
    -a ~{barcode} \
    -g ~{barcode} \
    -A ~{barcode} \
    -G ~{barcode} \
    -o ~{r2}.fq.gz \
    ~{umi_r2} 
    >>>

    runtime {
        cpu: 16
        memory: "6 GB" 
    }

    output {
        File cut_r2 = "~{r2}.fq.gz"
    }
}

# task FastQC {
#     input {
#         File fastqc_r1
#         File fastqc_r2
#     }
#     command <<<
#     eval "$(conda shell.bash hook)" 
#     conda activate eprint
#     fastqc -t 2 --extract -k 7 ~{fastqc_r1} -o .
#     fastqc -t 2 --extract -k 7 ~{fastqc_r2} -o .
#     >>>
#     runtime {
#         cpu: 3
#         memory: "5 GB"
#     }
# }

# task CutAdapt_round2 {
#     input {
#         File round1_left_r1
#         File round1_right_r2
#         String barcode
#     }

#     String round2_left_r1 = basename(round1_left_r1,'fq') + 'round2.fq'
#     String round2_right_r2 = basename(round1_right_r2,'fq') + 'round2.fq'

#     command <<<
#     eval "$(conda shell.bash hook)" 
#     conda activate eprint
#     cutadapt --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 6,6 -m 18 \
#     -a ~{barcode} \
#     -g ~{barcode} \
#     -A ~{barcode} \
#     -G ~{barcode} \
#     -o ~{round2_left_r1} \
#     -p ~{round2_right_r2} \
#     ~{round1_left_r1} \
#     ~{round1_right_r2}
#     >>>

#     runtime {
#         cpu: 3
#         memory: "6 GB"
#     }

#     output {
#         File result_round2_cutadapt_left = "${round2_left_r1}" 
#         File result_round2_cutadapt_right = "${round2_right_r2}"
#     }
# }

# task FastQC_round2 {
#     input {
#         File fastqc_round2_r1
#         File fastqc_round2_r2
        
#     }
#     command <<<
#     eval "$(conda shell.bash hook)" 
#     conda activate eprint
#     fastqc -t 2 --extract -k 7 ~{fastqc_round2_r1} -o .
#     fastqc -t 2 --extract -k 7 ~{fastqc_round2_r2} -o .
#     >>>
    
#     runtime {
#         cpu: 3
#         memory: "5 GB" 
#     }

# }

# task FastQ_sort {
#     input {
#         File fastq_sort_r1
#         File fastq_sort_r2
        
#     }
#     String sorted_r1 = basename(fastq_sort_r1,'.fq') + '.sorted.fq'
#     String sorted_r2 = basename(fastq_sort_r2,'.fq') + '.sorted.fq' 

#     command <<<
#     eval "$(conda shell.bash hook)" 
#     conda activate eprint
#     fastq-sort --id ~{fastq_sort_r1} > ~{sorted_r1}
#     fastq-sort --id ~{fastq_sort_r2} > ~{sorted_r2}
#     >>>
#     runtime {
#         cpu: 3
#         memory: "5 GB"
#     }
#     output {
#         File result_fastq_sort_left = "${sorted_r1}"
#         File result_fastq_sort_right = "${sorted_r2}"
#      }
# }


# task STAR_rmRep {
#     input {
#         File hg19_dup_tar
#         File fastq_starrep_r1
#         File fastq_starrep_r2  
#     }

#     String prefix = sub(basename(fastq_starrep_r1,'.fq'),'_r1','') + "_STAR_"

#     command <<<
#     mkdir RepElements
#     tar -xzf ~{hg19_dup_tar} -C RepElements
#     eval "$(conda shell.bash hook)" 
#     conda activate eprint
#     STAR \
#     --runMode alignReads \
#     --runThreadN 6 \
#     --genomeDir RepElements \
#     --genomeLoad NoSharedMemory \
#     --alignEndsType EndToEnd \
#     --outSAMunmapped Within \
#     --outFilterMultimapNmax 30 \
#     --outFilterMultimapScoreRange 1 \
#     --outFileNamePrefix ~{prefix} \
#     --outSAMtype BAM Unsorted \
#     --outFilterType BySJout \
#     --outBAMcompression 10 \
#     --outReadsUnmapped Fastx \
#     --outFilterScoreMin 10 \
#     --outSAMattrRGline ID:foo \
#     --outSAMattributes All \
#     --outSAMmode Full \
#     --outStd Log \
#     --readFilesIn ~{fastq_starrep_r1} ~{fastq_starrep_r2}
#     >>>
#     runtime {
#         cpu: 4
#         memory: "20 GB"
#     }
#     output {
#         File result_star_fq_r1 = "${prefix}Unmapped.out.mate1"
#         File result_star_fq_r2 = "${prefix}Unmapped.out.mate2"
#         File result_star_bam = "${prefix}Aligned.out.bam"
#     }
# }

# task FastQ_sort_STAR_unmapped {
#     input {
#         File unmapped_to_sort_r1
#         File unmapped_to_sort_r2
        
#     }
#     String sorted_r1 = basename(unmapped_to_sort_r1,'Unmapped.out.mate1') + 'r1_.fq'
#     String sorted_r2 = basename(unmapped_to_sort_r2,'Unmapped.out.mate2') + 'r2_.fq' 

#     command <<<
#     eval "$(conda shell.bash hook)" 
#     conda activate eprint
#     fastq-sort --id ~{unmapped_to_sort_r1} > ~{sorted_r1}
#     fastq-sort --id ~{unmapped_to_sort_r2} > ~{sorted_r2}
#     >>>
#     runtime {
#         cpu: 3
#         memory: "5 GB"
        
#     }
#     output {
#         File result_fastq_sort_after_rmRep_r1 = "${sorted_r1}"
#         File result_fastq_sort_after_rmRep_r2 = "${sorted_r2}"
#      }
# }

# task STAR_genome_map {
#     input {
#         File sorted_star_fq_r1
#         File sorted_star_fq_r2
#         File hg19_tar
#     }
#     String prefix = basename(sorted_star_fq_r1,'r1_.fq') + 'hg19'

#     command <<<
#     mkdir HG_19_DIR
#     tar -xzf ~{hg19_tar} -C HG_19_DIR
#     eval "$(conda shell.bash hook)" 
#     conda activate eprint
#     STAR \
#     --runMode alignReads \
#     --runThreadN 6 \
#     --genomeDir  HG_19_DIR \
#     --genomeLoad NoSharedMemory \
#     --readFilesIn ~{sorted_star_fq_r1} ~{sorted_star_fq_r2} \
#     --outSAMunmapped Within \
#     --outFilterMultimapNmax 1 \
#     --outFilterMultimapScoreRange 1 \
#     --outFileNamePrefix ~{prefix} \
#     --outSAMattributes All \
#     --outSAMtype BAM Unsorted \
#     --outFilterType BySJout \
#     --outReadsUnmapped Fastx \
#     --outFilterScoreMin 10 \
#     --outSAMattrRGline ID:foo \
#     --outStd Log \
#     --alignEndsType EndToEnd \
#     --outBAMcompression 10 \
#     --outSAMmode Full
#     >>>
#     runtime {
#         cpu: 4
#         memory: "30 GB"
#     }
#     output {
#         File result_star_hg19_fq_r1 = "${prefix}Unmapped.out.mate1"
#         File result_star_hg19_fq_r2 = "${prefix}Unmapped.out.mate2"
#         File result_star_hg19_bam = "${prefix}Aligned.out.bam"
#     }
# }