version 1.0

struct FastaSamples {
    File fastq_r1
    File fastq_r2
    String barcode    
}

workflow Eclip {
    
    input {
        Array[FastaSamples] samples
        File zipped_star_files
        File zipped_star_files_to_hg19
    }
    
    scatter (sample in samples) {
    call CutAdapt {
        input:
        fastq_r1 = sample.fastq_r1,
        fastq_r2 = sample.fastq_r2,
        barcode = sample.barcode
    }
    call FastQC_round1 {
        input:
        fastqc_r1 = CutAdapt.result_round1_cutadapt_left,
        fastqc_r2 = CutAdapt.result_round1_cutadapt_right
    }
    call CutAdapt_round2 {
        input:
        round1_left_r1 = CutAdapt.result_round1_cutadapt_left,
        round1_right_r2 = CutAdapt.result_round1_cutadapt_right,
        barcode = sample.barcode
    }
    call FastQC_round2 {
        input:
        fastqc_round2_r1 = CutAdapt_round2.result_round2_cutadapt_left,
        fastqc_round2_r2 = CutAdapt_round2.result_round2_cutadapt_right
    }
    call FastQ_sort {
        input:
        fastq_sort_r1 = CutAdapt_round2.result_round2_cutadapt_left,
        fastq_sort_r2 = CutAdapt_round2.result_round2_cutadapt_right
    }
    call STAR_rmRep {
        input:
        fastq_starrep_r1 = FastQ_sort.result_fastq_sort_left,
        fastq_starrep_r2 = FastQ_sort.result_fastq_sort_right,
        zipped_star_files = zipped_star_files
    }
    call FastQ_sort_STAR_unmapped {
        input:
        unmapped_to_sort_r1 = STAR_rmRep.result_star_fq_r1,
        unmapped_to_sort_r2 = STAR_rmRep.result_star_fq_r2
    }
    call STAR_genome_map {
        input:
        sorted_star_fq_r1 = FastQ_sort_STAR_unmapped.result_fastq_sort_after_rmRep_r1,
        sorted_star_fq_r2 = FastQ_sort_STAR_unmapped.result_fastq_sort_after_rmRep_r2,
        zipped_star_files_to_hg19 = zipped_star_files_to_hg19
    }
    }
}

task CutAdapt {
    
    input {
        File fastq_r1
        File fastq_r2
        String barcode
    }
    
    String left_r1 = basename(fastq_r1,'.gz')
    String right_r2 = basename(fastq_r2,'.gz')

    command <<<
    source /groups/cgsd/alexandre/miniconda3/etc/profile.d/conda.sh 
    conda activate stepbystep
    cutadapt --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 6,6 -m 18 \
    -a ~{barcode} \
    -g ~{barcode} \
    -A ~{barcode} \
    -G ~{barcode} \
    -o ~{left_r1} \
    -p ~{right_r2} \
    ~{fastq_r1} \
    ~{fastq_r2}
    >>>

    runtime {
        cpu: 3
        memory: "6 GB"
    }

    output {
        File result_round1_cutadapt_left = "${left_r1}"
        File result_round1_cutadapt_right = "${right_r2}"
    }

}

task FastQC_round1 {
    input {
        File fastqc_r1
        File fastqc_r2
    }
    command <<<
    source /groups/cgsd/alexandre/miniconda3/etc/profile.d/conda.sh 
    conda activate stepbystep
    fastqc -t 2 --extract -k 7 ~{fastqc_r1} -o .
    fastqc -t 2 --extract -k 7 ~{fastqc_r2} -o .
    >>>
    runtime {
        cpu: 3
        memory: "5 GB"
    }
}

task CutAdapt_round2 {
    input {
        File round1_left_r1
        File round1_right_r2
        String barcode
    }

    String round2_left_r1 = basename(round1_left_r1,'fq') + 'round2.fq'
    String round2_right_r2 = basename(round1_right_r2,'fq') + 'round2.fq'

    command <<<
    source /groups/cgsd/alexandre/miniconda3/etc/profile.d/conda.sh 
    conda activate stepbystep
    cutadapt --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 6,6 -m 18 \
    -a ~{barcode} \
    -g ~{barcode} \
    -A ~{barcode} \
    -G ~{barcode} \
    -o ~{round2_left_r1} \
    -p ~{round2_right_r2} \
    ~{round1_left_r1} \
    ~{round1_right_r2}
    >>>

    runtime {
        cpu: 3
        memory: "6 GB"
    }

    output {
        File result_round2_cutadapt_left = "${round2_left_r1}" 
        File result_round2_cutadapt_right = "${round2_right_r2}"
    }
}

task FastQC_round2 {
    input {
        File fastqc_round2_r1
        File fastqc_round2_r2
    }
    command <<<
    source /groups/cgsd/alexandre/miniconda3/etc/profile.d/conda.sh 
    conda activate stepbystep
    fastqc -t 2 --extract -k 7 ~{fastqc_round2_r1} -o .
    fastqc -t 2 --extract -k 7 ~{fastqc_round2_r2} -o .
    >>>
    
    runtime {
        cpu: 3
        memory: "5 GB"
    }

}

task FastQ_sort {
    input {
        File fastq_sort_r1
        File fastq_sort_r2
    }
    String sorted_r1 = basename(fastq_sort_r1,'.fq') + '.sorted.fq'
    String sorted_r2 = basename(fastq_sort_r2,'.fq') + '.sorted.fq' 

    command <<<
    source /groups/cgsd/alexandre/miniconda3/etc/profile.d/conda.sh 
    conda activate stepbystep
    fastq-sort --id ~{fastq_sort_r1} > ~{sorted_r1}
    fastq-sort --id ~{fastq_sort_r2} > ~{sorted_r2}
    >>>
    runtime {
        cpu: 3
        memory: "5 GB"
    }
    output {
        File result_fastq_sort_left = "${sorted_r1}"
        File result_fastq_sort_right = "${sorted_r2}"
     }
}


task STAR_rmRep {
    input {
        File zipped_star_files
        File fastq_starrep_r1
        File fastq_starrep_r2
    }

    String prefix = sub(basename(fastq_starrep_r1,'.fq'),'_r1','') + "_STAR_"

    command <<<
    mkdir RepElements
    tar -xzf ~{zipped_star_files} -C RepElements
    source /groups/cgsd/alexandre/miniconda3/etc/profile.d/conda.sh 
    conda activate stepbystep
    STAR \
    --runMode alignReads \
    --runThreadN 8 \
    --genomeDir RepElements \
    --genomeLoad NoSharedMemory \
    --alignEndsType EndToEnd \
    --outSAMunmapped Within \
    --outFilterMultimapNmax 30 \
    --outFilterMultimapScoreRange 1 \
    --outFileNamePrefix ~{prefix} \
    --outSAMtype BAM Unsorted \
    --outFilterType BySJout \
    --outBAMcompression 10 \
    --outReadsUnmapped Fastx \
    --outFilterScoreMin 10 \
    --outSAMattrRGline ID:foo \
    --outSAMattributes All \
    --outSAMmode Full \
    --outStd Log \
    --readFilesIn ~{fastq_starrep_r1} ~{fastq_starrep_r2}
    >>>
    runtime {
        cpu: 4
        memory: "20 GB"
    }
    output {
        File result_star_fq_r1 = "${prefix}Unmapped.out.mate1"
        File result_star_fq_r2 = "${prefix}Unmapped.out.mate2"
        File result_star_bam = "${prefix}Aligned.out.bam"
    }
}

task FastQ_sort_STAR_unmapped {
    input {
        File unmapped_to_sort_r1
        File unmapped_to_sort_r2
    }
    String sorted_r1 = basename(unmapped_to_sort_r1,'Unmapped.out.mate1') + 'r1_.fq'
    String sorted_r2 = basename(unmapped_to_sort_r2,'Unmapped.out.mate2') + 'r2_.fq' 

    command <<<
    source /groups/cgsd/alexandre/miniconda3/etc/profile.d/conda.sh 
    conda activate stepbystep
    fastq-sort --id ~{unmapped_to_sort_r1} > ~{sorted_r1}
    fastq-sort --id ~{unmapped_to_sort_r2} > ~{sorted_r2}
    >>>
    runtime {
        cpu: 3
        memory: "5 GB"
    }
    output {
        File result_fastq_sort_after_rmRep_r1 = "${sorted_r1}"
        File result_fastq_sort_after_rmRep_r2 = "${sorted_r2}"
     }
}

task STAR_genome_map {
    input {
        File sorted_star_fq_r1
        File sorted_star_fq_r2
        File zipped_star_files_to_hg19
    }
    String prefix = basename(sorted_star_fq_r1,'r1_.fq') + 'hg19'

    command <<<
    mkdir HG_19_DIR
    tar -xzf ~{zipped_star_files_to_hg19} -C HG_19_DIR
    source /groups/cgsd/alexandre/miniconda3/etc/profile.d/conda.sh 
    conda activate stepbystep
    STAR \
    --runMode alignReads \
    --runThreadN 8 \
    --genomeDir  HG_19_DIR \
    --genomeLoad NoSharedMemory \
    --readFilesIn ~{sorted_star_fq_r1} ~{sorted_star_fq_r2} \
    --outSAMunmapped Within \
    --outFilterMultimapNmax 1 \
    --outFilterMultimapScoreRange 1 \
    --outFileNamePrefix ~{prefix} \
    --outSAMattributes All \
    --outSAMtype BAM Unsorted \
    --outFilterType BySJout \
    --outReadsUnmapped Fastx \
    --outFilterScoreMin 10 \
    --outSAMattrRGline ID:foo \
    --outStd Log \
    --alignEndsType EndToEnd \
    --outBAMcompression 10 \
    --outSAMmode Full
    >>>
    runtime {
        cpu: 4
        memory: "30 GB"
    }
    output {
        File result_star_hg19_fq_r1 = "${prefix}Unmapped.out.mate1"
        File result_star_hg19_fq_r2 = "${prefix}Unmapped.out.mate2"
        File result_star_hg19_bam = "${prefix}Aligned.out.bam"
    }
}