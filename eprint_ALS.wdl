version 1.0

struct Samples {
    File fastq_r2
    Array[String] adapters
    String bc_pattern    
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
        bc_pattern = samples.bc_pattern
    }
    
    call CutAdapt {
        input:
        umi_r2 = UmiTools.umi_r2,
        adapters = samples.adapters
    }

    call FastQC {
        input:
        cut_r2 = CutAdapt.cut_r2
    }

    call FastQ_sort {
        input:
        cut_r2 = CutAdapt.cut_r2
    }

    call STAR_rmRep {
        input:
        sorted_cut_r2 = FastQ_sort.fastq_sort_r2,
        hg19_dup_tar = hg19_dup_tar
    }

    call FastQ_sort_STAR {
        input:
        star_r2 = STAR_rmRep.star_r2,
        star_bam = STAR_rmRep.star_bam
    }

    call STAR_genome_map {
        input:
        sorted_start_r2 = FastQ_sort_STAR.sorted_start_r2,
        star_bam = STAR_rmRep.star_bam,
        hg19_tar =hg19_tar
    }

}

task UmiTools {
    input {
        File fastq_r2
        String bc_pattern  
    }

    String r2 = basename(fastq_r2,'.fastq.gz')

    command <<<
    eval "$(conda shell.bash hook)"
    conda activate eprint
    umi_tools extract \
    --random-seed 1 \
    --bc-pattern ~{bc_pattern} \
    --log ~{r2}.metrics \
    --stdin ~{fastq_r2} \
    --stdout ~{r2}.fastq.gz
    >>>

    runtime  {
        cpu: 6
        memory: "16 GB"
    }

    output {
        File umi_r2 = "~{r2}.fastq.gz"
    }
}

task CutAdapt {
    
    input {
        File umi_r2
        Array[String] adapters
    }
    
    String r2 = basename(umi_r2,'.fastq.gz')

    command <<<
    eval "$(conda shell.bash hook)" 
    conda activate eprint
    cutadapt \
    --cores=16 \
    -g ~{sep="\\\n -g " adapters } \
    -o ~{r2}.fastq.gz \
    ~{umi_r2} 
    >>>

    runtime {
        cpu: 16
        memory: "6 GB" 
    }

    output {
        File cut_r2 = "~{r2}.fastq.gz"
    }
}

task FastQC {
    input {
        File cut_r2
    }

    String r2 = basename(cut_r2,'.fastq.gz')

    command <<<
    eval "$(conda shell.bash hook)" 
    conda activate eprint
    fastqc -t 2 --extract -k 7 ~{cut_r2} -o .
    >>>
    runtime {
        cpu: 3
        memory: "5 GB"
    }
}


task FastQ_sort {
    input {
        File cut_r2
    }

    String sorted_r2 = basename(cut_r2,'.fastq.gz')

    command <<<
    eval "$(conda shell.bash hook)" 
    conda activate eprint
    gzip -d ~{cut_r2} > ~{sorted_r2}.fastq
    fastq-sort --id ~{sorted_r2}.fastq > ~{sorted_r2}.fq
    >>>

    runtime {
        cpu: 6
        memory: "10 GB"
    }

    output {
        File fastq_sort_r2 = "~{sorted_r2}.fq"
     }
}


task STAR_rmRep {
    input {
        File sorted_cut_r2
        File hg19_dup_tar
    }

    String prefix = basename(sorted_cut_r2,'fq')

    command <<<
    mkdir RepElements
    tar -xf ~{hg19_dup_tar}
    eval "$(conda shell.bash hook)" 
    conda activate eprint
    STAR \
    --runMode alignReads \
    --runThreadN 40 \
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
    --readFilesIn ~{sorted_cut_r2}
    >>>
    runtime {
        cpu: 40
        memory: "60 GB"
    }
    output {
        File star_r2 = "~{prefix}Unmapped.out.mate1"
        File star_bam = "~{prefix}Aligned.out.bam"
    }
}

task FastQ_sort_STAR {
    input {
        File star_r2
        File star_bam
        
    }
    String sorted_r2 = basename(star_r2,'.Unmapped.out.mate1') + '_r2_.fq'

    command <<<
    eval "$(conda shell.bash hook)" 
    conda activate eprint
    fastq-sort --id ~{star_r2} > ~{sorted_r2}
    >>>
    
    runtime {
        cpu: 6
        memory: "10 GB"
    }
    
    output {
        File sorted_start_r2 = "~{sorted_r2}"
    }
}

task STAR_genome_map {
    input {
        File sorted_start_r2
        File star_bam
        File hg19_tar
    }
    String prefix = basename(sorted_start_r2,'_r2_.fq') + '_hg19.'

    command <<<
    mkdir HG_19_DIR
    tar -xf ~{hg19_tar}
    eval "$(conda shell.bash hook)" 
    conda activate eprint
    STAR \
    --runMode alignReads \
    --runThreadN 30 \
    --genomeDir  HG_19_DIR \
    --genomeLoad NoSharedMemory \
    --readFilesIn ~{sorted_start_r2} \
    --outSAMunmapped Within \
    --outFilterMultimapNmax 1 \
    --outFilterMultimapScoreRange 1 \
    --outFileNamePrefix ~{prefix} \
    --outSAMattributes All \
    --outSAMtype BAM Unsorted \
    --outFilterType BySJout \
    --outReadsUnmapped Fastx \
    --outFilterScoreMin 10 \
    --outSAMattrRGline ID:~{prefix} \
    --outStd Log \
    --alignEndsType EndToEnd \
    --outBAMcompression 10 \
    --outSAMmode Full
    >>>

    runtime {
        cpu: 30
        memory: "50 GB"
    }

    output {
        File star_hg19_r2 = "${prefix}.Unmapped.out.mate2"
        File star_hg19_bam = "${prefix}.Aligned.out.bam"
    }
}