version 1.0

workflow Call_Peaks {
    
    input {
        Array[File] samples
    }

    scatter (sample in samples) {
        
        call Sort_Index_Bam {
            input:
            sample = sample
        }

        call Clipper {
            input:
            call_peak_bam = Sort_Index_Bam.sorted_indexed_bam,
            call_peak_bai = Sort_Index_Bam.sorted_indexed_bai
        }

        }
}

task Sort_Index_Bam {
    
    input {
        File sample
    }

    String prefix = basename(sample,'_hg19.Aligned.out.bam')

    command <<<
    source /groups/cgsd/alexandre/miniconda3/etc/profile.d/conda.sh 
    conda activate eprint
    samtools sort -O bam -m 2G -@ 20 -o ~{prefix}.bam ~{sample} 
    samtools index -b ~{prefix}.bam
    >>>

    runtime {
        cpu: 20
        memory: "40 GB"
    }

    output {
        File sorted_indexed_bam = "~{prefix}.bam"
        File sorted_indexed_bai = "~{prefix}.bam.bai"
    }
}


task Clipper {
    
    input {
        File call_peak_bam
        File call_peak_bai
    }
    
    String bed_peak_intervals = basename(call_peak_bam,'.bam') + '.bed'
    
    command <<<
    clipper \
    --species hg19 \
    --bam ~{call_peak_bam} \
    --outfile ~{bed_peak_intervals}
    >>>
    
    runtime {
        cpu: 16
        memory: "60 GB"
        docker: "brianyee/clipper@sha256:094ede2a0ee7a6f2c2e07f436a8b63486dc4a072dbccad136b7a450363ab1876"
    }
}

