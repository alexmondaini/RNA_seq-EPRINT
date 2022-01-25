version 1.0

workflow BigWig {
    input {
        Array[File] bams
        File chr_size
    }
    scatter (bam in bams) {
    call Wigs {
        input:
        bam = bam,
        chr_size = chr_size
    }
    }
}

task Wigs {
    input {
        File bam
        File chr_size
        String? direction
    }
    String bw_pos = basename(bam,'.bam') + '_norm_pos.bw'
    String bw_neg = basename(bam,'.bam') + '_norm_neg.bw'

    command <<<
    makebigwigfiles \
    --bw_pos  ~{bw_pos} \
    --bw_neg ~{bw_neg}  \
    --bam  ~{bam} \
    --genome ~{chr_size} \
    --direction ~{default="r" direction}
    >>>

    runtime {
        cpu: 4
        memory: "30 GB"
        docker: "brianyee/makebigwigfiles@sha256:8d67afc36e388aa12f1b6d2bed8ea3b6ddaa9ec4296a93d5fa9f31a5b1ff16d4"
    }

    output {
        File result_bw_pos = "${bw_pos}"
        File result_bw_neg = "${bw_neg}"
    }
}