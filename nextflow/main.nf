params.reads          = "fastq/*.fastq.gz"
params.refFa          = "reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
params.refGtf         = "reference/Homo_sapiens.GRCh38.91.gtf"
params.outdir         = "results"
params.useAnnotations = true
params.minReads       = 10000   // QC: minimum total reads
params.minMapRate     = 50      // QC: minimum % reads mapped

process ALIGN {
    tag "$sample"
    cpus 8
    input:
    tuple val(sample), path(fastq)
    path fa
    output:
    tuple val(sample), path("${sample}.sam")
    script:
    def flags = sample.contains("directRNA") ? "-ax splice -uf -k14" : "-ax splice"
    """
    minimap2 ${flags} -t ${task.cpus} ${fa} ${fastq} > ${sample}.sam
    """
}

process SAM2BAM {
    tag "$sample"
    cpus 8
    publishDir "${params.outdir}/bam", mode: 'copy'
    input:
    tuple val(sample), path(sam)
    output:
    tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai")
    script:
    """
    samtools view -b -@ ${task.cpus} ${sam} | samtools sort -@ ${task.cpus} -o ${sample}.bam
    samtools index ${sample}.bam
    """
}

process QC {
    tag "$sample"
    publishDir "${params.outdir}/qc", mode: 'copy', pattern: "*.tsv"
    input:
    tuple val(sample), path(bam), path(bai)
    output:
    tuple val(sample), path(bam), path(bai), path("${sample}.status"), emit: checked
    path "${sample}_qc.tsv", emit: report
    script:
    """
    total=\$(samtools view -c -F 0x900 ${bam})
    mapped=\$(samtools view -c -F 0x904 ${bam})
    avglen=\$(samtools stats ${bam} | grep '^SN' | grep 'average length' | cut -f3)
    pct=\$(awk -v m=\$mapped -v t=\$total 'BEGIN{ if (t>0) printf "%.2f", 100*m/t; else print 0 }')
    status=\$(awk -v t=\$total -v p=\$pct -v mr=${params.minReads} -v mp=${params.minMapRate} 'BEGIN{ print (t>=mr && p>=mp) ? "PASS" : "FAIL" }')
    echo -e "sample\\ttotal_reads\\tmapped_reads\\tpct_mapped\\tavg_read_length\\tstatus" > ${sample}_qc.tsv
    echo -e "${sample}\\t\$total\\t\$mapped\\t\$pct\\t\$avglen\\t\$status" >> ${sample}_qc.tsv
    echo \$status > ${sample}.status
    """
}

process BAMBU {
    cpus 1
    publishDir "${params.outdir}/bambu", mode: 'copy'
    input:
    path bams
    path bais
    path fa
    path gtf
    output:
    path "*.{gtf,txt,rds}"
    script:
    """
    run_bambu.R ${params.useAnnotations} ${fa} ${gtf} ${bams}
    """
}

workflow {
    reads = channel.fromPath(params.reads, checkIfExists: true).map { f -> tuple(f.simpleName, f) }
    fa  = file(params.refFa,  checkIfExists: true)
    gtf = file(params.refGtf, checkIfExists: true)

    ALIGN(reads, fa)
    SAM2BAM(ALIGN.out)
    QC(SAM2BAM.out)
    QC.out.report.collectFile(name: 'qc_summary.tsv', keepHeader: true, skip: 1, storeDir: "${params.outdir}/qc")

    passed = QC.out.checked.filter { sample, bam, bai, status -> status.text.trim() == 'PASS' }
    BAMBU(passed.map { sample, bam, bai, status -> bam }.collect(),
          passed.map { sample, bam, bai, status -> bai }.collect(), fa, gtf)
}
