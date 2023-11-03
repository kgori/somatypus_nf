process individual_calling {
    cpus { 4 * 2**(task.attempt - 1) }
    executor 'lsf'
    queue 'normal'
    memory { 8.GB * 2**(task.attempt - 1) }
    time 12.h
    maxRetries 1

    input:
    tuple val(id), file(alignmentFiles), val(ref), file(referenceFiles)

    output:
    tuple val(id), file("platypusVariants_${id}_default.vcf")

    // publishDir 

    // # TODO: add $REGIONSARG and $EXTRA
    script:
    """
    platypus callVariants \
        --logFileName=${id}_default.log \
        --refFile=${referenceFiles[0]} \
        --bamFiles=${alignmentFiles[0]} \
        --minPosterior=0 \
        --minReads=3 \
        --nCPU=${task.cpus} \
        -o platypusVariants_${id}_default.vcf
    """
}

process alternative_calling {
    cpus { 4 * 2**(task.attempt - 1) }
    executor 'lsf'
    queue 'normal'
    memory { 8.GB * 2**(task.attempt - 1) }
    time 12.h
    maxRetries 1

    input:
    tuple val(id), file(alignmentFiles), val(ref), file(referenceFiles)

    output:
    tuple val(id), file("platypusVariants_${id}_alternative.vcf")

    // publishDir 

    // # TODO: add $REGIONSARG and $EXTRA
    script:
    """
    platypus callVariants \
        --logFileName=${id}_alternative.log \
        --refFile=${referenceFiles[0]} \
        --bamFiles=${alignmentFiles[0]} \
        --minPosterior=0 \
        --minReads=3 \
        --minFlank=0 \
        --trimReadFlank=10 \
        --nCPU=${task.cpus} \
        -o platypusVariants_${id}_alternative.vcf
    """
}

process split_calls {
    cpus 1
    executor 'lsf'
    queue 'normal'
    memory 1.GB
    time 2.h

    input:
    tuple val(id), file(vcf)

    output:
    tuple val(id), file("${vcf.baseName}.split.vcf")

    // publishDir

    script:
    """
    Somatypus_SplitMA-MNVs.py ${vcf} > ${id}_split.log
    """
}

process filter_calls {
    cpus 1
    executor 'lsf'
    queue 'normal'
    memory 1.GB
    time 2.h

    input:
    tuple val(id), file(vcf)

    output:
    tuple val(id), file("${vcf.baseName}.filtered.vcf")

    // publishDir

    script:
    """
    awk '!((\$7 ~ /badReads/) || (\$7 ~ /MQ/) || (\$7 ~ /strandBias/) || (\$7 ~ /SC/) || (\$7 ~ /QD/))' ${vcf} \
      > ${vcf.baseName}.filtered.vcf
    """
}

process indel_flag {
    cpus 1
    executor 'lsf'
    queue 'normal'
    memory { 8.GB * 2**(task.attempt - 1) }
    time 8.h
    maxRetries 2

    publishDir "${params.outputDir}/indel_flagged_SNVs"

    input:
    tuple val(id), path(vcf)

    output:
    file("${vcf.baseName}_indel_flagged_SNVs.txt")

    script:
    """
    ls -1 ${vcf} > list.txt
    Somatypus_IndelFlag.py list.txt ${vcf.baseName}_indel_flagged_SNVs.txt > ${vcf.baseName}_indel_flagged_SNVs.log
    """
}

process indel_flag2 {
    cpus 1
    executor 'lsf'
    queue 'normal'
    memory { 4.GB * 2**(task.attempt - 1) }
    time 2.h
    maxRetries 2

    input:
    tuple val(id), path(vcf)

    output:
    file("${vcf.baseName}_indel_flagged_SNVs.txt")

    script:
    """
    Somatypus_IndelFlag2.py ${vcf} ${vcf.baseName}_indel_flagged_SNVs.txt > ${vcf.baseName}_indel_flagged_SNVs.log
    sort ${vcf.baseName}_indel_flagged_SNVs.txt > tmp && mv tmp ${vcf.baseName}_indel_flagged_SNVs.txt
    """
}

process merge_indel_flagged {
    cpus 1
    executor 'lsf'
    queue 'normal'
    memory 16.GB
    time 2.h
    
    input:
    path(txt)

    output:
    file("indel_flagged_SNVs.txt")

    script:
    """
    sort -u -m ${txt} > indel_flagged_SNVs.txt
    """
}

process merge_filtered {
    cpus 1
    executor 'lsf'
    queue 'normal'
    memory 16.GB
    time 12.h

    input:
    path(vcf)
    path(filterlist)

    output:
    path "merged/*.vcf.gz"

    publishDir "${params.outputDir}/merged"

    script:
    """
    ls -1 ${vcf} > list.txt
    mkdir -p merged
    Somatypus_SNVmerge.py list.txt ${filterlist} merged > merged.log

    for f in merged/*.vcf; do
        cat <(echo "##fileformat=VCFv4.1") \
            <(echo "##source=Somatypus") \
            <(echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO") \
            <(sort -k1,1 -k2,2n -k4,4 -k5,5 \${f}) > \${f}.tmp
        mv \${f}.tmp \${f}
        bgzip \${f}
        tabix \${f}.gz
    done
    """
}

process extract_indels {
    cpus 1
    executor 'lsf'
    queue 'normal'
    memory 4.GB
    time 2.h

    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), file("${vcf.baseName}.extract.vcf")

    script:
    """
    Somatypus_IndelMerge2.py extract ${vcf} ${vcf.baseName}.extract.vcf > ${vcf.baseName}.extract.log
    """
}

process merge_extracted_indels {
    cpus 1
    executor 'lsf'
    queue 'normal'
    memory 16.GB
    time 2.h

    input:
    path vcf

    output:
    file "MergedIndels.vcf.gz"

    script:
    """
    Somatypus_IndelMerge2.py merge ${vcf} MergedIndels.vcf
    bgzip MergedIndels.vcf
    tabix MergedIndels.vcf.gz
    """
}

process genotype {
    // cpus 8
    // executor 'lsf'
    // queue 'long'
    // memory 48.GB
    // time 48.h

    input:
    tuple val(id), path(alignments), path(vcf)
    path reference

}
