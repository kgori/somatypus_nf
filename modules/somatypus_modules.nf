def handle_arg(arg, prefix) {
    if (arg) {
        return "${prefix}${arg}"
    } else {
        return ""
    }
}

regions = handle_arg(params.regions, "--regions=")
extra = handle_arg(params.extra, "")

if (extra =~ "--logFileName|--refFile|--bamFiles|--regions|--minPosterior|--minReads|--minFlank|--trimReadFlank|--source|--getVariantsFromBAMs|--nCPU|--output|-o") {
    error("ERROR: Additional Platypus options cannot include --logFileName, --refFile, --bamFiles, --regions, --minPosterior, --minReads, --minFlank, --trimReadFlank, --source, --getVariantsFromBAMs, --nCPU, --output, or -o.")
}

process individual_calling {
    input:
    tuple val(id), file(alignmentFiles), val(ref), file(referenceFiles)

    output:
    tuple val(id), file("platypusVariants_${id}_default.vcf")

    script:
    """
    platypus callVariants \
        --logFileName=${id}_default.log \
        --refFile=${referenceFiles[0]} \
        --bamFiles=${alignmentFiles[0]} \
        --minPosterior=0 \
        --minReads=3 \
        --nCPU=${task.cpus} \
        ${regions} \
        ${extra} \
        -o platypusVariants_${id}_default.vcf
    """
}

process alternative_calling {
    input:
    tuple val(id), file(alignmentFiles), val(ref), file(referenceFiles)

    output:
    tuple val(id), file("platypusVariants_${id}_alternative.vcf")

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
        ${regions} \
        ${extra} \
        -o platypusVariants_${id}_alternative.vcf
    """
}

process split_calls {
    input:
    tuple val(id), file(vcf)

    output:
    tuple val(id), file("${vcf.baseName}.split.vcf")

    script:
    """
    Somatypus_SplitMA-MNVs.py ${vcf} > ${id}_split.log
    """
}

process filter_calls {
    input:
    tuple val(id), file(vcf)

    output:
    tuple val(id), file("${vcf.baseName}.filtered.vcf")

    script:
    """
    awk '!((\$7 ~ /badReads/) || (\$7 ~ /MQ/) || (\$7 ~ /strandBias/) || (\$7 ~ /SC/) || (\$7 ~ /QD/))' ${vcf} \
      > ${vcf.baseName}.filtered.vcf
    """
}

process indel_flag2 {
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
    input:
    path(vcf)
    path(filterlist)

    output:
    path "merged/MergedSNVs_allele1.vcf.gz*"
    path "merged/MergedSNVs_allele2.vcf.gz*"
    path "merged/MergedSNVs_allele3.vcf.gz*"
    path "merged/IndelExcludedSNVs_allele1.vcf.gz*"
    path "merged/IndelExcludedSNVs_allele2.vcf.gz*"
    path "merged/IndelExcludedSNVs_allele3.vcf.gz*"

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
        tabix --csi \${f}.gz
    done
    """
}

process extract_indels {
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
    input:
    path vcf

    output:
    path "MergedIndels.vcf.gz*"

    script:
    """
    Somatypus_IndelMerge2.py merge ${vcf} MergedIndels.vcf
    bgzip MergedIndels.vcf
    tabix --csi MergedIndels.vcf.gz
    """
}

process genotype {
    input:
    tuple val(bamID), path (alignmentFiles),
      val(baiID), path (alignmentIndexFiles),
      val(refID), path (referenceFiles),
      val(vcfID), path(source_vcf)

    output:
    tuple val(vcfID), path("${source_vcf[0].getSimpleName()}.genotyped.vcf.gz*")

    script:
    """
    ls -1 ${alignmentFiles} > bamlist.txt
    Somatypus_CreateGenotypingRegions.py -o regions.txt -w 1000 ${source_vcf[0]}
    cat regions.txt | tr ':' '\t' | tr '-' '\t' | sort -c -k1,1 -k2n,2n -k3n,3n

    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ${source_vcf[0]} > positions_expected.txt

    platypus callVariants \
        --logFileName=${source_vcf[0].getSimpleName()}.log \
        --refFile=${referenceFiles[0]} \
        --bamFiles=bamlist.txt \
        --minPosterior=0 \
        --nCPU=${task.cpus} \
        --minReads=3 \
        --source=${source_vcf[0]} \
        --getVariantsFromBAMs=0 \
        --regions=regions.txt \
        -o ${source_vcf[0].getSimpleName()}.genotyped.first.vcf

    bgzip ${source_vcf[0].getSimpleName()}.genotyped.first.vcf
    tabix --csi ${source_vcf[0].getSimpleName()}.genotyped.first.vcf.gz

    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
      ${source_vcf[0].getSimpleName()}.genotyped.first.vcf.gz \
      > positions_obtained.txt

    grep -vxFf positions_obtained.txt positions_expected.txt > positions_missing.txt || [ \$? -eq 1 ]

    if [ -s positions_missing.txt ]; then
        echo -e "\nGenotyping missing calls\n"

        bcftools isec --complement -w1 -Oz -o ${source_vcf[0].getSimpleName()}.second_source.vcf.gz ${source_vcf[0]} ${source_vcf[0].getSimpleName()}.genotyped.first.vcf.gz
        tabix --csi ${source_vcf[0].getSimpleName()}.second_source.vcf.gz
        Somatypus_CreateGenotypingRegions.py -o missing_regions.txt -w 0 ${source_vcf[0].getSimpleName()}.second_source.vcf.gz

        platypus callVariants \
        --logFileName=${source_vcf[0].getSimpleName()}.second.log \
        --refFile=${referenceFiles[0]} \
        --bamFiles=bamlist.txt \
        --regions=missing_regions.txt \
        --minPosterior=0 \
        --nCPU=${task.cpus} \
        --minReads=3 \
        --source=${source_vcf[0].getSimpleName()}.second_source.vcf.gz \
        --getVariantsFromBAMs=0 \
        -o ${source_vcf[0].getSimpleName()}.genotyped.second.vcf

        bgzip ${source_vcf[0].getSimpleName()}.genotyped.second.vcf
        tabix --csi ${source_vcf[0].getSimpleName()}.genotyped.second.vcf.gz

        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ${source_vcf[0].getSimpleName()}.genotyped.second.vcf.gz > positions_obtained_2ndpass.txt
    fi

    cat positions_obtained* \
        | sort -k1,1 -k2,2n -k3,3 -k4,4 \
        | uniq > positions_obtained_final.txt

    grep -vxFf positions_obtained_final.txt positions_expected.txt > positions_still_missing.txt || exit_code=\$?

    if [[ \$exit_code -ne 1 ]]; then
        echo -e "\nERROR: Some positions are still missing from the genotyped VCF. Please check the log files.\n"
        exit 1
    fi

    bcftools view -h ${source_vcf[0].getSimpleName()}.genotyped.first.vcf.gz \
        > header.txt

    grep -h -v "^#" ${source_vcf[0].getSimpleName()}.genotyped.*.vcf \
        | awk '!((\$7 ~ /badReads/) || (\$7 ~ /MQ/) || (\$7 ~ /strandBias/) || (\$7 ~ /SC/) || (\$7 ~ /QD/))' \
        | sort -k1,1 -k2,2n -k4,4 -k5,5 \
        > variants.txt
    cat header.txt variants.txt | bgzip -c > ${source_vcf[0].getSimpleName()}.genotyped.vcf.gz
    tabix --csi ${source_vcf[0].getSimpleName()}.genotyped.vcf.gz
    rm header.txt variants.txt
    bgzip *.first.vcf 
    if [ -s ${source_vcf[0].getSimpleName()}.genotyped.second.vcf ]; then
        bgzip *.second.vcf
    fi
    """
}

process split_vcf_into_n_files {
    input:
    tuple val(id), path(vcf)
    val n

    output:
    path "split_out/${vcf[0].getSimpleName()}-*"

    script:
    """
    mkdir split_out
    Somatypus_SplitVCF.py --vcf ${vcf[0]}  -n ${n} -s Files -o split_out/${vcf[0].getSimpleName()}
    """
}

process split_vcf_into_n_rows {
    input:
    tuple val(id), path(vcf)
    val r

    output:
    path "split_out/${vcf[0].getSimpleName()}-*"

    script:
    """
    mkdir split_out
    Somatypus_SplitVCF.py --vcf ${vcf[0]}  -r ${r} -s Rows -o split_out/${vcf[0].getSimpleName()}
    """
}

process merge_vcf {
    label 'bcfmerge'

    input:
    val(id)
    path vcf_files
    path index_files

    output:
    tuple val(id), path("${id}.vcf.gz*")

    script:
    """
    bcftools concat --threads ${task.cpus} -a -D -Oz -o ${id}.vcf.gz ${vcf_files}
    bcftools index --csi ${id}.vcf.gz
    """
}

process merge_filter_indelflagged {
    label 'bcfmerge'

    input:
    path vcf_files
    path index_files

    output:
    tuple val("merged.VAFfilt"), path("merged.VAFfilt.vcf.gz*")

    script:
    """
    bcftools concat --threads ${task.cpus} -a -D -Ov -o merged.vcf ${vcf_files}
    Somatypus_IndelRescuedFilter.py merged.vcf > VAFfilter.log
    bgzip merged.VAFfilt.vcf
    tabix --csi merged.VAFfilt.vcf.gz
    rm merged.vcf
    """
}

process finalise_snvs {
    label 'bcfmerge'

    input:
    path vcf_files
    path index_files

    output:
    path "Somatypus_SNVs_final.vcf.gz*"

    publishDir "${params.outputDir}", mode: 'copy'

    script:
    """
    bcftools concat --threads ${task.cpus} -a -D -Ov -o merged.vcf ${vcf_files}
    Somatypus_VAFfilter.py merged.vcf > VAFfilter.log
    bgzip -c merged.VAFfilt.vcf > Somatypus_SNVs_final.vcf.gz
    tabix --csi Somatypus_SNVs_final.vcf.gz
    rm merged.vcf merged.VAFfilt.vcf
    """
}

process finalise_indels {
    input:
    path vcf_file
    path index_file

    output:
    path "Somatypus_Indels_final.vcf.gz*"

    publishDir "${params.outputDir}", mode: 'copy'

    script:
    """
    gunzip -c ${vcf_file} > ${vcf_file.baseName}
    Somatypus_VAFfilter.py ${vcf_file.baseName} > VAFfilter.log
    bgzip -c "${vcf_file.getBaseName(2)}.VAFfilt.vcf" > Somatypus_Indels_final.vcf.gz
    tabix --csi Somatypus_Indels_final.vcf.gz
    rm ${vcf_file.baseName} ${vcf_file.getBaseName(2)}.VAFfilt.vcf
    """
}
