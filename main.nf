params.outputDir
params.inputDir
params.reference
params.regions = ""
params.extra = ""
params.vcfSplitSize = 2000000

def pair_vcfs_with_index(ch) {
    // ch is a flat channel of VCFs and index files
    // produced by a path glob like '*.vcf.gz*'
    // Output: a channel of tuples (vcf, index)
    def branched = ch.flatten().branch {
        index: it =~ /vcf.gz.(csi|tbi)$/
        vcf: it =~ /vcf.gz$/
    }
    return branched.vcf.concat(branched.index) |
        map { it -> tuple ( it.getSimpleName(), it ) } |
        groupTuple
}

def remove_duplicate_filepair_keys(primary_ch, secondary_ch) {
    // primary_ch and secondary_ch are channels produced by
    // `fromFilePairs`. If any keys are present in both channels,
    // remove the key from the secondary channel.
    // This will avoid duplicating work on a bam file if it has
    // a bai and a csi index.
    return primary_ch.concat(secondary_ch).groupTuple() |
        map { it -> tuple(it[0], it[1][0]) }
}

/* STEP 1 */    include { individual_calling }                                 from "./modules/somatypus_modules.nf"
/* STEP 2 */    include { alternative_calling }                                from "./modules/somatypus_modules.nf"
/* STEP 3 */    include { split_calls }                                        from "./modules/somatypus_modules.nf"
/* STEP 4 */    include { filter_calls }                                       from "./modules/somatypus_modules.nf"
/* STEP 5 */    include { indel_flag2 }                                        from "./modules/somatypus_modules.nf"
/* STEP 5 */    include { merge_indel_flagged }                                from "./modules/somatypus_modules.nf"
/* STEP 6 */    include { merge_filtered }                                     from "./modules/somatypus_modules.nf"
/* STEP 7 */    include { extract_indels }                                     from "./modules/somatypus_modules.nf"
/* STEP 7 */    include { merge_extracted_indels }                             from "./modules/somatypus_modules.nf"
/* STEPS 8-16 */include { genotype as genotype_merged_snvs_allele1 }           from "./modules/somatypus_modules.nf"
/* STEPS 8-16 */include { genotype as genotype_merged_snvs_allele2 }           from "./modules/somatypus_modules.nf"
/* STEPS 8-16 */include { genotype as genotype_merged_snvs_allele3 }           from "./modules/somatypus_modules.nf"
/* STEPS 8-16 */include { genotype as genotype_indel_excluded_snvs_allele1 }   from "./modules/somatypus_modules.nf"
/* STEPS 8-16 */include { genotype as genotype_indel_excluded_snvs_allele2 }   from "./modules/somatypus_modules.nf"
/* STEPS 8-16 */include { genotype as genotype_indel_excluded_snvs_allele3 }   from "./modules/somatypus_modules.nf"
/* STEPS 8-16 */include { genotype as genotype_indels }                        from "./modules/somatypus_modules.nf"
/* STEPS 8-16 */include { split_vcf_into_n_rows as split_merged_snvs }         from "./modules/somatypus_modules.nf"
/* STEPS 8-16 */include { split_vcf_into_n_rows as split_indel_excluded_snvs } from "./modules/somatypus_modules.nf"
/* STEPS 8-16 */include { merge_vcf as merge_genotyped_merged_snvs }           from "./modules/somatypus_modules.nf"
/* STEPS 8-16 */include { merge_vcf as merge_genotyped_indel_excluded_snvs }   from "./modules/somatypus_modules.nf"
/* STEP 17 */   include { merge_filter_indelflagged }                          from "./modules/somatypus_modules.nf"
/* STEP 17 */   include { finalise_snvs }                                      from "./modules/somatypus_modules.nf"
/* STEP 17 */   include { finalise_indels }                                    from "./modules/somatypus_modules.nf"

workflow {
    bam_bai_inputs   = Channel.fromFilePairs("${params.inputDir}/*.bam{,.bai}")
    bam_csi_inputs   = Channel.fromFilePairs("${params.inputDir}/*.bam{,.csi}")
    bam_inputs = remove_duplicate_filepair_keys(bam_csi_inputs, bam_bai_inputs)
    cram_inputs = Channel.fromFilePairs("${params.inputDir}/*.cram{,.crai}")
    inputs = bam_inputs.concat(cram_inputs)
    
    refGlob = "${params.reference}" + "{,.fai}"
    reference = Channel.fromFilePairs("${refGlob}", checkIfExists: true)
    
    // Start running the pipeline
    // Steps 1-4 can be run on each input in parallel
    indiv_calls   = individual_calling(inputs.combine(reference))
    alt_calls     = alternative_calling(inputs.combine(reference))
    split_result  = split_calls(indiv_calls.mix(alt_calls))
    filter_result = filter_calls(split_result)

    // Steps 5-8 begin to merge the results together 
    indel_flagged = indel_flag2(split_result).collect() | merge_indel_flagged
    snv_merged    = merge_filtered(filter_result.map { it -> it[1] } | collect,
                                   indel_flagged)
    indel_extracted = extract_indels(indiv_calls.mix(alt_calls))
    merged_indels = merge_extracted_indels(indel_extracted.map { it -> it[1] } | collect)

    // Steps 9-16 to do the genotyping
    bams = inputs.map { it -> it[1][0] }.collect().map { it -> tuple("bams", it) }
    bais = inputs.map { it -> it[1][1] }.collect().map { it -> tuple("bais", it) }
     
    merged_snvs_allele1 = split_merged_snvs(snv_merged[0] | pair_vcfs_with_index, "${params.vcfSplitSize}")
        .flatten()
        .map { it -> tuple ( it.baseName - ~/.gz/, it ) } | groupTuple
    merged_snvs_allele2 = snv_merged[1] | pair_vcfs_with_index
    merged_snvs_allele3 = snv_merged[2] | pair_vcfs_with_index
    indel_excluded_allele1 = split_indel_excluded_snvs(snv_merged[3] | pair_vcfs_with_index, "${params.vcfSplitSize}")
        .flatten()
        .map { it -> tuple ( it.baseName - ~/.gz/, it ) } | groupTuple
    indel_excluded_allele2 = snv_merged[4] | pair_vcfs_with_index
    indel_excluded_allele3 = snv_merged[5] | pair_vcfs_with_index
    indels = merged_indels | pair_vcfs_with_index

    genotyped_snvs1_tmp = genotype_merged_snvs_allele1(bams
        .combine(bais)
        .combine(reference)
        .combine(merged_snvs_allele1))
    genotyped_snvs1 = merge_genotyped_merged_snvs(
        Channel.value("MergedSNVs_allele1"),
        genotyped_snvs1_tmp.map{ it[1][0] } | collect,
        genotyped_snvs1_tmp.map{ it[1][1] } | collect)

    genotyped_snvs2 = genotype_merged_snvs_allele2(bams
      .combine(bais)
      .combine(reference)
      .combine(merged_snvs_allele2))

    genotyped_snvs3 = genotype_merged_snvs_allele3(bams
      .combine(bais)
      .combine(reference)
      .combine(merged_snvs_allele3))

    genotyped_iex1_tmp = genotype_indel_excluded_snvs_allele1(bams
      .combine(bais)
      .combine(reference)
      .combine(indel_excluded_allele1))
    genotyped_iex1 = merge_genotyped_indel_excluded_snvs(
        Channel.value("IndelExcludedSNVs_allele1"),
        genotyped_iex1_tmp.map{ it[1][0] } | collect,
        genotyped_iex1_tmp.map{ it[1][1] } | collect)

    genotyped_iex2 = genotype_indel_excluded_snvs_allele2(bams
      .combine(bais)
      .combine(reference)
      .combine(indel_excluded_allele2))

    genotyped_iex3 = genotype_indel_excluded_snvs_allele3(bams
      .combine(bais)
      .combine(reference)
      .combine(indel_excluded_allele3))

    genotyped_indels = genotype_indels(bams
      .combine(bais)
      .combine(reference)
      .combine(indels))

    // Steps 17-18 to merge the results
    iex = genotyped_iex1.mix(genotyped_iex2).mix(genotyped_iex3)
        | map { it[1] } 
        | flatten 
        | branch {
            index: it =~ /vcf.gz.(csi|tbi)$/
            vcf: it =~ /vcf.gz$/
        } 
    //     .combine(genotyped_iex2 | map { it[1] })
    //     .combine(genotyped_iex3 | map { it[1] })
    // brn = tmp
    //     .branch{
    //         index: it =~ /vcf.gz.(csi|tbi)$/
    //         vcf: it =~ /vcf.gz$/
    //     }
    filtered = merge_filter_indelflagged(iex.vcf.collect(), iex.index.collect())

    snvs = genotyped_snvs1.mix(genotyped_snvs2).mix(genotyped_snvs3).mix(filtered)
        | map { it[1] } 
        | flatten 
        | branch {
            index: it =~ /vcf.gz.(csi|tbi)$/
            vcf: it =~ /vcf.gz$/
        }
    indels = genotyped_indels
        | map { it[1] } 
        | flatten 
        | branch {
            index: it =~ /vcf.gz.(csi|tbi)$/
            vcf: it =~ /vcf.gz$/
        }
    finalise_snvs(snvs.vcf.collect(), snvs.index.collect())
    finalise_indels(indels.vcf, indels.index)
}
