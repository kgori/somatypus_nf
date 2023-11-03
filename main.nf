params.outputDir
params.inputDir
params.reference
params.regions = ""
params.extra = ""

/* STEP 1 */include { individual_calling }     from "./modules/somatypus_modules.nf"
/* STEP 2 */include { alternative_calling }    from "./modules/somatypus_modules.nf"
/* STEP 3 */include { split_calls }            from "./modules/somatypus_modules.nf"
/* STEP 4 */include { filter_calls }           from "./modules/somatypus_modules.nf"
/* STEP 5 */include { indel_flag2 }            from "./modules/somatypus_modules.nf"
/* STEP 5 */include { merge_indel_flagged }    from "./modules/somatypus_modules.nf"
/* STEP 6 */include { merge_filtered }         from "./modules/somatypus_modules.nf"
/* STEP 7 */include { extract_indels }         from "./modules/somatypus_modules.nf"
/* STEP 7 */include { merge_extracted_indels } from "./modules/somatypus_modules.nf"

workflow {
    bam_inputs = Channel.fromFilePairs("${params.inputDir}/*.bam{,.bai}")
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
}
