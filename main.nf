#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ===================================
// Global configuration and parameters
// ===================================

log.info """\
         PRIMATE CONSTRAINT PIPELINE
         ==========================
         HAL file: ${params.hal_file}
         Roulette table: ${params.roulette_table}
         Species file: ${params.species_file}
         Output directory: ${params.outdir}
         Max memory: ${params.max_memory}
         Max CPUs: ${params.max_cpus}
         Max time: ${params.max_time}
         """
         .stripIndent()

// ===================================
// Data preparation processes
// ===================================

process PREP_FILES_BED {
    executor 'local'
    cpus 1
    memory '4 GB'
    time '1h'

    output:
    path 'Homo_sapiens.GRCh38.113.chr.gff3', emit: gff_file
    path 'Homo_sapiens.GRCh38.regulatory_features.v113.gff3', emit: regulatory_gff_file

    script:
    """
    bash ${projectDir}/scripts/prep_files_bed.sh
    """
}

process MAKE_BED {
    executor 'local'
    cpus 10
    memory '40 GB'
    time '1h'

    input:
    path gff_file
    path regulatory_gff_file

    output:
    path 'features_by_loci.bed', emit: bed_output

    script:
    """
    python ${projectDir}/scripts/make_bed.py \
        --gff_file $gff_file \
        --regulatory_gff_file $regulatory_gff_file \
        --output features_by_loci.bed
    """
}

process MAKE_BED_PT2 {
    executor 'local'
    cpus {params.max_cpus}
    memory {params.max_memory}
    time {params.max_time}

    output:
    path('features_by_loci.bed'), emit: bed_output

    script:
    """
    bash ${projectDir}/scripts/generate_mask.sh
    """
}

process CHECK_BASESPACE_AUTH {
    executor 'local'
    errorStrategy 'finish' 

    output:
    val('authorize'), emit: auth

    script:
    """
    #!/bin/bash -ue
    if [ -d "\$HOME/.basespace" ]; then
        exit 0
    else
        echo "BaseSpace CLI is not authenticated. Please run 'bs auth'" > auth_check.txt
        exit 2
    fi
    """
}

// ===================================
// VCF handling processes (VCF2MT module)
// ===================================

process GET_VCF_IDS {
    executor 'local'

    input:
    val auth

    output:
    path "vcf_ids.txt", emit: vcf_ids_file

    script:
    """
    echo $HOME
    bash ${params.scriptsDir}/get_vcf_ids.sh
    """
}

process DOWNLOAD_VCFS {
    storeDir "${projectDir}/vcfs"
    executor 'local'
    cpus 1
    memory '1 GB'
    time '1h'
    maxForks 10
    
    input:
    tuple val(id), val(file)

    output:
    path "${file}", emit: vcf_file

    script:
    """
    bash ${params.scriptsDir}/download_vcf.sh $id .
    """
}

process GET_CHROM_FORMAT_VCF {
    executor 'local'
    maxForks 6

    input:
    tuple val(species), path(vcf), path(alias)

    output:
    tuple val(species), path(vcf), path(alias), env(format), emit: format
    
    script:
    """
    format=\$(python ${params.scriptsDir}/getformat.py $vcf $alias)
    """
}

process VCF2MT {
    storeDir "${projectDir}/mts"
    cpus { 12 }
    memory { '40 G' }
    maxRetries 2

    input:
    tuple val(species), path(vcf_file), val(vcf_format), path(alias_file), val(hal_format), path(primate_fasta), path(primate_index)

    output:
    tuple val(species), path("${vcf_file.baseName.replaceAll(/\.SNV\.vcf$/, '')}.mt"), emit: mt

    script:
    vcf_bname=vcf_file.baseName.replaceAll(/\.SNV\.vcf$/, '')
    """
    local_tmp=\$(mktemp -d \${PWD}/tmp_XXXXXX)
    memory=\$(echo "${task.memory}" | awk '/GB/ { gsub(" GB", "", \$1); print \$1 "g" } /TB/ { gsub(" TB", "", \$1); printf "%.0f", \$1 * 1000; print "g" }')
    python3 ${params.scriptsDir}/vcf2mt.py --cpus "${task.cpus}" \
                                           --species "${species}" \
                                           --memory "\${memory}" \
                                           --tmpdir "\$local_tmp" \
                                           --vcf_file "${vcf_file}" \
                                           --alias_file "${alias_file}" \
                                           --vcf_format "${vcf_format}" \
                                           --hal_format "${hal_format}" \
                                           --primate_fasta_path "${primate_fasta}" \
                                           --primate_index_path "${primate_index}" \
                                           --output .
    sleep 60
    trap "rm -rf \${local_tmp}" EXIT
    """
}

// ===================================
// Reference and alias handling processes 
// ===================================

process MAKE_ALIAS_MANIFEST {
    executor 'local'

    output:
    path "alias_manifest.txt", emit: alias_manifest

    script:
    """
    bash ${params.scriptsDir}/make_alias_manifest.sh
    """
}

process GET_ALIAS_FILES {
    executor 'local'
    maxForks 5

    input:
    tuple val(species), val(gca_id)

    output:
    tuple val(species), path("*.txt"), emit: alias_file

    script:
    """
    bash ${params.scriptsDir}/get_alias_file.sh ${species} ${gca_id} .
    """
}

process MAKE_SAMPLE2REF_MAP {
    executor 'local'

    output:
    path "sample2ref.txt", emit: sample2ref_file

    script:
    """
    bash ${params.scriptsDir}/make_sample2ref_map.sh
    """
}

process GET_REFERENCE {
    executor 'local'
    maxForks 3
    cpus 1
    memory '4 GB'
    time '1h'

    input:
    tuple val(species), path(hal)

    output:
    tuple val(species), path("*.fa"), path("*.fa.fai"), emit: fasta

    script:
    """
    bash ${params.scriptsDir}/getreference.sh $species $hal    
    """
}

process GET_REFERENCE_HUMAN {
    executor 'local'
    maxForks 3
    cpus 1
    memory '4 GB'
    time '1h'

    input:
    path(hal)

    output:
    tuple val('Homo_sapiens'), path("*.fa"), path("*.fa.fai"), emit: fasta

    script:
    """
    bash ${params.scriptsDir}/getreference.sh Homo_sapiens $hal
    """
}

process GET_CHROM_FORMAT_HAL {
    executor 'local'
    maxForks 10

    input:
    tuple val(species), path(alias), path(hal_)

    output:
    tuple val(species), env(format), emit: format

    script:
    """
    format=\$(python ${params.scriptsDir}/getformat_hal.py $species $alias $hal_)
    """
}

process GET_NUM_CONTIGS {
    executor 'local'
    maxForks 1
    memory '40 GB'
    cpus 10

    input:
    tuple val(species), path(fasta), path(indexes)

    output:
    tuple val(species), env(contigs), emit: contig_count

    script:
    """
    local_tmp=\$(mktemp -d \${PWD}/tmp_XXXXXX) 
    memory=\$(echo "${task.memory}" | awk '/GB/ { gsub(" GB", "", \$1); print \$1 "g" } /TB/ { gsub(" TB", "", \$1); printf "%.0f", \$1 * 1000; print "g" }')
    contigs=\$(python ${params.scriptsDir}/get_contig_count.py --cpus "${task.cpus}" --memory "\${memory}" --tmpdir "\$local_tmp" --species $species --species_fasta_path $fasta --species_index_path $indexes)
    sleep 60
    trap "rm -rf \${local_tmp}" EXIT
    """
}

// ===================================
// VEP processing module
// ===================================

process MAKE_POS_BED {
    executor 'local'
    cpus 1
    memory '4 GB'
    time '1h'
    maxRetries 3

    input:
    path(bed_file)

    output:
    path("features_by_position.parquet"), emit: bedxpos

    script:
    """
    python ${params.scriptsDir}/make_pos_bed.py ${bed_file} .
    """
}

process MAKE_HG38_BEDS {
    executor 'local'
    cpus 1
    memory '4 GB'
    time '1h'
    
    input:
    path(bedxpos)
    
    output:
    path("hg38_*_part_*.bed"), emit:hg38beds

    script:
    """
    python ${params.scriptsDir}/make_hg38_beds.py $bedxpos hg38
    """
}

process MAKE_VEP_CONFIG {
    cpus 1
    memory '2 GB'
    time '30min'
    maxRetries 3
    
    output:
    path("vep_config.json"), emit: config

    script:
    """
    ${params.scriptsDir}/vep_config.sh 
    """
}

process MAKE_REGIONSXVAR {
    storeDir "${projectDir}/regionsxvar"
    maxForks 1
    cpus { params.max_cpus }
    memory { params.max_memory }
    time { params.max_time }
    maxRetries 2

    input:
    tuple path(position_table_path), path(human_fasta_path), path(human_index_path)

    output:
    path("regionsxvar.ht"), emit: regionsxvar

    script:
    """
    local_tmp=\$(mktemp -d \${PWD}/tmp_XXXXXX) 
    memory=\$(echo "${task.memory}" | awk '/GB/ { gsub(" GB", "", \$1); print \$1 "g" } /TB/ { gsub(" TB", "", \$1); printf "%.0f", \$1 * 1000; print "g" }')
    python3 ${params.scriptsDir}/regionsxvar.py --cpus "${task.cpus}" \
                                        --memory "\${memory}" \
                                        --tmpdir "\$local_tmp" \
                                        --position_table "${position_table_path}" \
                                        --human_fasta_path "${human_fasta_path}" \
                                        --human_index_path "${human_index_path}" \
                                        --output .
    sleep 60
    trap "rm -rf \${local_tmp}" EXIT
    """
}

process ANNOTATE_VEP {
    storeDir "${projectDir}/vep"
    maxForks 1
    container = "${projectDir}/annotation.sif"
    cpus { params.max_cpus }
    memory { params.max_memory }
    time { params.max_time }
    maxRetries 2

    input:
    tuple path(regionsxvar_path), path(vep_config)

    output:
    path("vep.ht"), emit: vep

    script:
    """
    local_tmp=\$(mktemp -d \${PWD}/tmp_XXXXXX) 
    memory=\$(echo "${task.memory}" | awk '/GB/ { gsub(" GB", "", \$1); print \$1 "g" } /TB/ { gsub(" TB", "", \$1); printf "%.0f", \$1 * 1000; print "g" }')
    python3 ${params.scriptsDir}/vep.py --cpus "${task.cpus}" \
                                        --memory "\${memory}" \
                                        --tmpdir "\$local_tmp" \
                                        --regionsxvar_table "${regionsxvar_path}" \
                                        --vep_config "${vep_config}" \
                                        --output .
    sleep 60
    trap "rm -rf \${local_tmp}" EXIT
    """
}

process PROCESS_VEP {  
    storeDir "${projectDir}/vep"
    cpus { params.max_cpus } 
    memory { params.max_memory } 
    time { '24h' }
    maxRetries 2

    input:
    path(vep_ht)

    output:
    path("vep_cdsinfo.ht"), emit: vep_processed

    script:
    """
    local_tmp=\$(mktemp -d \${PWD}/tmp_XXXXXX)
    memory=\$(echo "${task.memory}" | awk '/GB/ { gsub(" GB", "", \$1); print \$1 "g" } /TB/ { gsub(" TB", "", \$1); printf "%.0f", \$1 * 1000; print "g" }')    
    python3 ${params.scriptsDir}/processvep.py \
        --vep_ht "${vep_ht}" \
        --cpus "${task.cpus}" \
        --memory "\${memory}" \
        --tmpdir "\$local_tmp" \
        --output .
    sleep 60
    trap "rm -rf \${local_tmp}" EXIT
    """
}

// ===================================
// Liftover module
// ===================================

process MAKE_LIFTOVER_BED {
    executor 'local'
    maxForks 1
    cpus { params.max_cpus }
    memory { params.max_memory }
    time { params.max_time }
    maxRetries 0

    when:
    !file("${projectDir}/liftover/${species}").exists()

    input:
    tuple val(species), path(hal_file), path(bed_files)

    output:
    tuple val(species), path("*_final.bed"), emit: lifted_beds

    script:
    """
    local_tmp=\$(mktemp -d \${PWD}/tmp_XXXXXX)
    parallel --tmpdir "\$local_tmp" -j \$((${task.cpus} - 10)) bash ${params.scriptsDir}/liftover_bed.sh {} $species $hal_file -o . ::: $bed_files
    sleep 60
    trap "rm -rf \${local_tmp}" EXIT
    """
}

process MAKE_POSINFO {
    storeDir "${projectDir}/liftover/${species}"
    executor 'local'
    maxForks 1
    cpus { params.max_cpus }
    memory { params.max_memory }
    time { params.max_time }
    maxRetries 0

    input:
    tuple val(species), path(species_fasta_path), path(species_index_path), path(human_fasta_path), path(human_index_path), path(species_bed_path), path(regionsxvar_path)

    output:
    tuple val(species), path("${species}_pos_info.ht"), emit: pos_info

    script:
    """
    local_tmp=\$(mktemp -d \${PWD}/tmp_XXXXXX) 
    memory=\$(echo "${task.memory}" | awk '/GB/ { gsub(" GB", "", \$1); print \$1 "g" } /TB/ { gsub(" TB", "", \$1); printf "%.0f", \$1 * 1000; print "g" }')    
    python3 ${params.scriptsDir}/bed2position.py \
        --species "$species" \
        --human_fasta_path "$human_fasta_path" \
        --human_index_path "$human_index_path" \
        --species_fasta_path "$species_fasta_path" \
        --species_index_path "$species_index_path" \
        --species_bed_path "$species_bed_path" \
        --regionsxvar_table "$regionsxvar_path" \
        --cpus "${task.cpus}" \
        --memory "\${memory}" \
        --tmpdir "\$local_tmp" \
        --output "\${PWD}" 
    sleep 60
    trap "rm -rf \${local_tmp}" EXIT
    """ 
}

process MAKE_LIFTOVER_TABLE {
    storeDir "${projectDir}/liftover/${species}"
    executor 'local'
    maxForks 1
    cpus { params.max_cpus }
    memory { params.max_memory }
    time { params.max_time }
    maxRetries 0

    input:
    tuple val(species), path(pos_info), path(vep_path), path(species_fasta_path), path(species_index_path)

    output:
    tuple val(species), path("${species}_var_info.ht"), emit: lifttable

    script:
    """
    local_tmp=\$(mktemp -d \${PWD}/tmp_XXXXXX) 
    memory=\$(echo "${task.memory}" | awk '/GB/ { gsub(" GB", "", \$1); print \$1 "g" } /TB/ { gsub(" TB", "", \$1); printf "%.0f", \$1 * 1000; print "g" }')    
    python3 ${params.scriptsDir}/pos2liftover.py \
        --species "$species" \
        --species_fasta_path "$species_fasta_path" \
        --species_index_path "$species_index_path" \
        --pos_ht "$pos_info" \
        --vep_annotations "$vep_path" \
        --cpus "${task.cpus}" \
        --memory "\${memory}" \
        --tmpdir "\$local_tmp" \
        --output "\${PWD}" 
    sleep 90
    trap "rm -rf \${local_tmp}" EXIT
    """ 
}

// ===================================
// Species merging module
// ===================================

process MERGEMTS {
    storeDir "${projectDir}/liftover/${species}"
    maxForks 6  // Increased from 1
    cpus { Math.min(32, params.max_cpus) }  // Limit to 32 cores per instance
    memory { '50 GB' }  
    time { params.max_time }
    maxRetries 2

    input:
    tuple val(species), path(mts), path(pos_ht), path(lifttable), val(contig_count), path(sample_mapping)

    output:
    tuple val(species), path("${species}_lifted_vars.ht"), path("${species}_colinfo_vars.ht"), emit: speciesht

    script:
    """
    local_tmp=\$(mktemp -d \${PWD}/tmp_XXXXXX) 
    memory=\$(echo "${task.memory}" | awk '/GB/ { gsub(" GB", "", \$1); print \$1 "g" } /TB/ { gsub(" TB", "", \$1); printf "%.0f", \$1 * 1000; print "g" }')
    python3 ${params.scriptsDir}/liftvars.py \
        --species "$species" \
        --sample_mts "$mts" \
        --species_ht "$lifttable" \
        --pos_ht "$pos_ht" \
        --contig_count "$contig_count" \
        --mapping "$sample_mapping" \
        --cpus "${task.cpus}" \
        --memory "\${memory}" \
        --tmpdir "\$local_tmp" \
        --outdir "\${PWD}"
    sleep 60
    trap "rm -rf \${local_tmp}" EXIT
    """
}

process MERGE_SPECIES_HTS {
    storeDir "${projectDir}/merged_vars"
    maxForks 1
    cpus { params.max_cpus }
    memory { params.max_memory }
    time { params.max_time }
    maxRetries 2

    input:
    tuple path(regionsxvar), path(human_fa), path(human_faidx), path(var_hts)

    output:
    path("aggregated_primate_counts.ht"), emit: all_species

    script:
    """
    local_tmp=\$(mktemp -d \${PWD}/tmp_XXXXXX) 
    memory=\$(echo "${task.memory}" | awk '/GB/ { gsub(" GB", "", \$1); print \$1 "g" } /TB/ { gsub(" TB", "", \$1); printf "%.0f", \$1 * 1000; print "g" }')
    python3 ${params.scriptsDir}/merge_all.py \
        --regionsxvar "$regionsxvar" \
        --human_fa "$human_fa" \
        --human_faidx "$human_faidx" \
        --var_hts "$var_hts" \
        --cpus "${task.cpus}" \
        --memory "\${memory}" \
        --tmpdir "\$local_tmp" \
        --outdir "\${PWD}"
    sleep 60
    trap "rm -rf \${local_tmp}" EXIT
    """
}

// ===================================
// Main workflow
// ===================================

workflow {
    CHECK_BASESPACE_AUTH()
    /*
    // Data preparation
    PREP_FILES_BED()
    
    MAKE_BED(
        PREP_FILES_BED.out.gff_file,
        PREP_FILES_BED.out.regulatory_gff_file
    )
    */
    
    MAKE_BED_PT2()

    // VCF handling
    GET_VCF_IDS(CHECK_BASESPACE_AUTH.out.auth.collect())
    vcf_ids_ch = GET_VCF_IDS.out.vcf_ids_file
        .splitText()
        .map { line ->
            def fields = line.trim().split('\t')
            return tuple(fields[0], fields[1])
        }
    DOWNLOAD_VCFS(vcf_ids_ch)
    
    // Alias and reference handling
    MAKE_ALIAS_MANIFEST()
    alias_manifest_ch = MAKE_ALIAS_MANIFEST.out.alias_manifest
         .splitCsv(sep: '\t')
         .map { line -> tuple(line[0].trim(), line[1].trim()) }
    GET_ALIAS_FILES(alias_manifest_ch)
    MAKE_SAMPLE2REF_MAP()
    
    // HAL file handling
    Channel.fromPath(params.hal_file)
        .set { hal }
    
    // Get sample to species mapping
    MAKE_SAMPLE2REF_MAP.out.sample2ref_file
        .map { path ->
            def species2ref = [:]
            path.text.eachLine { line ->
                def (id, ref_species) = line.split(/\s+/).collect { it.trim() }
                species2ref[id] = ref_species
            }
            return species2ref
        }
        .combine(DOWNLOAD_VCFS.out.vcf_file)
        .map { species2ref, vcf_file ->
            def id = vcf_file.getBaseName().tokenize( '.' )[0]
            def ref_species = species2ref[id]
            return tuple(ref_species, vcf_file)
        }
        .combine(GET_ALIAS_FILES.out.alias_file, by: 0)
        .set { sp_vcf_alias_ch }
    
    // VEP related processing
    //MAKE_POS_BED(MAKE_BED.out.bed_output)
    MAKE_POS_BED(MAKE_BED_PT2.out.bed_output)
    MAKE_HG38_BEDS(MAKE_POS_BED.out.bedxpos)
    
    // Extract species for reference generation
    sp_vcf_alias_ch
        .map { tuple -> tuple[0] }
        .unique()
        .set { sp_ch }
    
    sp_ch
        .combine(hal)
        .set { sp_hal_ch }
    
    sp_hal_ch
        .combine(MAKE_HG38_BEDS.out.hg38beds.toList())
        .set { sp_hal_hg38beds_ch }
      
    // Process liftover and references
    sp_hal_hg38beds_ch
    .filter { species, hal_file, bed_files ->
        file("${projectDir}/liftover/${species}/${species}_pos_info.ht").exists() &&
        file("${projectDir}/liftover/${species}/${species}_var_info.ht").exists()
    }
    .map { species, hal_file, bed_files ->
        def pos_ht = file("${projectDir}/liftover/${species}/${species}_pos_info.ht")
        def var_ht = file("${projectDir}/liftover/${species}/${species}_var_info.ht")
        tuple(species, pos_ht, var_ht)
    }
    .set { existing_lifttable_ch }

    sp_hal_hg38beds_ch
        .filter { species, hal_file, bed_files ->
            !file("${projectDir}/liftover/${species}/${species}_var_info.ht").exists()
        }
        | MAKE_LIFTOVER_BED

    
    GET_REFERENCE(sp_hal_ch)
    GET_REFERENCE_HUMAN(hal)
    
    // Get chromosome formats

    skip_formatting_ch = sp_vcf_alias_ch
        .map { species, vcf_file, alias_f ->
            // Check if the expected output .mt file for this VCF does not exist
            def mt_base_name = vcf_file.baseName.replaceAll(/\.SNV\.vcf$/, '')
            mt_file = file("${projectDir}/mts/${mt_base_name}.mt")
            done = mt_file.exists()
            return tuple(species, mt_base_name, mt_file,vcf_file, alias_f, done)
        }

    GET_CHROM_FORMAT_VCF(
        skip_formatting_ch
        .filter { species, mt_base_name, mt_file, vcf_file, alias_f, done -> !done }
        .map { species, mt_base_name, mt_file, vcf_file, alias_f, done ->
            tuple(species, vcf_file, alias_f)
        }
    )
    
    GET_ALIAS_FILES.out.alias_file
        .combine(hal)
        .set{ sp_alias_hal_ch }
    
    GET_CHROM_FORMAT_HAL(sp_alias_hal_ch)
    
    // Prepare for VCF2MT
    GET_CHROM_FORMAT_VCF.out.format
        .filter { it[3] != 'Other' }
        .map { species, vcf, chromAlias, format -> 
            tuple(species, vcf, format)
        }
        .combine(GET_ALIAS_FILES.out.alias_file, by: 0)
        .combine(GET_CHROM_FORMAT_HAL.out.format, by: 0)
        .combine(GET_REFERENCE.out.fasta, by: 0)
        .set{ sp_vcf_format_alias_halformat_ref_ch }
    
    // Run VCF2MT
    VCF2MT(sp_vcf_format_alias_halformat_ref_ch)
    
    // VEP processing
    MAKE_VEP_CONFIG()
    
    MAKE_POS_BED.out.bedxpos
        .combine(
             GET_REFERENCE_HUMAN.out.fasta
             .map{ species, fasta, fai ->
                 tuple(fasta, fai)
             }
        )
        .set{ regionsxvar_ch }

    MAKE_REGIONSXVAR(regionsxvar_ch)

    MAKE_REGIONSXVAR.out.regionsxvar
        .combine(MAKE_VEP_CONFIG.out.config)
        .set{ vep_ch }
    
    ANNOTATE_VEP(vep_ch)
    PROCESS_VEP(ANNOTATE_VEP.out.vep)
    
    // Get contig counts
    GET_NUM_CONTIGS(GET_REFERENCE.out.fasta)
    
    // Liftover tables
    GET_REFERENCE.out.fasta
        .combine(
            GET_REFERENCE_HUMAN.out.fasta
            .map{ species, fasta, fai ->
                tuple(fasta, fai)
            }
         )
        .combine(MAKE_LIFTOVER_BED.out.lifted_beds, by: 0)
        .combine(MAKE_REGIONSXVAR.out.regionsxvar)
        .set{ sp_fastas_liftbed_vep_ch }
    
    MAKE_POSINFO(sp_fastas_liftbed_vep_ch)

    MAKE_POSINFO.out.pos_info
        .combine(PROCESS_VEP.out.vep_processed)
        .combine(GET_REFERENCE.out.fasta, by: 0)
        .set{ sp_posinfo_vep_fasta_ch }

    MAKE_LIFTOVER_TABLE(sp_posinfo_vep_fasta_ch)
    
    Channel.fromPath(params.species_file)
        .set{ sample_mapping_ch }

    MAKE_POSINFO.out.pos_info
        .combine(MAKE_LIFTOVER_TABLE.out.lifttable, by: 0)
        .mix(existing_lifttable_ch)
        .set { all_lifttables_ch }

    // Merge MTs
    VCF2MT.out.mt
        .mix(
            skip_formatting_ch
            .filter { species, mt_base_name, mt_file, vcf_file, alias_f, done -> done }
            .map { species, mt_base_name, mt_file, vcf_file, alias_f, done ->
                tuple(species, mt_file)
            }
        )
        .groupTuple()
        .combine(all_lifttables_ch, by:0) 
        .combine(GET_NUM_CONTIGS.out.contig_count, by: 0)
        .combine(sample_mapping_ch)
        .set{ sp_mts_lifttable_ch }
    
    MERGEMTS(sp_mts_lifttable_ch)

    MAKE_REGIONSXVAR.out.regionsxvar
        .combine(
            GET_REFERENCE_HUMAN.out.fasta
            .map{ species, fasta, fai ->
                tuple(fasta, fai)
            })
        .flatten()
        .concat(
            MERGEMTS.out.speciesht
            .map { it[1] }
            .collect()
        )
        .toList()
        .set{ vep_fasta_varhts_ch }
    
    MERGE_SPECIES_HTS(vep_fasta_varhts_ch)
}
