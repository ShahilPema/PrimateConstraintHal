#!/usr/bin/env nextflow

process getSystemCpus {
    cache false
    output:
    env(cpus), emit:cpus

    script:
    """
    cpus=\$(nproc)
    """
}

process PREP_FILES_BED {
    executor 'local'
    //clusterOptions '--cpus-per-task=10 --mem=20G --time=01:00:00'

    output:
    path 'MANE_transcript.txt', emit: mane_transcript
    path 'Homo_sapiens.GRCh38.112.chr.gff3', emit: gff_file

    script:
    """
    bash ${projectDir}/scripts/prep_files_bed.sh
    """
}

process MAKE_BED {
    executor 'local'
    //clusterOptions '--cpus-per-task=1 --mem=4G --time=01:00:00'

    input:
    path mane_transcript
    path gff_file

    output:
    path 'features_by_loci.bed', emit: bed_output

    script:
    """
    python ${projectDir}/scripts/make_bed.py \
        --mane_transcript $mane_transcript \
        --gff_file $gff_file \
        --output features_by_loci.bed
    """
}

process CHECK_BASESPACE_AUTH {
    executor 'local'
    errorStrategy 'terminate'
    containerOptions "--bind $HOME:$HOME"

    output:
    path 'auth_check.txt', emit: auth_status

    script:
    """
    if [ -d "\$HOME/.basespace" ]; then
        echo "BaseSpace authentication found" > auth_check.txt
    else
        echo "Error: BaseSpace authentication not found." >&2
        echo "Please install Illumina BaseSpace CLI and run 'bs auth' to authenticate before proceeding." >&2
        exit 1
    fi
    """
}

process GET_VCF_IDS {
    executor 'local'

    containerOptions "--bind $HOME:$HOME"

    input:

    output:
    path "vcf_ids.txt", emit: vcf_ids_file

    script:
    """
    bash ${projectDir}/scripts/get_vcf_ids.sh
    """
}

process DOWNLOAD_VCFS {
    containerOptions "--bind $HOME:$HOME"

    //executor 'slurm'
    //clusterOptions "--cpus-per-task=1 --mem=1G --time=1:00:00 --output=local_logs/%j.out --error=local_logs/%j.err"
    
    executor 'local'
    cpus 1
    memory '1 GB'
    maxForks 10
    
    input:
    tuple val(id), val(file)

    output:
    path "*.vcf.gz", emit: vcf_file

    script:
    """
    data_dir=${projectDir}/vcfs
    mkdir -p "\${data_dir}"
    if [ ! -f "\${data_dir}/${file}" ]; then
        bash ${projectDir}/scripts/download_vcf.sh $id "\${data_dir}"
    fi
    ln -s "\${data_dir}/${file}" "${file}"
    """
}

process MAKE_ALIAS_MANIFEST {
    executor 'local'

    input:

    output:
    path "alias_manifest.txt", emit: alias_manifest

    script:
    """
    bash ${projectDir}/scripts/make_alias_manifest.sh
    """
}

process GET_ALIAS_FILES {
    executor 'local'

    input:
    tuple val(species), val(gca_id)

    output:
    tuple val(species), path("*.txt"), emit: alias_file

    script:
    """
    bash ${projectDir}/scripts/get_alias_file.sh ${species} ${gca_id} .
    """
}

process MAKE_SAMPLE2REF_MAP {
    executor 'local'

    output:
    path "sample2ref.txt", emit: sample2ref_file

    script:
    """
    bash ${projectDir}/scripts/make_sample2ref_map.sh
    """
}

process MAKE_POS_BED {
    executor 'local'

    input:
    path(bed_file)

    output:
    path("features_by_position.parquet"), emit: bedxpos

    script:
    """
    python ${projectDir}/scripts/make_pos_bed.py ${bed_file} .
    """

}

process MAKE_HG38_BEDS {
    executor 'local'
    
    input:
    path(bedxpos)
    
    output:
    path("hg38_*_part_*.bed"), emit:hg38beds

    script:
    """
    python ${projectDir}/scripts/make_hg38_beds.py $bedxpos hg38
    """
}

process MAKE_LIFTOVER_BED {
    executor 'local'
    maxForks 3

    input:
    tuple val(species), path(hal_file), path(bed_files)

    output:
    tuple val(species), path("${species}_lifted.bed"), emit: lifttable

    script:
    """
    data_dir=${projectDir}/lift_tables/${species}
    mkdir -p "\$data_dir"
    if [ ! -f "\${data_dir}/${species}_lifted.bed" ]; then
        tmp_missing_files=\$(ls ${bed_files} | while read file; do
            base_name=\$(basename "\$file" .bed)
            if [ ! -f "\${data_dir}/\${base_name}_final.bed" ]; then
                echo "\$file"
            fi
        done)
        
        missing_files=\$(echo "\$tmp_missing_files" | tr ' ' '\\n')
    
        if [ ! -z "\$missing_files" ]; then
            echo "\$missing_files" | parallel -j \$((${task.cpus} - 40)) bash ${projectDir}/scripts/liftover_bed.sh {} $species $hal_file -o \${data_dir}
        fi
    
        cat "\${data_dir}"/*_final.bed > "\${data_dir}/${species}_lifted.bed"
    fi
    ln -sf "\${data_dir}/${species}_lifted.bed" "${species}_lifted.bed"
    """
}

process GET_REFERENCE {
    executor 'local'
    maxForks 3

    input:
    tuple val(species), path(hal)

    output:
    tuple val(species), path("*.fa"), path("*.fa.fai"), emit: fasta

    script:
    """
    bash ${projectDir}/scripts/getreference.sh $species $hal    
    """
}

process GET_REFERENCE_HUMAN {
    executor 'local'
    maxForks 3

    input:
    path(hal)

    output:
    tuple val('Homo_sapiens'), path("*.fa"), path("*.fa.fai"), emit: fasta

    script:
    """
    bash ${projectDir}/scripts/getreference.sh Homo_sapiens $hal
    """
}

process MAKE_VEP_CONFIG {
    
    output:
    path("vep_config.json"), emit:config

    script:
    """
    ${projectDir}/scripts/vep_config.sh 
    """
}

process ANNOTATE_VEP {
    maxForks 1
    container = "${projectDir}/annotation.sif"

    input:
    tuple val(cpus), path(position_table_path), path(human_fasta_path), path(human_index_path), path(vep_config)

    output:
    path("vep.ht"), emit: vep

    script:
    """
    data_dir=${projectDir}/vep
    mkdir -p "\${data_dir}"
    if [ ! -d "\${data_dir}/vep.ht" ]; then
       python3 ${projectDir}/scripts/vep.py --cpus "${cpus}" \
                                    --position_table "${position_table_path}" \
                                    --human_fasta_path "${human_fasta_path}" \
                                    --human_index_path "${human_index_path}" \
                                    --vep_config "${vep_config}" \
                                    --output "\${data_dir}"
    fi
    ln -sf "\${data_dir}/vep.ht" "vep.ht"
    """
}

process PROCESS_VEP {  
    memory '115 GB'   

    input:
    path(vep_ht)

    output:
    path("vep_cdsinfo.ht"), emit: vep_processed

    script:
    """
    data_dir=${projectDir}/vep
    mkdir -p "\${data_dir}"
    if [ ! -d "\${data_dir}/vep_cdsinfo.ht" ]; then
       python3 ${projectDir}/scripts/processvep.py --vep_ht "${vep_ht}" --output "\${data_dir}"
    fi
    ln -sf "\${data_dir}/vep_cdsinfo.ht" "vep_cdsinfo.ht"
    """
}

process MAKE_LIFTOVER_TABLE {
    executor 'local'
    maxForks 1

    input:
    tuple val(species), path(species_fasta_path), path(species_index_path), path(human_fasta_path), path(human_index_path), path(species_bed_path),  path(vep_path), val(all_cpus)

    output:
    tuple val(species), path("${species}_pos_info.ht"), path("${species}_var_info.ht"), emit: lifttable

    script:
    """ 
    data_dir=${projectDir}/lift_tables/${species}
    mkdir -p "\$data_dir"
    if [ ! -d "\${data_dir}/${species}_pos_info.ht" ] || [ ! -d "\${data_dir}/${species}_var_info.ht" ]; then    
        python3 ${projectDir}/scripts/bed2liftover_ht.py \
            --species "$species" \
            --human_fasta_path "$human_fasta_path" \
            --human_index_path "$human_index_path" \
            --species_fasta_path "$species_fasta_path" \
            --species_index_path "$species_index_path" \
            --species_bed_path "$species_bed_path" \
            --vep_annotations "$vep_path" \
            --all_cpus "${all_cpus}" \
            --output "\${data_dir}" 
    fi  
    ln -sf "\${data_dir}/${species}_pos_info.ht" "${species}_pos_info.ht"
    ln -sf "\${data_dir}/${species}_var_info.ht" "${species}_var_info.ht"
    """ 
}


process GET_CHROM_FORMAT_VCF {
    executor 'local'
    maxForks 3

    input:
    tuple val(species), path(vcf), path(alias)

    output:
    tuple val(species), path(vcf), path(alias), env(format), emit: format
    
    script:
    """
    format=\$(python ${projectDir}/scripts/getformat.py $vcf $alias)
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
    format=\$(python ${projectDir}/scripts/getformat_hal.py $species $alias $hal_)
    """
}

process VCF2MT {
    maxForks 12

    input:
    tuple val(species), path(vcf_file), val(vcf_format), path(alias_file), val(hal_format), path(primate_fasta), path(primate_index), val(cpus)

    output:
    tuple val(species), path("*.mt"), emit: mt

    script:
    vcf_bname=vcf_file.baseName.replaceAll(/\.SNV\.vcf$/, '')
    """
    file="${vcf_bname}.mt"
    data_dir=${projectDir}/mts
    mkdir -p "\${data_dir}"
    if [ ! -d "\${data_dir}/\${file}" ]; then
       python3 ${projectDir}/scripts/vcf2mt.py --species "${species}" \
                                       --vcf_file "${vcf_file}" \
                                       --alias_file "${alias_file}" \
                                       --vcf_format "${vcf_format}" \
                                       --hal_format "${hal_format}" \
                                       --primate_fasta_path "${primate_fasta}" \
                                       --primate_index_path "${primate_index}" \
                                       --cpus "${cpus}" \
                                       --output "\${data_dir}"
    fi
    ln -s "\${data_dir}/\${file}" "\${file}"
    """
}

process GET_NUM_CONTIGS {
    executor 'local'
    maxForks 1

    input:
    tuple val(species), path(fasta), path(indexes)

    output:
    tuple val(species), env(contigs), emit: contig_count

    script:
    """
    contigs=\$(python ${projectDir}/scripts/get_contig_count.py --species $species --species_fasta_path $fasta --species_index_path $indexes)
    """
}

process MERGEMTS {
    maxForks 1

    input:
    tuple val(species), path(mts), path(lifttable), val(contig_count), val(all_cpus)

    output:
    tuple val(species), path("${species}_lifted_vars.ht"), path("${species}_colinfo_vars.ht"), emit: speciesht

    script:
    """
    data_dir=${projectDir}/lift_tables/${species}
    mkdir -p "\$data_dir"
    if [ ! -d "\${data_dir}/${species}_lifted_vars.ht" ] || [ ! -d "\${data_dir}/${species}_colinfo_vars.ht" ]; then
        python3 ${projectDir}/scripts/liftvars.py \
            --species "$species" \
            --sample_mts "$mts" \
            --full_ht "$lifttable" \
            --contig_count "$contig_count" \
            --all_cpus "${all_cpus}" \
            --outdir "\${data_dir}"
    fi
    ln -sf "\${data_dir}/${species}_lifted_vars.ht" "${species}_lifted_vars.ht"
    ln -sf "\${data_dir}/${species}_colinfo_vars.ht" "${species}_colinfo_vars.ht"
    """
}

process MERGE_SPECIES_HTS {
    maxForks 1

    input:
    tuple path(vep), path(mapping), path(human_fa), path(human_faidx), path(var_hts), path(pos_hts)

    output:
    path("aggregated_primate_counts_with_position_info.ht"), emit: all_species

    script:
    """
    python3 ${projectDir}/scripts/merge_all.py \
         --vep "$vep" \
         --mapping "$mapping" \
         --human_fa "$human_fa" \
         --human_faidx "$human_faidx" \
         --var_hts "$var_hts" \
         --pos_hts "$pos_hts" \
    """
}

workflow {
    getSystemCpus()
    PREP_FILES_BED()
    MAKE_BED(PREP_FILES_BED.out.mane_transcript, PREP_FILES_BED.out.gff_file)
    CHECK_BASESPACE_AUTH()
    GET_VCF_IDS()
    vcf_ids_ch = GET_VCF_IDS.out.vcf_ids_file.splitText().map { line ->
        def fields = line.trim().split('\t')
        tuple(fields[0], fields[1])
    }
    DOWNLOAD_VCFS(vcf_ids_ch)
    MAKE_ALIAS_MANIFEST()
    alias_manifest_ch = MAKE_ALIAS_MANIFEST.out.alias_manifest
         .splitCsv(sep: '\t')
         .map { line -> tuple(line[0].trim(), line[1].trim()) }
    GET_ALIAS_FILES(alias_manifest_ch)
    MAKE_SAMPLE2REF_MAP()
    MAKE_POS_BED(MAKE_BED.out.bed_output)
    Channel.fromPath('447-mammalian-2022v1.hal')
    .set { hal }
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
        def id = vcf_file.getBaseName().tokenize( '.' )[0]  // Extract the sample ID from the filename (before the first dot)
        def ref_species = species2ref[id]  // Lookup the reference species using the ID
        return tuple(ref_species, vcf_file)
    }
    .combine(GET_ALIAS_FILES.out.alias_file, by: 0)
    .set { sp_vcf_alias_ch }
    MAKE_HG38_BEDS(MAKE_POS_BED.out.bedxpos)
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
    MAKE_LIFTOVER_BED(sp_hal_hg38beds_ch)
    GET_REFERENCE(sp_hal_ch)
    GET_REFERENCE_HUMAN(hal)
    GET_CHROM_FORMAT_VCF(sp_vcf_alias_ch) 
    GET_ALIAS_FILES.out.alias_file
    .combine(hal)
    .set{ sp_alias_hal_ch }
    GET_CHROM_FORMAT_HAL(sp_alias_hal_ch)
    GET_CHROM_FORMAT_VCF.out.format
    .filter { it[3] != 'Other' }
    .map { species, vcf, chromAlias, format -> 
        tuple(species, vcf, format) //remove after modifying
    }
    .combine(GET_ALIAS_FILES.out.alias_file, by: 0)
    .combine(GET_CHROM_FORMAT_HAL.out.format, by: 0)
    .combine(GET_REFERENCE.out.fasta, by: 0)
    .combine(getSystemCpus.out.cpus)
    .set{ sp_vcf_format_alias_halformat_ref_ch }
    //If you modify filtering in VCF2MT, you will need to rework this script to remove those site from the AN count
    //primad vcfs already have variants filtered according to the criteria in this script
    VCF2MT(sp_vcf_format_alias_halformat_ref_ch)
    MAKE_VEP_CONFIG()
    getSystemCpus.out.cpus
    .combine(MAKE_POS_BED.out.bedxpos)
    .combine(
         GET_REFERENCE_HUMAN.out.fasta
         .map{ species, fasta, fai ->
             tuple(fasta, fai)
         }
    )
    .combine(MAKE_VEP_CONFIG.out.config)
    .set{ vep_ch }
    ANNOTATE_VEP(vep_ch)
    PROCESS_VEP(ANNOTATE_VEP.out.vep)
    GET_REFERENCE.out.fasta
    .combine(
        GET_REFERENCE_HUMAN.out.fasta
        .map{ species, fasta, fai ->
            tuple(fasta, fai)
        }
     )
    .combine(MAKE_LIFTOVER_BED.out.lifttable, by: 0)
    .combine(PROCESS_VEP.out.vep_processed)
    .combine(getSystemCpus.out.cpus)
    .set{ sp_fastas_liftbed_vep_ch }
    MAKE_LIFTOVER_TABLE(sp_fastas_liftbed_vep_ch)
    GET_NUM_CONTIGS(GET_REFERENCE.out.fasta)
    VCF2MT.out.mt
    .groupTuple()
    .combine(
        MAKE_LIFTOVER_TABLE.out.lifttable
        .map{ species, pos_table, var_table ->
            tuple(species, var_table)
        },
        by: 0) 
    .combine(GET_NUM_CONTIGS.out.contig_count, by: 0)
    .combine(getSystemCpus.out.cpus)
    .set{ sp_mts_lifttable_ch }
    MERGEMTS(sp_mts_lifttable_ch)
    Channel.fromPath('233primates_s1.csv')
    .set{ sample_mapping_ch }
    ANNOTATE_VEP.out.vep
    .combine(sample_mapping_ch)
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
    .concat(
        MAKE_LIFTOVER_TABLE.out.lifttable
        .map{ species, pos_table, var_table ->
            tuple(pos_table)
        }
        .collect()
    )
    .toList()
    .set{ vep_varhts_poshts_ch }
    MERGE_SPECIES_HTS(vep_varhts_poshts_ch)
}
