#!/usr/bin/env nextflow

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

process VCF2PARQ {
    //executor 'slurm'
    //clusterOptions "--partition=short --time=04:00:00 --cpus-per-task=4 --mem=10G --output=slurm_logs/%j.out --error=slurm_logs/%j.err"
    
    executor 'local'
    memory '40 GB'
    
    errorStrategy 'finish'
 
    input:
    tuple val(species), path(vcf_file), path(alias_file), path(hal_file)

    output:
    tuple val(species), path("*.parquet"), emit: sample_parq_file

    script:
    """
    data_dir=${projectDir}/parquets
    file=\$(basename ${vcf_file} .SNV.vcf.gz)
    mkdir -p "\${data_dir}"
    
    # Check if this is a known mismatch_counts case
    if grep -q "\${file}" "\${data_dir}/mismatch_counts_files.txt" 2>/dev/null; then
        touch "mismatch_counts.parquet"
        exit 0
    fi
    
    # Check if regular parquet file already exists
    if [ -f "\${data_dir}/\${file}.parquet" ]; then
        ln -s "\${data_dir}/\${file}.parquet" "\${file}.parquet"
        exit 0
    fi
    
    # Continue with normal processing if file doesn't exist
    python ${projectDir}/scripts/process_vcf.py $species $vcf_file $alias_file $hal_file "."
    
    # Check output and handle accordingly
    if [ -f "mismatch_counts.parquet" ]; then
        echo "\${file}" >> "\${data_dir}/mismatch_counts_files.txt"
    else
        mv "\${file}.parquet" "\${data_dir}/"
        ln -s "\${data_dir}/\${file}.parquet" "\${file}.parquet"
    fi
    """
}

process FILTER_VARIANTS {
    maxForks 60

    input:
    tuple val(species), path(parquet_file)

    output:
    tuple val(species), path("*_filtered.parquet"), emit: sample_parq_filt_file

    script:
    """
    python ${projectDir}/scripts/filter_variants.py $species $parquet_file
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

process MAKE_LIFTOVER_TABLE {
    executor 'local'
    //maxForks 1
    //cpus 180
    maxForks 3

    input:
    tuple val(species), path(hal_file), path(bed_files)

    output:
    tuple val(species), path("${species}_alignment.parquet"), emit: lifttable

    script:
    """
    data_dir=${projectDir}/lift_tables/${species}
    mkdir -p "\$data_dir"
    tmp_missing_files=\$(ls ${bed_files} | while read file; do
        base_name=\$(basename "\$file" .bed)
        if [ ! -f "\${data_dir}/\${base_name}_lifted_final.bed" ]; then
            echo "\$file"
        fi
    done)
    
    missing_files=\$(echo "\$tmp_missing_files" | tr ' ' '\\n')

    if [ ! -z "\$missing_files" ]; then
        echo "\$missing_files" | parallel -j \$((${task.cpus} - 40)) bash ${projectDir}/scripts/liftover_bed.sh {} $species $hal_file -o \${data_dir}
    fi

    python ${projectDir}/scripts/merge_lifted_beds.py $species \${data_dir}/*.bed

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

process ALIGNMENT_BY_VARS {
    executor 'local'
    maxForks 15

    input:
    tuple val(species), path(alignment_file), path(ref_vars), path(fasta), path(index)

    output:
    tuple val(species), path("${species}_alignment_vars.parquet"), emit: lifttable2

    script:
    """
    basename=\$(basename ${alignment_file} .parquet)
    file="\${basename}_vars.parquet"
    data_dir=${projectDir}/lifttablesvars
    mkdir -p "\${data_dir}"
    if [ ! -f "\${data_dir}/\${file}" ]; then
        python ${projectDir}/scripts/alignment_by_vars.py $species $alignment_file $ref_vars $fasta \${data_dir}/\${file}
    fi
    ln -s "\${data_dir}/\${file}" "\${file}"
    """

}

process GET_SAMPLE_COUNTS {
    //executor 'slurm'
    //clusterOptions "--partition=short --time=01:00:00 --cpus-per-task=4 --mem=20G --output=slurm_logs/%j.out --error=slurm_logs/%j.err"
    executor 'local'
    cpus 16
    memory '20GB'

    input:
    tuple val(species), path(lifttable), path(sample_parquet)

    output:
    tuple val(species), path("*_counts.parquet"), emit: sample_counts

    script:
    """
    python ${projectDir}/scripts/getsamplecounts.py $species $lifttable $sample_parquet .
    """
}

process MERGE_COUNTS_X_SPECIES {
    //executor 'slurm'
    //clusterOptions "--partition=short --time=01:00:00 --cpus-per-task=4 --mem=20G --output=slurm_logs/%j.out --error=slurm_logs/%j.err"
    executor 'local'
    maxForks 3

    input:
    tuple val(species), path(sample_parquets), path(bedxpos)

    output:
    tuple val(species), path("${species}_counts.parquet"), emit: species_counts

    script:
    """
    python ${projectDir}/scripts/mergecountsxspecies.py $species ${sample_parquets.join(',')} $bedxpos .
    """
}


process MERGE_ALL_SPECIES_COUNTS {
    //executor 'slurm'
    //clusterOptions "--partition=short --time=01:00:00 --cpus-per-task=5 --mem=10G --output=slurm_logs/%j.out --error=slurm_logs/%j.err"
    executor 'local'

    input:
    path(species_counts)

    output:
    path("merged_counts_*.parquet"), emit: merged_species_counts

    script:
    """
    python ${projectDir}/scripts/mergeallspecies.py ${species_counts.join(',')}
    """
}

workflow {
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
    .combine(hal)
    .set { sp_vcf_alias_hal_ch }
    VCF2PARQ(sp_vcf_alias_hal_ch)
    VCF2PARQ.out.sample_parq_file.filter { tuple ->
        def (species, file) = tuple
        !file.name.endsWith('mismatch_counts.parquet')
    }
    .set { variant_parq_ch }
    FILTER_VARIANTS(variant_parq_ch)
    MAKE_HG38_BEDS(MAKE_POS_BED.out.bedxpos)
    sp_vcf_alias_hal_ch
    .map { tuple -> tuple[0] }
    .unique()
    .combine(hal)
    .combine(MAKE_HG38_BEDS.out.hg38beds.toList())
    .set { sp_hg38beds_ch }
    MAKE_LIFTOVER_TABLE(sp_hg38beds_ch)
    sp_vcf_alias_hal_ch
    .map { species, vcf, alias, hal ->
        return tuple(species, hal)
    }
    .unique()
    .set { sp_hal_ch }
    GET_REFERENCE(sp_hal_ch)
    mil127_vars_ch = Channel.fromPath('127vars.parquet')
    FILTER_VARIANTS.out.sample_parq_filt_file
    .combine(MAKE_LIFTOVER_TABLE.out.lifttable, by:0)
    .map { species, sample_file, alignment ->
        return tuple(species, alignment, sample_file)
    }
    .filter { species, alignment, sample_file ->
        species == 'Macaca_mulatta'
    }
    .set { sp_alignment_call_ch }
    sp_alignment_call_ch.view()
    /*
    GET_SAMPLE_COUNTS(sp_alignment_call_ch)
    mil127_vars_ch
    .concat(
        GET_SAMPLE_COUNTS.out.sample_counts
        .map { it[1] }
    )
    .toList()
    .set { sample_counts_ch }
    sample_counts_ch.view()
    //MERGE_ALL_SPECIES_COUNTS(sample_counts_ch)
    */
}
