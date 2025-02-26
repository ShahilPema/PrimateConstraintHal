#!/bin/bash

# Define the file name
output_file="vep_config.json"

# Write the content to the file
cat << 'EOF' > $output_file
{
  "command": [
    "/opt/vep/src/ensembl-vep/vep",
    "--format", "vcf",
    "--json",
    "--no_stats",
    "--cache", "--offline",
    "--assembly", "GRCh38",
    "--dir", "/opt/vep/.vep",
    "--species", "homo_sapiens",
    "--plugin", "TSSDistance",
    "--plugin", "SpliceRegion",
    "--plugin", "LoF,loftee_path:/plugins,gerp_bigwig:/opt/vep/.vep/gerp_conservation_scores.homo_sapiens.GRCh38.bw,human_ancestor_fa:/opt/vep/.vep/human_ancestor.fa.gz,conservation_file:/opt/vep/.vep/loftee.sql",
    "-o", "STDOUT"
  ],
  "env": {
  },
  "vep_json_schema":
    "Struct{input:String,id:String,end:Int32,transcript_consequences:Array[Struct{allele_num:Int32,cdna_end:Int32,cdna_start:Int32,gene_id:String,protein_end:Int32,protein_start:Int32,variant_allele:String,codons:String,consequence_terms:Array[String],strand:Int32,transcript_id:String,impact:String,cds_start:Int32,cds_end:Int32,amino_acids:String,flags:Array[String],distance:Int32,lof:String,lof_flags:String,lof_filter:String,lof_info:String,spliceregion:Array[String],tssdistance:Int32}],intergenic_consequences:Array[Struct{allele_num:Int32,variant_allele:String,consequence_terms:Array[String],impact:String}],motif_feature_consequences:Array[Struct{allele_num:Int32,cdna_end:Int32,cdna_start:Int32,gene_id:String,protein_end:Int32,protein_start:Int32,variant_allele:String,codons:String,consequence_terms:Array[String],strand:Int32,transcript_id:String,impact:String,cds_start:Int32,cds_end:Int32,amino_acids:String,flags:Array[String]}],regulatory_feature_consequences:Array[Struct{allele_num:Int32,cdna_end:Int32,cdna_start:Int32,gene_id:String,protein_end:Int32,protein_start:Int32,variant_allele:String,codons:String,consequence_terms:Array[String],strand:Int32,transcript_id:String,impact:String,cds_start:Int32,cds_end:Int32,amino_acids:String,flags:Array[String],biotype:String,regulatory_feature_id:String}],assembly_name:String,seq_region_name:String,most_severe_consequence:String,start:Int32,allele_string:String,strand:Int32}"
}
EOF

echo "JSON configuration file '$output_file' has been created."

