#!/bin/bash

# File name
output_file="ref_species.txt"

# Array of species
species=(
    "Ateles_geoffroyi_b"
    "Carlito_syrichta"
    "Cebus_albifrons"
    "Cercocebus_atys"
    "Cercopithecus_albogularis"
    "Chlorocebus_aethiops"
    "Colobus_guereza"
    "Daubentonia_madagascariensis"
    "Erythrocebus_patas"
    "Galago_moholi"
    "Gorilla_gorilla"
    "Loris_tardigradus"
    "Macaca_mulatta"
    "Mandrillus_sphinx"
    "Microcebus_murinus"
    "Nomascus_siki_b"
    "Otolemur_garnettii"
    "Pan_troglodytes"
    "Papio_anubis"
    "Pithecia_pithecia"
    "Pongo_abelii"
    "Pongo_pygmaeus"
    "Propithecus_coquerelli"
    "Rhinopithecus_roxellana"
    "Saguinus_midas"
    "Sapajus_apella"
    "Theropithecus_gelada"
    "hg38"
)

# Write species to file
printf "%s\n" "${species[@]}" > "$output_file"

