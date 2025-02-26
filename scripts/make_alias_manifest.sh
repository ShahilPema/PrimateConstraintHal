#!/bin/bash

# Define the output file
output_file="alias_manifest.txt"

# Create the file with the specified content
cat <<EOL > $output_file
Aotus_nancymaae	GCA_000952055.2
Ateles_geoffroyi_b	GCA_023783555.1
Carlito_syrichta	GCA_000164805.1
Cebus_albifrons	GCA_023783575.1
Cercocebus_atys	GCA_000955945.1
Cercopithecus_albogularis	GCA_023783535.1
Chlorocebus_aethiops	GCA_023783515.1
Colobus_guereza	GCA_021498455.1
Daubentonia_madagascariensis	GCA_023783475.1
Erythrocebus_patas	GCA_023783455.1
Galago_moholi	GCA_023783435.1
Gorilla_gorilla	GCA_900006655.3
Loris_tardigradus	GCA_023783135.1
Macaca_mulatta	GCA_003339765.3
Mandrillus_sphinx	GCA_023783085.1
Microcebus_murinus	GCA_000165445.3
Nomascus_siki_b	GCA_023783065.1
Otolemur_garnettii	GCA_000181295.3
Pan_troglodytes	GCA_002880755.3
Papio_anubis	GCA_000264685.2
Pithecia_pithecia	GCA_023779675.1
Pongo_abelii	GCA_002880775.3
Pongo_pygmaeus	GCA_023767775.1
Propithecus_coquerelli	GCA_000956105.1
Rhinopithecus_roxellana	GCA_007565055.1
Saguinus_midas	GCA_021498475.1
Sapajus_apella	GCA_023762875.2
Theropithecus_gelada	GCA_003255815.1
EOL

# Inform the user
echo "File '$output_file' has been created."

