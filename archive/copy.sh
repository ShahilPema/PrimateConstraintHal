#!/bin/bash

# Define the destination directory
destination_directory="/home/ubuntu/PrimateConstraintHal/Cebus_albifrons_parquets"

# List of files to copy
files=(
    "/home/ubuntu/PrimateConstraintHal/work/73/ec49822788cdb75fcfe4ee7abc86c8/PD_0077.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/ba/2b243f77fb44cd8825edb5c59e5d28/PD_0007.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/f2/8e97ef48e638b668240e0aaa812f4f/PD_0141.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/df/33dc00a8698b37488209ca9a85be9d/PD_0013.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/9d/9ea6fa5544730442b0174686dac931/PD_0011.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/44/72a07d8637fe682d339596c239b388/PD_0079.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/f1/945a22af4f0484dcc664186b88e8a1/PD_0012.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/2a/7c6074710dbe75f7eb869cc55d0c38/PD_0374.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/71/b18c2daf1162c16a98e5b181319d4a/PD_0370.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/6b/95462a4af86606731a6e307f580fef/PD_0369.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/9b/7cbe5c5c6994a8ec71374cd07f46e4/PD_0373.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/67/d360b93758da2746902437ad55f128/PD_0361.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/d3/a1bc8d26f444d13ac9dfffc4047254/PD_0360.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/d0/8e486c7348543087bb53293b3446e6/PD_0358.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/7f/ea7449854f97a738037061bf308963/PD_0078.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/e5/79e982009d569533b9cb54f2d1e762/PD_0355.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/32/d8538e81d291d3d4cd1353500f0e09/PD_0311.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/92/909e3f3fffba2b422f32b6cf13f05d/PD_0310.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/8f/eb10c56fd76231c27cfdcba9ece236/PD_0375.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/7f/72104e0dfdedc438986eabe127d8f4/PD_0363.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/51/b39af82ec53c2f84975244cc00046b/PD_0402.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/90/a0381f65eaae585f0b275ef5e7f8f9/PD_0368.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/8e/8b8f8945e7c0906948cd34363e0e90/PD_0372.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/83/5f601787870918463bb5ab8faa5a27/PD_0354.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/b5/ae4e2a63f0967431230170a13e2418/PD_0357.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/8d/019e0237d8953a07ff982eb33a90a4/PD_0356.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/5d/780a1cb4259b945864b4265d364769/PD_0371.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/0f/533427af959b61359c9bda6fa073c5/PD_0367.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/8e/401e22e50f9e1258daa214bdea360f/PD_0359.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/77/fceaa62b68035556b34fc76d0dbe3a/PD_0353.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/1d/ef911f0147426773868823b523e7a3/PD_0362.parquet"
    "/home/ubuntu/PrimateConstraintHal/work/78/91d37193b3dd62aa8d0636248c785d/PD_0366.parquet"
)

# Export the destination directory
export destination_directory

# Calculate the number of files
num_files=${#files[@]}

# Use xargs to copy files in parallel
echo "${files[@]}" | tr ' ' '\n' | xargs -I {} -P "$num_files" cp {} "$destination_directory"

echo "All files copied to $destination_directory"
