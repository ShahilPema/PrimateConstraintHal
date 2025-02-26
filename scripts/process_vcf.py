import subprocess
import os
import polars as pl
import sys

def main(species, alias_file, vcf_file, hal_file, outdir):
    try:
        # Step 1: Get Contigs from HAL File
        hal_command = f"halStats --sequences {species} {hal_file}"
        result = subprocess.run(hal_command, shell=True, capture_output=True, text=True, check=True)
        hal_contigs = result.stdout.split(',')
        hal_data = pl.DataFrame(hal_contigs, schema=pl.Schema({'column_1': pl.Utf8}))
    except subprocess.CalledProcessError as e:
        print(f"Error executing halStats command: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error occurred while getting contigs from HAL file: {e}")
        sys.exit(1)

    try:
        # Step 2: Create Alias Dictionary
        alias_data = pl.read_csv(alias_file, separator='\t')
        alias_data.columns = [item.replace(" ", "").replace("#", "") for item in alias_data.columns]
        alias_dict = {}
        for column in alias_data.columns:
            format = alias_data.select(column)
            format = format.with_columns(
                pl.lit(column).alias('format')
            )
            alias_dict = {**alias_dict, **dict(format.iter_rows())}
    except FileNotFoundError:
        print(f"Alias file not found: {alias_file}")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error occurred while creating alias dictionary: {e}")
        sys.exit(1)

    try:
        # Step 3: Determine Assemblies for HAL and VCF
        vcf_data = pl.read_csv(vcf_file, separator='\t', comment_prefix='#', has_header=False)
        vcf_assembly = vcf_data.select(pl.col('column_1')).to_series().replace_strict(alias_dict, default=None).value_counts().sort('count', descending=True).select('column_1').to_series().item(0)
        hal_assembly = hal_data.select(pl.col('column_1')).to_series().replace_strict(alias_dict, default=None).value_counts().sort('count', descending=True).select('column_1').to_series().item(0)
        print(f"vcf_assembly:{vcf_assembly}\nhal_assembly:{hal_assembly}")
    except FileNotFoundError:
        print(f"VCF file not found: {vcf_file}")
        sys.exit(1)
    except Exception as e:
        print(f"Error determining assemblies for HAL and VCF: {e}")
        sys.exit(1)

    class HALAssemblyNotFoundError(Exception):
        pass

    class VCFAssemblyNotFoundError(Exception):
        pass
    
    try:
        # Step 4: Prepare VCF and Remap Contigs if Necessary
        if hal_assembly == None:
            raise HALAssemblyNotFoundError("HAL assembly contigs not found in alias file.")
        if vcf_assembly == None:
            with open(f'{outdir}/mismatch_counts.parquet', 'w') as f:
                pass
            raise VCFAssemblyNotFoundError("VCF file was called with reference different from HAL assembly.")
        if hal_assembly == vcf_assembly:
            print("Contigs in HAL and VCF match. No remapping necessary.")
        else:
            print("Remapping contigs...")
            mapping_data = alias_data.select(vcf_assembly, hal_assembly)
            mapping_dict = dict(mapping_data.iter_rows())

            vcf_data = vcf_data.with_columns(
                pl.col('column_1').replace(mapping_dict).alias('column_1')
            )
        vcf_data = vcf_data.with_columns(
            pl.col('column_10').alias('column_9b'),
            (pl.col('column_1').cast(pl.Utf8) + pl.lit(':') + pl.col('column_2').cast(pl.Utf8)).alias('column_9c')
        ).drop('column_10')
        vcf_data = vcf_data.select(sorted(vcf_data.columns))
        vcf_data = vcf_data.filter(
            (pl.col('column_4').str.len_chars() == 1) & (pl.col('column_5').str.len_chars() == 1)
        )
    except HALAssemblyNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except VCFAssemblyNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(0)    
    except KeyError as e:
        print(f"Error during remapping: missing key in alias dictionary - {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error during remapping contigs: {e}")
        sys.exit(1)

    try:
        # Step 5: Write BED File
        vcf_basename = os.path.basename(vcf_file).split('.')[0]
        output_parquet = (f"{outdir}/{vcf_basename}.parquet")
        column_dict = {
            "column_1": f"{species}_contig",
            "column_2": f"{species}_position",
            "column_3": "ID",
            "column_4": "REF",
            "column_5": "ALT",
            "column_6": "QUAL",
            "column_7": "FILTER",
            "column_8": "INFO",
            "column_9": "FORMAT",
            "column_9b": "SAMPLE",
            "column_9c": "POS_ORIGIN"
        }
        vcf_data = vcf_data.rename(column_dict)
        vcf_data.write_parquet(output_parquet)
    except Exception as e:
        print(f"Error writing BED file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    try:
        species = sys.argv[1]
        vcf_file = sys.argv[2]
        alias_file = sys.argv[3]
        hal_file = sys.argv[4]
        outdir = sys.argv[5]
        
        main(species, alias_file, vcf_file, hal_file, outdir)
    except IndexError:
        print("Usage: python script.py <species> <alias_file> <vcf_file> <hal_file> <output>")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error in main execution: {e}")
        sys.exit(1)
