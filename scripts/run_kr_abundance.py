import os
import glob
import argparse
import pandas as pd
import sys
sys.path.append(os.getcwd())
from Metagenomics_pipeline.kraken_abundance_pipeline import process_sample, aggregate_kraken_results, generate_abundance_plots
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

def create_sample_id_df(input_dir):
    """
    Create a DataFrame with sample IDs based on the input FASTQ file names.
    """
    sample_ids = [os.path.basename(f).replace("_R1.fastq.gz", "").replace("_R1.fastq", "") for f in glob.glob(os.path.join(input_dir, "*_R1.fastq*"))]
    sample_id_df = pd.DataFrame(sample_ids, columns=["Sample_IDs"])
    return sample_id_df

def main():
    parser = argparse.ArgumentParser(description="Pipeline for Trimmomatic trimming, Bowtie2 host depletion (optional), and Kraken2 taxonomic classification.")
    parser.add_argument("--kraken_db", required=True, help="Path to Kraken2 database.")
    parser.add_argument("--bowtie2_index", help="Path to Bowtie2 index (optional).")
    parser.add_argument("--output_dir", required=True, help="Directory to save output files.")
    parser.add_argument("--input_dir", required=True, help="Directory containing input FASTQ files.")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads to use for Trimmomatic, Bowtie2, and Kraken2.")
    parser.add_argument("--metadata_file", help="Path to the metadata CSV file (optional).")
    parser.add_argument("--no_metadata", action='store_true', help="Use sample IDs as metadata instead of a metadata file.")
    parser.add_argument("--read_count", type=int, default=0, help="Minimum read count threshold.")
    parser.add_argument("--top_N", type=int, default=None, help="Select the top N most common viruses or bacteria.")
    parser.add_argument("--no_bowtie2", action='store_true', help="Skip Bowtie2 host depletion.")
    parser.add_argument("--bacteria", action='store_true', help="Generate bacterial abundance plots.")
    parser.add_argument("--virus", action='store_true', help="Generate viral abundance plots.")
    parser.add_argument("--use_precomputed_reports", action='store_true', help="Use precomputed Kraken reports instead of running Kraken2.")
    
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # Check Kraken database existence
    if not os.path.isdir(args.kraken_db):
        logging.error(f"Kraken database directory '{args.kraken_db}' not found.")
        sys.exit(1)


    run_bowtie = not args.no_bowtie2 and args.bowtie2_index is not None

    for forward in glob.glob(os.path.join(args.input_dir, "*_R1.fastq*")):
        base_name = os.path.basename(forward).replace("_R1.fastq.gz", "").replace("_R1.fastq", "")
 
        reverse = os.path.join(args.input_dir, f"{base_name}_R2.fastq.gz") if forward.endswith(".gz") else os.path.join(args.input_dir, f"{base_name}_R2.fastq")
        
        if not os.path.isfile(reverse):
            reverse = None

        process_sample(forward, reverse, base_name, args.bowtie2_index, args.kraken_db, args.output_dir, args.threads, run_bowtie,args.use_precomputed_reports)


    # Load metadata or create sample ID DataFrame
    if args.no_metadata:
        sample_id_df = create_sample_id_df(args.input_dir)
        logging.info("Using sample IDs as metadata.")
        sample_id_df.to_csv(os.path.join(args.output_dir, "sample_ids.csv"), index=False)
        merged_tsv_path = aggregate_kraken_results(args.output_dir, sample_id_df=sample_id_df, read_count=args.read_count)
    else:
        if not args.metadata_file:
            raise ValueError("Metadata file must be provided if --no_metadata is not specified.")
        elif not os.path.isfile(args.metadata_file):
            logging.error(f"Metadata file '{args.metadata_file}' not found.")
            sys.exit(1)
        merged_tsv_path = aggregate_kraken_results(args.output_dir, metadata_file=args.metadata_file, read_count=args.read_count)

    # Process each sample in the input directory
    for forward in glob.glob(os.path.join(args.input_dir, "*_R1.fastq*")):
        base_name = os.path.basename(forward).replace("_R1.fastq.gz", "").replace("_R1.fastq", "")
        reverse = os.path.join(args.input_dir, f"{base_name}_R2.fastq.gz") if forward.endswith(".gz") else os.path.join(args.input_dir, f"{base_name}_R2.fastq")
        
        if not os.path.isfile(reverse):
            logging.warning(f"Reverse file {reverse} not found. Skipping sample {base_name}.")
            continue

        try:
            process_sample(forward, reverse, base_name, args.bowtie2_index, args.kraken_db, args.output_dir, args.threads, run_bowtie, args.use_precomputed_reports)
        except Exception as e:
            logging.error(f"Error processing sample {base_name}: {e}")

    # Generate abundance plots based on provided flags
    if merged_tsv_path and os.path.isfile(merged_tsv_path):
        if args.virus:
            logging.info("Generating viral abundance plots.")
            generate_abundance_plots(merged_tsv_path, args.top_N)
        elif args.bacteria:
            logging.info("Generating bacterial abundance plots.")
            generate_abundance_plots(merged_tsv_path, args.top_N)
        else:
            logging.warning("No plot type specified. Use --virus or --bacteria to generate plots.")

if __name__ == "__main__":
    main()
