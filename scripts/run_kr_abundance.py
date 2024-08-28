#import sys
#sys.path.append('/home/harouna/ARSNACAdata/bamfiles/mypipeline/Metagenomics_pipeline/')
import sys
import os
sys.path.append(os.getcwd())
from Metagenomics_pipeline.kraken_abundance_pipeline import process_sample, aggregate_kraken_results, generate_abundance_plots
import argparse
import glob
def main():
    parser = argparse.ArgumentParser(description="Pipeline for Trimmomatic trimming, Bowtie2 host depletion, and Kraken2 taxonomic classification.")
    parser.add_argument("--kraken_db", required=True, help="Path to Kraken2 database.")
    parser.add_argument("--bowtie2_index", required=True, help="Path to Bowtie2 index.")
    parser.add_argument("--output_dir", required=True, help="Directory to save output files.")
    parser.add_argument("--input_dir", required=True, help="Directory containing input FASTQ files.")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads to use for Trimmomatic, Bowtie2, and Kraken2.")
    parser.add_argument("--metadata_file", required=True, help="Path to the metadata CSV file.")
    parser.add_argument("--read_count", type=int, default=0, help="Minimum read count threshold.")
    parser.add_argument("--top_N", type=int, default=None, help="Select the top N most common viruses or bacteria.")
    parser.add_argument("--virus", action='store_true', help="Focus on viruses.")
    parser.add_argument("--bacteria", action='store_true', help="Focus on bacteria.")

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    for forward in glob.glob(os.path.join(args.input_dir, "*_R1.fastq*")):
        base_name = os.path.basename(forward).replace("_R1.fastq.gz", "").replace("_R1.fastq", "")
        reverse = os.path.join(args.input_dir, f"{base_name}_R2.fastq.gz") if forward.endswith(".gz") else os.path.join(args.input_dir, f"{base_name}_R2.fastq")
        if not os.path.isfile(reverse):
            reverse = None

        process_sample(forward, reverse, base_name, args.bowtie2_index, args.kraken_db, args.output_dir, args.threads)

    merged_tsv_path = aggregate_kraken_results(args.output_dir, args.metadata_file, args.read_count, args.top_N, args.virus, args.bacteria)
    generate_abundance_plots(merged_tsv_path, args.virus, args.bacteria, args.top_N)

if __name__ == "__main__":
    main()
