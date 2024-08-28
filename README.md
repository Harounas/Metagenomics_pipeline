# kraken_and_abundanceplots.py is used to perform taxonomic classification using Kraken and visualize viral and bacterial abundances according to information available in the provided metadata. This pipeline assumes the first column in the metadata corresponds to the sample IDs

run_kr_abundance --input_dir /home/harouna/ARSNACAdata/bamfiles/mypipeline --output_dir kraken_summary_files --kraken_db /home/harouna/Mydatabses/Standard --bowtie2_index /home/harouna/240405_VH01476_5_AAC3HYYHV_batch2/GRCh38_noalt_as --metadata_file /home/harouna/ARSNACAdata/bamfiles/METADATAARSNACAB.09.07.2024.csv --threads 8 --virus --read_count 10 --top_N 100

* --kraken_db: /path/to/kraken_db
* --bowtie2_index: /path/to/bowtie2_index
* --output_dir /path/to/output
* --input_dir /path/to/input_fastq_files
* --metadata_file: /path/to/metadata.csv
* --read_count: minimum read count
* --top_N N: select top N viral or bacterial species
* --virus: we are considering only virus. Use --bacteria if you need only bacteria


  #Installation
  git clone https://github.com/Harounas/Metagenomics_pipeline.git
  cd Metagenomics_pipeline
   sudo pip install -e .







