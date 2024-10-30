def aggregate_kraken_results(kraken_dir, metadata_file=None, sample_id_df=None, read_count=0):
    try:
        # Check if metadata_file is provided and exists
        if metadata_file and os.path.exists(metadata_file):
            metadata = pd.read_csv(metadata_file, sep=",")
            sample_id_col = metadata.columns[0]  # Assume the first column is the sample ID
        elif sample_id_df is not None:
            metadata = sample_id_df
            sample_id_col = metadata.columns[0]  # Assume the first column is the sample ID
        else:
            raise ValueError("No metadata provided, and sample_id_df is None.")

        # Dictionary to store aggregated results
        aggregated_results = {}
        sampleid = []
        
        # Iterate over each Kraken report file
        for file_name in os.listdir(kraken_dir):
            if file_name.endswith("_report.txt"):
                with open(os.path.join(kraken_dir, file_name), 'r') as f:
                    for line in f:
                        fields = line.strip().split('\t')
                        perc_frag_cover = fields[0]
                        nr_frag_cover = fields[1]
                        nr_frag_direct_at_taxon = int(fields[2])
                        rank_code = fields[3]
                        ncbi_ID = fields[4]
                        scientific_name = fields[5]
                        parts = file_name.split('_')
                        extracted_part = '_'.join(parts[:-1])
                        sampleandtaxonid = extracted_part + str(ncbi_ID)
                        sampleid.append(extracted_part)
                        if rank_code == 'S' and nr_frag_direct_at_taxon >= read_count:
                            if extracted_part in metadata[sample_id_col].unique():
                                sample_metadata = metadata.loc[metadata[sample_id_col] == extracted_part].iloc[0].to_dict()
                                aggregated_results[sampleandtaxonid] = {
                                    'Perc_frag_cover': perc_frag_cover,
                                    'Nr_frag_cover': nr_frag_cover,
                                    'Nr_frag_direct_at_taxon': nr_frag_direct_at_taxon,
                                    'Rank_code': rank_code,
                                    'NCBI_ID': ncbi_ID,
                                    'Scientific_name': scientific_name,
                                    'SampleID': extracted_part,
                                    **sample_metadata
                                }

        # Output aggregated results to a TSV file
        merged_tsv_path = os.path.join(kraken_dir, "merged_kraken1.tsv")
        with open(merged_tsv_path, 'w') as f:
            # Write headers dynamically
            headers = ['Perc_frag_cover', 'Nr_frag_cover', 'Nr_frag_direct_at_taxon', 'Rank_code', 'NCBI_ID', 'Scientific_name', 'SampleID'] + metadata.columns[1:].tolist()
            f.write("\t".join(headers) + "\n")
            for sampleandtaxonid, data in aggregated_results.items():
                f.write("\t".join(str(data[col]) for col in headers) + "\n")

        # Save the sampleid list as a CSV file
        sampleid_df = pd.DataFrame(sampleid, columns=['Sample_IDs'])
        sampleid_csv_path = os.path.join(kraken_dir, "sample_ids.csv")
        sampleid_df.to_csv(sampleid_csv_path, index=False)

        return merged_tsv_path

    except Exception as e:
        print(f"Error aggregating Kraken results: {e}")
        return None
