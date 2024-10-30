import pandas as pd
import random
from collections import defaultdict
import plotly.express as px
import plotly.io as pio
import os
from .trimmomatic import run_trimmomatic
from .bowtie2 import run_bowtie2
from .kraken2 import run_kraken2

def process_sample(forward, reverse, base_name, bowtie2_index, kraken_db, output_dir, threads, run_bowtie, use_precomputed_reports):
    try:
        if not use_precomputed_reports:
            # Step 1: Run Trimmomatic
            trimmed_forward, trimmed_reverse = run_trimmomatic(forward, reverse, base_name, output_dir, threads)

            # Step 2: Optionally run Bowtie2 to deplete host genome reads
            if run_bowtie:
                unmapped_r1, unmapped_r2 = run_bowtie2(trimmed_forward, trimmed_reverse, base_name, bowtie2_index, output_dir, threads)
            else:
                unmapped_r1, unmapped_r2 = trimmed_forward, trimmed_reverse

            # Step 3: Run Kraken2 with the reads
            kraken_report = run_kraken2(unmapped_r1, unmapped_r2, base_name, kraken_db, output_dir, threads)
        else:
            # Use the precomputed Kraken2 report
            kraken_report = os.path.join(output_dir, f"{base_name}_report.txt")
            if not os.path.exists(kraken_report):
                raise FileNotFoundError(f"Precomputed Kraken2 report not found: {kraken_report}")

        return kraken_report

    except Exception as e:
        print(f"Error processing sample {base_name}: {e}")
        return None



def aggregate_kraken_results(kraken_dir, metadata_file, read_count):
    try:
        metadata = pd.read_csv(metadata_file, sep=",")
        sample_id_col = metadata.columns[0]  # Assume the first column is the sample ID

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


def generate_abundance_plots(merged_tsv_path, top_N):
    try:
        df = pd.read_csv(merged_tsv_path, sep="\t")
        df.columns = df.columns.str.replace('/', '_').str.replace(' ', '_')
        df = df.apply(lambda col: col.map(lambda x: x.strip() if isinstance(x, str) else x))
        df = df[df['Scientific_name'] != 'Homo sapiens']  # Remove human reads

        # Generate both viral and bacterial abundance plots
        for focus, filter_str, plot_title in [
            ('Virus_Type', 'Virus', 'Viral'),
            ('Bacteria_Type', 'Virus', 'Bacterial')
        ]:
            if focus == 'Bacteria_Type':
                df_focus = df[~df['Scientific_name'].str.contains(filter_str, case=False, na=False)]
            else:
                df_focus = df[df['Scientific_name'].str.contains(filter_str, case=False, na=False)]
            df_focus = df_focus.rename(columns={'Scientific_name': focus})

            if top_N:
                top_N_categories = df_focus[focus].value_counts().head(top_N).index
                df_focus = df_focus[df_focus[focus].isin(top_N_categories)]

            categorical_cols = df_focus.select_dtypes(include=['object']).columns.tolist()
            categorical_cols.remove(focus)

            for col in categorical_cols:
                grouped_sum = df_focus.groupby([focus, col])['Nr_frag_direct_at_taxon'].mean().reset_index()

                colordict = defaultdict(int)
                random_colors = ["#{:06x}".format(random.randint(0, 0xFFFFFF)) for _ in range(len(grouped_sum[col].unique()))]
                for target, color in zip(grouped_sum[focus].unique(), random_colors):
                    colordict[target] = color

                plot_width = 1100 + 5 * len(grouped_sum[col].unique())
                plot_height = 800 + 5 * len(grouped_sum[col].unique())
                font_size = max(10, 14 - len(grouped_sum[col].unique()) // 10)

                fig = px.bar(
                    grouped_sum,
                    x=col,
                    y='Nr_frag_direct_at_taxon',
                    color=focus,
                    color_discrete_map=colordict,
                    title=f"{plot_title} Abundance by {col}"
                )

                fig.update_layout(
                    xaxis=dict(tickfont=dict(size=font_size), tickangle=45),
                    yaxis=dict(tickfont=dict(size=font_size)),
                    title=dict(text=f'Average {plot_title} Abundance by {col}', x=0.5, font=dict(size=16)),
                    bargap=0.5,
                    legend=dict(
                        font=dict(size=font_size),
                        x=1,
                        y=1,
                        traceorder='normal',
                        orientation='v',
                        itemwidth=30,
                        itemsizing='constant',
                        itemclick='toggleothers',
                        itemdoubleclick='toggle'
                    ),
                    width=plot_width,
                    height=plot_height
                )

                fig.write_image(f"{plot_title}_Abundance_by_{col}.png", format='png', scale=3)

    except Exception as e:
        print(f"Error generating abundance plots: {e}")
