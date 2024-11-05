import pandas as pd
import random
from collections import defaultdict
import plotly.express as px
import plotly.io as pio
import os
from .trimmomatic import run_trimmomatic
from .bowtie2 import run_bowtie2
from .kraken2 import run_kraken2
import distinctipy
import numpy as np
import matplotlib.pyplot as plt
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



def generate_sample_ids_csv(kraken_dir):
    """
    Generates a CSV file containing sample IDs extracted from Kraken report filenames.

    Parameters:
    - kraken_dir (str): Path to the directory containing Kraken report files.

    Returns:
    - str: Path to the generated sample_ids.csv file.
    """
    try:
        # Extract sample IDs from Kraken report filenames
        sample_ids = []
        for file_name in os.listdir(kraken_dir):
            if file_name.endswith("_report.txt"):
                sample_id = '_'.join(file_name.split('_')[:-1])
                sample_ids.append(sample_id)

        # Save sample IDs to CSV
        sampleid_df = pd.DataFrame(sample_ids, columns=['Sample_IDs'])
        sampleid_csv_path = os.path.join(kraken_dir, "sample_ids.csv")
        sampleid_df.to_csv(sampleid_csv_path, index=False)
        
        print(f"Sample IDs saved to {sampleid_csv_path}")
        return sampleid_csv_path

    except Exception as e:
        print(f"Error generating sample_ids.csv: {e}")
        return None
        
def aggregate_kraken_results(kraken_dir, metadata_file=None, sample_id_df=None, read_count=10):
    """
    Aggregates Kraken results, merging metadata or using sample IDs if metadata is not provided.

    Parameters:
    - kraken_dir (str): Path to the directory containing Kraken report files.
    - metadata_file (str, optional): Path to the metadata CSV file. Defaults to None.
    - sample_id_df (DataFrame, optional): DataFrame of sample IDs. Used if metadata_file is not provided.
    - read_count (int): Minimum read count threshold for filtering results.
    
    Returns:
    - str: Path to the generated merged TSV file.
    """
    try:
        # Load metadata from file first, fall back to sample_id_df if not provided
        if metadata_file:
            metadata = pd.read_csv(metadata_file, sep=",")
            print("Using metadata from the provided metadata file.")
        elif sample_id_df is not None:
            metadata = sample_id_df
            print("Using sample IDs as metadata.")
        else:
            raise ValueError("Either metadata_file or sample_id_df must be provided.")

        sample_id_col = metadata.columns[0]  # Assume the first column is the sample ID

        # Dictionary to store aggregated results
        aggregated_results = {}

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

                        # Check if rank code is species-level and meets the read count threshold
                        if (rank_code == 'S' or rank_code == 'S1' or rank_code == 'S2' or rank_code == 'S3') and nr_frag_direct_at_taxon >= read_count:
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
        merged_tsv_path = os.path.join(kraken_dir, "merged_kraken.tsv")
        with open(merged_tsv_path, 'w') as f:
            # Write headers dynamically
            headers = ['Perc_frag_cover', 'Nr_frag_cover', 'Nr_frag_direct_at_taxon', 'Rank_code', 'NCBI_ID', 'Scientific_name', 'SampleID'] + metadata.columns[1:].tolist()
            f.write("\t".join(headers) + "\n")
            for sampleandtaxonid, data in aggregated_results.items():
                f.write("\t".join(str(data[col]) for col in headers) + "\n")

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
                # Create a color mapping based on unique values in the 'focus' column
                #colordict = dict(zip(grouped_sum[focus].unique(), distinctipy.get_colors(len(grouped_sum[focus].unique()))))
                colordict = defaultdict(int)
                #random_colors = ["#{:06X}".format(random.randint(0, 0xFFFFFF)) for _ in range(len(grouped_sum[focus].unique()))]
                random_colors = ['#{:06x}'.format(random.randint(0xff0000, 0xffff00)) for _ in range(len(grouped_sum[focus].unique()))]
                
                #for target, color in zip(grouped_sum[focus].unique(), random_colors):
                for target, color in zip(grouped_sum[focus].unique(), random_colors):
                    colordict[target] = color
               
                #colordict=distinctipy.get_colors(len(grouped_sum[col].unique()))
                #colordict = dict(zip(grouped_sum[focus].unique(), distinctipy.get_colors(len(grouped_sum[focus].unique()))))
                # Generate a unique color for each unique item in the 'focus' column
                #random_colors = distinctipy.get_colors(len(grouped_sum[col].unique()))
                #colordict = {focus: color for category, color in zip(grouped_sum[col].unique(), random_colors)}
                def color_distance(c1, c2):
                  return np.sqrt(sum((a - b) ** 2 for a, b in zip(c1, c2)))

                # Function to generate colors with a minimum distance
               # Function to generate unique colors with a minimum distance
                def generate_distant_colors(num_colors, min_distance=30):
                   colors = []
                   seen_colors = set()  # To track already generated colors in tuple format
                   while len(colors) < num_colors:
                      new_color = (random.randint(0, 255), random.randint(0, 255), random.randint(0, 255))
        
                       # Check if the new color is sufficiently distant and not already seen
                      if all(color_distance(new_color, existing) >= min_distance for existing in colors) and new_color not in seen_colors:
                         colors.append(new_color)
                         seen_colors.add(new_color)  # Add the new color to the set
                   # Convert RGB colors to hex format
                   hex_colors = [f"#{r:02x}{g:02x}{b:02x}" for r, g, b in colors]
                   return hex_colors
                #num_categories = len(grouped_sum[focus].unique())
                # Generate distinct colors for each category
               # colors = generate_distant_colors(num_categories, min_distance=90)
                #colors = base_colors[:len(unique_targets)]
                #colors = plt.get_cmap('tab20').colors 
                #colordict = dict(zip(grouped_sum[focus].unique(), colors))
                # Generate colors using 'tab20' colormap
                #colors = plt.get_cmap('tab20').colors

                  # Ensure there are enough colors for unique items in grouped_sum[focus]
                #unique_targets = grouped_sum[focus].unique()
                #if len(unique_targets) > len(base_colors):
                  #raise ValueError("Not enough unique colors in 'tab20' colormap for the number of unique targets.")
                  #colors = plt.get_cmap('tab20')(np.linspace(0, 1, len(unique_targets)))
                #else:
                # Use base 'tab20' colors if they are sufficient
                 # colors = base_colors[:len(unique_targets)]
                 # Map each unique target to a color from the colormap
                #olordict = dict(zip(grouped_sum[focus].unique(), colors[:len(unique_targets)]))
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
