#### The script takes the merged splitting summary table and maps the number annotations back to the original names. Furthermore, it filters out the
#### clusters that are in fact no CTC clutsers put primary tumor samples 

import pandas as pd
from pathlib import Path






def load_node_description(sample_name, color_annotation):
    file_path = Path('/home/jovyan/work/ctc-data/WES_experiment/splitting_summaries/') / Path(sample_name) / f"{sample_name}_samples_nodeDescription.tsv"
    df = pd.read_csv(file_path, sep='\t', header = None)
    return df

def search_rows_with_color_annotation(df, color_annotation):
    color_annotation = f"={color_annotation},"
    return ",".join(df[df.iloc[:, -1].str.contains(color_annotation, na=False)].iloc[:, 0].values.tolist())

def remove_primary_tumor_samples(line, master_table):
    samples = line.split(',')
    sample_constituents = []
    for s in samples:
        if not 'CTC' in master_table.loc[master_table['sample_id'] == s, 'sample_type'].to_string():
            continue
        else:
            sample_constituents.append(s)
    remaining_samples = ','.join(sample_constituents)
    return remaining_samples


def filter_out_non_CTC_clusters(splitting_summary, master_table):
    splitting_summary_filtered = splitting_summary.copy(deep=True) 
    splitting_summary_filtered['sample_names'] = splitting_summary_filtered['sample_names'].map(lambda row: remove_primary_tumor_samples(row,master_table))
    rows_to_keep = [idx for idx, row in splitting_summary_filtered.iterrows() if len(row['sample_names'])>0]
    splitting_summary_filtered = splitting_summary_filtered.loc[rows_to_keep,]

    return splitting_summary_filtered

def compute_number_of_cells(sample_names, master_table, WBC = False):
    all_sample_names = sample_names.split(',')

    if WBC:
        cell_count_column = 'n_wbc_attached'
    else:
        cell_count_column = 'n_cells'
    n_cells = 0
    for sample in all_sample_names:
        n_cells += int(master_table.loc[master_table['sample_id'] == sample, cell_count_column].iloc[0])
    return n_cells


def main():
    splitting_summary = pd.read_csv('/home/jovyan/work/ctc-data/WES_experiment/splitting_summaries/splittingSummary_full.tsv', sep = '\t', header = 0)
    splitting_summary['sample_names'] = splitting_summary.apply(
            lambda row: search_rows_with_color_annotation(load_node_description(row['Sample Name'], row['Color']), row['Color']),
            axis=1
        )

    master_table = pd.read_csv('/home/jovyan/work/CTC-SCITE/experiments/preprocessing/gDNA_sample_annotation_final.tsv', sep = '\t')

    splitting_summary['n_cells'] = splitting_summary.apply(
        lambda row: compute_number_of_cells(row['sample_names'],master_table),
        axis = 1
    )

    splitting_summary['n_wbcs'] = splitting_summary.apply(
        lambda row: compute_number_of_cells(row['sample_names'],master_table, WBC = True),
        axis = 1
    )


    splitting_summary = filter_out_non_CTC_clusters(splitting_summary,master_table)

    splitting_summary = splitting_summary.loc[splitting_summary['Sample Name'] != 'Br26',:]
    splitting_summary.to_csv('splittingSummary_full_with_sample_names.tsv', sep = '\t', index = False)

if __name__ == '__main__':
    main()
