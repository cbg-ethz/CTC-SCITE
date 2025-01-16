from pathlib import Path
import logging

import numpy as np
import pandas as pd
import math
from matplotlib import pyplot as plt
from scipy.stats import combine_pvalues

# Generate a list of barcode counts files

logging.basicConfig(level=logging.INFO)


def load_cluster_data(path, pattern, cluster2tumor):
    #files = [str(file) for file in Path(path).glob('*.csv') if pattern in str(file)]
    files = cluster2tumor.loc[cluster2tumor['tumor_sample'] == pattern, 'basename']
    # Create empty dictionary to store data frames
    df_dict = {}
    cluster_sizes = []
    # Loop through all sample barcode counts files
    for idx,file in enumerate(files):
        basename = Path(f'{file}.csv').stem
        #pattern = next((p for p in ["910", "905", "904", "903", "902", "141", "140"] if p in basename), None)
        
        cluster_size = cluster2tumor.loc[cluster2tumor['basename'] == file, 'cluster_size'].iloc[0]
    
        
        try:
            # Load counts into data frame
            # Use pd.read_csv with chunksize for faster processing
            df = pd.read_csv(Path(path) / f'{file}.csv', sep='\t', header=None, dtype={0: 'str'})
            logging.info(f'File {file} loaded')
        except Exception as e:
            logging.error(f"Error reading {file} - skipping\n")
            logging.error(e)
            continue
        df.columns = ['Barcode ID', '1', 'Barcode Count', '3']
             


        # Sort by decreasing barcode counts
        df_sorted = df.sort_values(by=df.columns[2], ascending=False)
        df_sorted['Mouse ID'] = pattern
        df['Cluster Size'] = cluster_size
        df['CTC Cluster ID'] = id ### The choice of cluster ID is completely arbitrary and only helps to distinguish between different clusters
        # Calculate fraction of total counts for each Barcode
        df_sorted['prop_col'] = df_sorted.iloc[:, 2] / df_sorted.iloc[:, 2].sum()
        
        # Calculate cumulative barcode proportions
        df_sorted['cumprop_col'] = df_sorted['prop_col'].cumsum()
        
        # Find the largest row number for which cumprop_col < 0.9
        
        n_clones = len(df_sorted[df_sorted['cumprop_col'] < 0.9].index) + 1
        df_sorted['Clone present in cluster'] = False
        df_sorted.loc[df_sorted.index[:n_clones], 'Clone present in cluster'] = True

        # do a quality check to exclude all clusters for which the number of clones exceeds the number of cells in the cluster


        if not df_sorted['Clone present in cluster'].sum() > cluster_size:
            df_dict[basename] = df_sorted
            cluster_sizes.append(cluster_size)
            logging.info(f'File {file} added to dictionary')

    return df_dict, list(set(cluster_sizes))


def load_primary_data(path):
    file_list = [filename for filename in list(path.glob("combined_*_filter_merge.rds.csv"))]
    #file_list = ["combined_140_order_filter_merge.rds.csv", "combined_141_order_filter_merge.rds.csv", "combined_902_order_filter_merge.rds.csv", "combined_903_order_filter_merge.rds.csv", "combined_904_order_filter_merge.rds.csv", "combined_905_order_filter_merge.rds.csv", "combined_910_order_filter_merge.rds.csv"]
    primary_dict = {}
    for file in file_list:
        mouse_ID = file.stem.split('_')[1]
        primary_dict[mouse_ID] = pd.read_csv(file, sep=',', header=0)

        primary_dict[mouse_ID]['observations'] = primary_dict[mouse_ID]['observations'].astype(str)
        primary_dict[mouse_ID]['prop_av'] = primary_dict[mouse_ID]['prop_av'].astype(float)
    return primary_dict


def preprocess_primary_data(primary_data):
    # Determine the cutoff for each mouse_ID

    for mouse_ID, data in primary_data.items():
        cutoff = data['prop_av'].quantile(0.01)
        data_sorted = data.sort_values(by="prop_av", ascending=False)
        
        # Calculate cumulative barcode proportions
        data_sorted['cumulative_proportions'] = data_sorted['prop_av'].cumsum()
        
        # Find the largest row number for which cumprop_col < 0.9
        
        n_clones = len(data_sorted[data_sorted['prop_av'] > cutoff].index)+1
        data_sorted = data_sorted.iloc[:n_clones,:]
        primary_data.update({mouse_ID: data_sorted})

    return primary_data

    # Calculate fraction of total counts for each Barcode


def merge_data(files, primary_data):
    for mouse_ID, data in primary_data.items():
        for cluster_ID, cluster_data in files.items():
            if str(cluster_data['Mouse ID'].iloc[0]) == mouse_ID:
                # Filter out barcode IDs in cluster_data that are not in primary_data
                
                cluster_data = cluster_data[cluster_data['Barcode ID'].isin(data['observations'])]

                # Check for barcode IDs with "Clone present in cluster" True but not in observations
                missing_barcodes = cluster_data[(cluster_data['Clone present in cluster'] == True) & (~cluster_data['Barcode ID'].isin(data['observations']))]
                if not missing_barcodes.empty:
                    logging.warning(f"Mouse ID {mouse_ID}, Cluster ID {cluster_ID} has barcodes with 'Clone present in cluster' True but not in observations: {missing_barcodes['Barcode ID'].tolist()}")

                # Merge the prop_av values from primary_data into cluster_data
                cluster_data = cluster_data.merge(data[['observations', 'prop_av']], left_on='Barcode ID', right_on='observations', how='left')

                cluster_data.drop(columns=['observations'], inplace=True)
 
                files.update({cluster_ID: cluster_data})
    return files


def compute_G_score(clones_in_cluster, n_cells):
    return 2*clones_in_cluster["prop_av"].map(lambda x: np.log(1/(1-(1-x) ** n_cells))).sum()



def simulate_G_scores(proportions, n_cells, n_simulations):
            G_scores = []
            probabilities = proportions/np.sum(proportions)

            for _ in range(n_simulations):
                simulated_counts = np.random.multinomial(n_cells, probabilities)
                simulated_data = pd.DataFrame({
                    'simulated_counts': simulated_counts,
                    'prop_av': proportions
                })
                simulated_data['Clone present in cluster'] = simulated_data['simulated_counts'].apply(lambda x: 1 if x > 0 else 0)
                simulated_data = simulated_data[simulated_data['simulated_counts'] > 0]
                G_score = compute_G_score(simulated_data, n_cells)
                G_scores.append(G_score)
            
            plt.hist(G_scores, bins=60, edgecolor='k', alpha=0.7)
            plt.xlabel('G Score')
            plt.ylabel('Frequency')
            plt.title('Histogram of Simulated G Scores')
            plt.show()
            
            return G_scores


def compute_test(all_clones, n_cells_in_cluster, simulations = None):

    resolution_of_simulation = 10000
    if simulations is None:
        G_scores = simulate_G_scores(all_clones['prop_av'], n_cells_in_cluster, resolution_of_simulation)
    else:
        G_scores = simulations
    clones_in_cluster = all_clones[all_clones['Clone present in cluster'] == True]
    G = compute_G_score(clones_in_cluster, n_cells_in_cluster)
    p_value = sum(G_scores >= G)/resolution_of_simulation
    return(p_value)
    #p_value = clones_in_cluster["prop_av"].prod()* math.factorial(n_cells_in_cluster)/math.factorial(clones_in_cluster.shape[0])




if __name__ == '__main__':
    path = '/home/jovyan/work/ctc-data/barcoding_experiment/Cluster data'
    primary_data_path = '/home/jovyan/work/ctc-data/barcoding_experiment/combined_primary_cluster'
    primary_data = load_primary_data(Path(primary_data_path))
    primary_data = preprocess_primary_data(primary_data)
    cluster2tumor = pd.read_csv('/home/jovyan/work/ctc-data/barcoding_experiment/summary_df_filter_final_all.csv')
    tumor_samples = list(set(cluster2tumor['tumor_sample']))
    # Create empty dictionary to store data frames
    p_values = []
    mouse_models = []
    number_of_simulations = 10000
    for pattern in tumor_samples:
        files, cluster_sizes = load_cluster_data(path, pattern, cluster2tumor)

        merged_data = merge_data(files, primary_data)

        first_dataset = next(iter(merged_data.values()))
        
        simulated_G_scores = pd.DataFrame(np.zeros((number_of_simulations, len(cluster_sizes))), columns=[str(cluster_size) for cluster_size in cluster_sizes])
        
        for cell_number in cluster_sizes:
            if simulated_G_scores[str(cell_number)].sum() == 0:
                logging.info(f"Simulating null distribution for cluster size:  {cell_number}")
                simulated_G_scores[str(cell_number)] = simulate_G_scores(first_dataset['prop_av'], cell_number, number_of_simulations)

        for cluster_id, cluster_data in merged_data.items():
            n_cells = cluster_id.split('_')[1]
            p_values.append(compute_test(cluster_data, int(n_cells), simulated_G_scores[n_cells])+1e-16)
            mouse_models.append(pattern)
            logging.info(f"P value for cluster {cluster_id}: {p_values[-1]}")

    p_value_summary = pd.DataFrame({'P value': p_values, 'Mouse Model': mouse_models})
    logging.info(p_value_summary)
    p_value_summary.to_csv('/home/jovyan/work/ctc-data/barcoding_experiment/')

    # Combine p-values using Fisher's method
    combined_p_value = combine_pvalues(p_values, method='fisher')[1]
    print(f'Combined p-value: {combined_p_value}')
    plt.hist(p_values, bins=30, edgecolor='k', alpha=0.7)
    plt.xlabel('P values')
    plt.ylabel('Frequency')
    plt.title('Histogram of p_values')
    plt.show()
    # Plot the p_values stratified by Mouse model
    plt.figure(figsize=(10, 6))
    for mouse_model in p_value_summary['Mouse Model'].unique():
        subset = p_value_summary[p_value_summary['Mouse Model'] == mouse_model]
        plt.hist(subset['P value'], bins=30, alpha=0.5, label=f'Mouse Model {mouse_model}')

    plt.xlabel('P values')
    plt.ylabel('Frequency')
    plt.title('Histogram of P values stratified by Mouse Model')
    plt.legend()
    plt.show()

    """
    for cluster_id, cluster_data in merged_data.items():
        for idx in range(15):
            simulate_G_scores(cluster_data['prop_av'], idx, 1000)
        n_cells = int(cluster_id.split('_')[1])
        clones_in_cluster = cluster_data[cluster_data['Clone present in cluster'] == True]
        G_score = compute_G_score(clones_in_cluster, n_cells)
        print(f'Cluster ID: {cluster_id}, G_score: {G_score}')
        simulate_G_scores(cluster_data['prop_av'], n_cells, 100000)
    """