import pandas as pd
import random


def get_simulation_files(wildcards):
    clusterSizes = pd.load_csv(
        checkpoints.get_cluster_sizes.get(SAMPLE=wildcards.SAMPLE).output, header=0
    )
    simulations = []
    for _, row in clusterSizes.iterrows():
        simulations.append(
            SIMULATION_FOLDER
            / f"{wildcards.SAMPLE}_{row[0]}"
            / f"{wildcards.SAMPLE}_{row[0]}_postSampling.tsv"
        )
    return simulations
