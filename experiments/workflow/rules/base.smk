checkpoint get_cluster_sizes:
    input:
        nodeDescription=INPUT_FOLDER
        / "{SAMPLE}"
        / "{SAMPLE}_samples_nodeDescription.tsv",
    output:
        INPUT_FOLDER / "{SAMPLE}" / "{SAMPLE}_clusterSizes.csv",
    log:
        PROJECT_DIR / "logs" / "get_cluster_sizes.{SAMPLE}.log",
    conda:
        Path(workflow.basedir) / "envs" / "R.yml"
    shell:
        """
        RScript {SCRIPT_DIR}/createInputSummary.R -i {input.nodeDescription} -n {wildcards.SAMPLE}
        """


rule simulate_clusters:
    input:
        nodeDescription=INPUT_FOLDER
        / "{SAMPLE}"
        / "{SAMPLE}_samples_nodeDescription.tsv",
    output:
        SIMULATION_FOLDER
        / "{SAMPLE}_{SIMSIZE}"
        / "{SAMPLE}_{SIMSIZE}_samples_nodeDescription.tsv",
    log:
        PROJECT_DIR / "logs" / "simulate_clusters.{SAMPLE}_{SIMSIZE}.log",
    conda:
        Path(workflow.basedir) / "envs" / "R.yml"
    shell:
        """
        RScript {SCRIPT_DIR}/simulateClusters.R -i {input.nodeDescription} -n {wildcards.SAMPLE} -s {SIMSIZE}
        """


####TODO: replace the parameters for -p
###TODO: replace the parameter for -o
###### UNFINISHED rule run_CTC_SCITE:
#    input:
#        variants=INPUT_FOLDER / "{SAMPLE}" / "{SAMPLE}.txt",
#        nodeDescription=INPUT_FOLDER
#        / "{SAMPLE}"
#        / "{SAMPLE}_samples_nodeDescription.tsv",
#    output:
#        INPUT_FOLDER / "{SAMPLE}" / "{SAMPLE}_1M_{iterator}_seed{seed}_postSampling.tsv",
#    params:
#        executable= PROJECT_DIR / ".." / "CTC_SCITE",
#        chainlength = config.get(chainlength, 1000000),
#    log:
#        PROJECT_DIR / "logs" / "run_CTC_SCITE.{SAMPLE}.log",
#    shell:
#        """
#        {params.executable}/./CTC_SCITE -i {input.variants} - r 1 -l {params.chainlength} -g 1 -seed {wildcards.seed} -description {input.nodeDescription} -e 0.2 -p {params.chainlength}/5000 0.2* {params.chainlength} -samples {input.nodeDescription} -o {output}
#        """


###TODO: define the input function and the shell directive
rule postprocess_tree_sampling:
    input:
        get_CTC_SCITE_runs,
    output:
        INPUT_FOLDER / "{SAMPLE}" / "{SAMPLE}_postSampling.tsv",


rule tree_sampling:
    input:
        variants=INPUT_FOLDER / "{SAMPLE}" / "{SAMPLE}.txt",
        trees=INPUT_FOLDER / "{SAMPLE}" / "{SAMPLE}_postSampling.tsv",
        nodeDescription=INPUT_FOLDER
        / "{SAMPLE}"
        / "{SAMPLE}_samples_nodeDescription.tsv",
    output:
        INPUT_FOLDER / "{SAMPLE}" / "{SAMPLE}_postSampling_splits.txt",
    resources:
        mem_mb=4096,
        runtime=2880,
    params:
        code_dir=SCRIPT_DIR / "CTC_parseSampling",
    log:
        PROJECT_DIR / "logs" / "treeSampling.{SAMPLE}.log",
    shell:
        """
        {params.code_dir}/./CTC_parseSampling -i {input.variants} -samples {input.trees} -description {input.nodeDescription}
        """


rule get_coverage_scores:
    input:
        readCounts=INPUT_FOLDER / "{SAMPLE}" / "{SAMPLE}.txt",
    output:
        INPUT_FOLDER / "{SAMPLE}" / "{SAMPLE}_covScore.txt",
    resources:
        mem_mb=2048,
        runtime=60,
    params:
        code_dir=SCRIPT_DIR,
    log:
        PROJECT_DIR / "logs" / "get_coverage_scores.{SAMPLE}.log",
    shell:
        """
        python {params.code_dir}/getPositionCoverageScore.py {input.readCounts}
        """


# rule prepare_markdown_file:
#     input:
#         PROJECT_DIR / "workflow" / "resources" / "template.Rmd",
#     output:
#         PROJECT_DIR / "data" / "markdowns" / "{SAMPLE}_treeSampling.Rmd",
#     resources:
#         mem_mb=1024,
#         runtime=10,
#     params:
#         sampling_depth=1000,
#         script_dir={markdown_helper_functions},
#         author=config["author"],
#         input_dir=input_folder,
#         simulation_dir=simulations_folder,
#     log:
#         PROJECT_DIR / "logs" / "prepare_markdown_file.{SAMPLE}.log",
#     shell:
#         """
#         ( sed -e "s/__tree__/{wildcards.SAMPLE}/g" \
#         -e "s/__nSamplingEvents__/{params.sampling_depth}/g" \
#         -e "s/__date__/$(date +'%Y-%m-%d')/g" \
#         -e "s/__functionsScript__/{params.script_dir}/g" \
#         -e "s/__author__/{params.author}/g" \
#         -e "s/__input_dir__/{params.input_dir}/g" \
#         -e "s/__sim_dir__/{params.simulation_dir}/g" \
#         {input} > {output} ) &> {log}
#         """


rule render_markdown_file:
    input:
        PROJECT_DIR / "workflow" / "resources" / "template.Rmd",
        simulations=get_simulation_files,
    output:
        PROJECT_DIR / "data" / "htmls" / "{SAMPLE}_treeSampling.html",
    params:
        sampling_depth=1000,
        script_dir=RESOURCES_DIR / "functions.R",
        author=config["author"],
        input_dir=INPUT_FOLDER,
        simulation_dir=SIMULATION_FOLDER,
    resources:
        mem_mb_per_cpu=4096,
        runtime=2880,
    threads: 16
    conda:
        Path(workflow.basedir) / "envs" / "R.yml"
    log:
        PROJECT_DIR / "logs" / "render_markdown_file.{SAMPLE}.log",
    shell:
        """
        ( Rscript -e "rmarkdown::render('{input}', output_file = '{output}', params = list(inputFolder = '{params.input_dir}', nSamplingEvents={params.sampling_depth}, simulationInputFolder = '{params.simulation_dir}', treeName = '{wildcards.SAMPLE}', functionsScript = '{params.script_dir}'))" ) &> {log}
        """


rule render_topSeparators_markdowns:
    input:
        markdown=PROJECT_DIR / "data" / "markdowns" / "{SAMPLE}_topSeparators.Rmd",
        coverageScore=INPUT_FOLDER / "{SAMPLE}" / "{SAMPLE}_covScore.txt",
        clustersplits=INPUT_FOLDER / "{SAMPLE}" / "{SAMPLE}_postSampling_splits.txt",
    output:
        PROJECT_DIR / "data" / "htmls" / "{SAMPLE}_topSeparators.html",
    resources:
        mem_mb=6144,
        runtime=2880,
    params:
        markdowndir=MARKDOWNS,
    threads: 4
    conda:
        Path(workflow.basedir) / "envs" / "R.yml"
    log:
        PROJECT_DIR / "logs" / "render_topSeparators_markdowns.{SAMPLE}.log",
    shell:
        """
        cd {params.markdowndir}
        ( Rscript -e "rmarkdown::render('{input.markdown}', output_file = '{output}')" ) &> {log}
        """
