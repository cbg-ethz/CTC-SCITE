rule tree_sampling:
    input:
        variants= INPUT_FOLDER / "{SAMPLE}" / "{SAMPLE}.txt",
        trees= INPUT_FOLDER / "{SAMPLE}" / "{SAMPLE}_postSampling.tsv",
        nodeDescription= INPUT_FOLDER / "{SAMPLE}" / "{SAMPLE}_samples_nodeDescription.tsv",
    output:
         INPUT_FOLDER / "{SAMPLE}" / "{SAMPLE}_postSampling_splits.txt" 
    resources:
        mem_mb=4096,
        runtime=2880,
    params:
        code_dir=  SCRIPT_DIR / "CTC_parseSampling",
    log:
        PROJECT_DIR / "logs" / "treeSampling.{SAMPLE}.log",
    shell:
        """
        {params.code_dir}/./CTC_parseSampling -i {input.variants} -samples {input.trees} -description {input.nodeDescription}
        """

rule get_coverage_scores:
    input:
        readCounts=INPUT_FOLDER / "{SAMPLE}" / "{SAMPLE}.txt"
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

rule prepare_markdown_file:
    input:
        PROJECT_DIR / "workflow" / "resources" / "template.Rmd",
    output:
        PROJECT_DIR / "data" / "markdowns" / "{SAMPLE}_treeSampling.Rmd",
    resources:
        mem_mb=1024,
        runtime=10,
    params:
        sampling_depth=1000,
        script_dir={markdown_helper_functions},
        author=config["author"],
        input_dir=input_folder,
        simulation_dir=simulations_folder,
    log:
        PROJECT_DIR / "logs" / "prepare_markdown_file.{SAMPLE}.log",
    shell:
        """
        ( sed -e "s/__tree__/{wildcards.SAMPLE}/g" \
        -e "s/__nSamplingEvents__/{params.sampling_depth}/g" \
        -e "s/__date__/$(date +'%Y-%m-%d')/g" \
        -e "s/__functionsScript__/{params.script_dir}/g" \
        -e "s/__author__/{params.author}/g" \
        -e "s/__input_dir__/{params.input_dir}/g" \
        -e "s/__sim_dir__/{params.simulation_dir}/g" \
        {input} > {output} ) &> {log}
        """


rule render_markdown_file:
    input:
        PROJECT_DIR / "data" / "markdowns" / "{SAMPLE}_treeSampling.Rmd",
    output:
        PROJECT_DIR / "data" / "htmls" / "{SAMPLE}_treeSampling.html",
    resources:
        mem_mb=2048,
        runtime=2880,
    threads: 16
    conda:
        Path(workflow.basedir) / "envs" / "R.yml"
    log:
        PROJECT_DIR / "logs" / "render_markdown_file.{SAMPLE}.log",
    shell:
        """
        ( Rscript -e "rmarkdown::render('{input}', output_file = '{output}')" ) &> {log}
        """


rule render_topSeparators_markdowns:
    input:
        markdown = PROJECT_DIR / "data" / "markdowns" / "{SAMPLE}_topSeparators.Rmd",
        coverageScore = INPUT_FOLDER / "{SAMPLE}" / "{SAMPLE}_covScore.txt",
        clustersplits = INPUT_FOLDER / "{SAMPLE}" / "{SAMPLE}_postSampling_splits.txt",
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
