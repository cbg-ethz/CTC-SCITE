rule prepare_markdown_file:
    input:
        PROJECT_DIR / 'workflow' / 'resources' / 'template.Rmd'
    output:
        PROJECT_DIR / 'data' / 'markdowns' / '{SAMPLE}.Rmd',
    resources:
        mem_mb = 1024,
        runtime = 10,
    params:
        sampling_depth = 1000,
        script_dir =  {markdown_helper_functions},
        author = config['author'],
        input_dir = input_folder,
        simulation_dir = simulations_folder,
    log:
        PROJECT_DIR / 'logs' / 'prepare_markdown_file.{SAMPLE}.log',
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
        PROJECT_DIR / 'data' / 'markdowns' / '{SAMPLE}.Rmd',
    output:
        PROJECT_DIR / 'data' / 'htmls' / '{SAMPLE}.html',
    threads: 16
    conda:
        Path(workflow.basedir) / 'envs' / 'R.yml',
    log:
        PROJECT_DIR / 'logs' / 'render_markdown_file.{SAMPLE}.log',
    shell:
        """
        ( Rscript -e "rmarkdown::render('{input}', output_file = '{output}')" ) &> {log}
        """

rule render_topSeparators_markdowns:
    input:
        PROJECT_DIR / 'data' / 'markdowns' / '{SAMPLE}_top_Separators.Rmd',
    output:
        PROJECT_DIR / 'data' / 'htmls' / '{SAMPLE}_top_Separators.html',
    threads: 4
    log:
        PROJECT_DIR / 'logs' / 'render_topSeparators_markdowns.{SAMPLE}.log',
    shell:
        """
        ( Rscript -e "rmarkdown::render('{input}', output_file = '{output}')" ) &> {log}
        """

        