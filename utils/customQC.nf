process customQC {
    tag "customQC ${run_name}" 
    label 'customQC'
    shell = ['/bin/bash', '-euo', 'pipefail']
    container = '/hpc/local/CentOS7/cog_bioinf/singularity_cache/tidyverse_latest.sif'

    input:
	val(run_name)
        file(rmd_file)
        file(qc_files:"*")

    output:
        path("*", emit: qc_report)

    script:
        """
        QC_RNASeq.R --input \$PWD --threads ${task.cpus} --verbose
        """
}
