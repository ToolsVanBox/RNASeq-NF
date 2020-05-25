process customQC {
    tag "customQC ${run_name}" 
    label 'customQC'
    shell = ['/bin/bash', '-euo', 'pipefail']
    //container = 'library://icaoberg/default/r-base:latest'

    input:
	val(run_name)
        file(rmd_file)
        file(qc_files:"*")

    output:
        path("*", emit: qc_report)

    script:
        """
        module load R/3.5.1
        Rscript -e "rmarkdown::render('QC.Rmd')" QC_RNASeq.R --input \$PWD --threads ${task.cpus} --verbose
        """
}
