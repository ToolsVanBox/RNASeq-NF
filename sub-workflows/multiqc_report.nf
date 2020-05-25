include MultiQC from '../NextflowModules/MultiQC/1.8/MultiQC.nf' params( optional:params.options.MultiQC )
include customQC from '../utils/customQC.nf' params(params)

workflow multiqc_report {
    take:
      title
      fastqc_logs
      trim_logs
      sortmerna_logs
      star_logs
      post_mapping_qc_logs
      fc_logs
      salmon_logs
      
    main: 
     if (params.rmd_template) {
          rmd_file = Channel
              .fromPath(params.rmd_template, checkIfExists: true)
              .ifEmpty { exit 1, "QC rmd template not found: ${params.rmd_template}"}
      }
      qc_files = Channel.empty().mix( fastqc_logs, trim_logs, sortmerna_logs, star_logs, post_mapping_qc_logs, fc_logs, salmon_logs ).collect()
      MultiQC( qc_files.map { qc_files -> [title, qc_files] } )
      customQC(title, rmd_file, qc_files )
}
