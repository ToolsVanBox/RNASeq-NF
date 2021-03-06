### 2.2 Configuration files

**Analysis configuration**

Default [settings](./settings.md) are stored in `./conf/base.config` and will be automatically set upon execution of the pipeline. These defaults can be overwriten by either appending the setting of interest directly as command-line arguments or storing them in a seperate configuration file. An example of such a configuration file can be found in  `./conf/test_run.config`. 

Runtime specific resources (memory, cpu's) should be sufficent for most jobs but can be always be altered if required. 

```
  withLabel : SortMeRNA_4_2_0 {
      time = '24h'
      penv = 'threaded'
      cpus = 4
      memory = '15G'
      publishDir.path = "${params.out_dir}/QC/"
      publishDir.mode = 'copy'
      publishDir.saveAs = {filename ->
                     if (filename.indexOf("_rRNA_report.txt") > 0) "SorteMeRNA/$filename"
                     else if (filename.indexOf("_filtered_rRNA.fastq.gz") > 0) "SorteMeRNA/rRNA-reads/$filename"
                     else null }

  }
```

The new configuration can be appended with the  `-c` option. For example;

*  `./nextflow run ./RNASeq-NF/main.nf -c ./RNASeq-NF/conf/my-run.config ....` (configuration file)
*  `./nextflow run ./RNASeq-NF/main.nf --skipCount=true.....`(without configuration file)

**Resource configuration**

Genomic resources files should always be included as a separate configuration file via the the `--genome_config` parameter, either directly or within a configuration file as described above. 

In it's most basic form, the pipeline requires the following resource parameters.

* `--genome_fasta` - path to genome sequence (.fasta). The index (.fai) should be stored in the same directory.
* `--genome_gtf` - path to genome annotation (.gtf)

Transcript quantification with Salmon requires transcript sequences in .fasta format (In case no Salmon index is provided)
* `--transcripts_fasta` - path to transcript sequences (.fasta).

GATK requires a sequence dictionary (.dict) and fasta index (.fai) to be present in the same directionary as the reference sequence. Please prepare them before running the pipeline if applicable as described in step 2 and 3 of the [GATK reference preparation guide](https://gatkforums.broadinstitute.org/gatk/discussion/2798/howto-prepare-a-reference-for-use-with-bwa-and-gatk)

Known variants in vcf format (optional) should be included in the genome config as shown below.

```
 genome_known_sites = ['/hpc/cog_bioinf/common_dbs/GATK_bundle/1000G_phase1.indels.b37.vcf',
  '/hpc/cog_bioinf/common_dbs/GATK_bundle/dbsnp_137.b37.vcf',
  '/hpc/cog_bioinf/common_dbs/GATK_bundle/Mills_and_1000G_gold_standard.indels.b37.vcf']
```

All remaining files will be created at run-time and stored in the output folder (reference_resources). Since buidling the STAR index requires a significant amount of memory, ensure that sufficient memory is available on the host system and allocated in process.config. 

**Note** 

***Generating resource files at run-time should only be done upon running the pipeline for the first time with a specific genome build. Include the resource files in your genome config for subsequent analysis to avoid re-building them and putting a strain on the host-system. Example genome configuration files can be found in `./resources/test_run`.*** 



**Nextflow configuration**

Nextflow's base configuration (executors, containerization etc.) are stored in the nextflow.config file.

Set the cacheDir to a location with enough space to store the singularity images. 

`cacheDir = '/path/to/my/cache/dir`


