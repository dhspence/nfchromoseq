# dhslab/nfchromoseq pipeline parameters

Nextflow ChromoSeq tumor-only WGS pipeline

## Dragen align inputs

Required input files for ChromoSeq dragen alignments.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `dragen_staging_path` | Staging path for dragen intermediate files | `string` | /staging/intermediate-results-dir |  |  |

## Input/output options

Define where the pipeline should find input data and save output data.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `input` | Path to Mastersheet. <details><summary>Help</summary><small>You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row.</small></details>| `string` |  | True |  |
| `outdir` | Batch directory for output. | `string` |  | True |  |
| `rundir` | Illumina run folder path. | `string` |  |  |  |
| `demuxdir` | Location of demux fastq files. | `string` | null/demux_fastq |  |  |
| `realign` | Option to realign provided cram files before analysis. N.B. crams are reverted and aligned automatically if other inputs are provided for that sample. | `boolean` |  |  |  |
| `publish_dragen_outputs` | Option to copy dragen outputs to outdir if analysis is performed. | `boolean` |  |  |  |
| `run_demux` | Run demux (default is true. mainly for debugging) | `boolean` | True |  |  |
| `run_align` | Run align (default is true. mainly for debugging) | `boolean` | True |  |  |
| `run_analysis` | Run analysis (default is true. mainly for debugging) | `boolean` | True |  |  |
| `debug` | Option to print channels for debuggin. | `boolean` |  |  |  |
| `tracedir` | Location of pipeline trace info | `string` | null/pipeline_info |  |  |
| `email` | Email address for completion summary. <details><summary>Help</summary><small>Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.</small></details>| `string` |  |  |  |
| `multiqc_title` | MultiQC report title. Printed as page header, used for filename if not otherwise specified. | `string` | ChromoSeq QC |  |  |

## Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `help` | Display help text. | `boolean` |  |  | True |
| `version` | Display version and exit. | `boolean` |  |  | True |
| `publish_dir_mode` | Method used to save pipeline results to output directory. <details><summary>Help</summary><small>The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.</small></details>| `string` | copy |  | True |
| `email_on_fail` | Email address for completion summary, only when pipeline fails. <details><summary>Help</summary><small>An email address to send a summary email to when the pipeline is completed - ONLY sent if the pipeline does not exit successfully.</small></details>| `string` |  |  | True |
| `plaintext_email` | Send plain-text email instead of HTML. | `boolean` |  |  | True |
| `max_multiqc_email_size` | File size limit when attaching MultiQC reports to summary emails. | `string` | 25.MB |  | True |
| `monochrome_logs` | Do not use coloured log outputs. | `boolean` |  |  | True |
| `hook_url` | Incoming hook URL for messaging service <details><summary>Help</summary><small>Incoming hook URL for messaging service. Currently, MS Teams and Slack are supported.</small></details>| `string` |  |  | True |
| `multiqc_config` | Custom config file to supply to MultiQC. | `string` |  |  | True |
| `multiqc_logo` | Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file | `string` |  |  | True |
| `multiqc_methods_description` | Custom MultiQC yaml file containing HTML including a methods description. | `string` |  |  |  |
| `validate_params` | Boolean whether to validate parameters against the schema at runtime | `boolean` | True |  | True |
| `validationShowHiddenParams` | Show all params when using `--help` <details><summary>Help</summary><small>By default, parameters set as _hidden_ in the schema are not shown on the command line when a user runs with `--help`. Specifying this option will tell the pipeline to show all parameters.</small></details>| `boolean` |  |  | True |
| `validationFailUnrecognisedParams` | Validation of parameters fails when an unrecognised parameter is found. <details><summary>Help</summary><small>By default, when an unrecognised parameter is found, it returns a warinig.</small></details>| `boolean` |  |  | True |
| `validationLenientMode` | Validation of parameters in lenient more. <details><summary>Help</summary><small>Allows string values that are parseable as numbers or booleans. For further information see [JSONSchema docs](https://github.com/everit-org/json-schema#lenient-mode).</small></details>| `boolean` |  |  | True |
| `job_group_name` | job group name for Wash U RIS compute | `string` | /dspencer/chromoseq | True |  |
| `bcl_first_tile` | Option to demux first tile only. | `boolean` |  |  |  |
| `mastersheet_fields` | Mastersheet fields. | `string` | ['lanes', 'id', 'index', 'exceptions', 'mrn', 'accession', 'dob', 'sex', 'fastq_list', 'demux_path', 'dragen_path'] |  | True |
| `schema_ignore_params` | Ignore these parameters in the schema. | `string` | publish_dragen_outputs,bcl_first_tile,run_analysis,mastersheet_fields,job_group_name,queue,user_group,dragen_container,rundir,demuxdir,genomes,chromoseq_parameters,chromoseq_inputs,dragen_inputs,dragen_staging_path,input_cram_reference,realign |  | True |
| `show_hidden_params` |  | `boolean` |  |  | True |

## Reference genome options

Reference genome related files and options required for the workflow.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `fasta` | Path to FASTA genome file. <details><summary>Help</summary><small>This parameter is *mandatory* if `--genome` is not specified. If you don't have a BWA index available this will be generated for you automatically. Combine with `--save_reference` to save BWA index for future runs.</small></details>| `string` |  | True |  |
| `input_cram_reference` | Path to reference fasta for crams. N.B provide if different from genome/fasta provided. | `string` | /storage1/fs1/dspencer/Active/spencerlab/refdata/hg38/sequence/hg38_mgi_patch.fa |  |  |
| `genome` | Name of iGenomes reference. <details><summary>Help</summary><small>If using a reference genome configured in the pipeline using iGenomes, use this parameter to give the ID for the reference. This is then used to build the full paths for all required reference genome files e.g. `--genome GRCh38`. <br><br>See the [nf-core website docs](https://nf-co.re/usage/reference_genomes) for more details.</small></details>| `string` |  |  | True |
| `igenomes_base` | Location of S3 bucket for iGenomes. | `string` | s3://ngi-igenomes/igenomes |  | True |
| `genomes` | iGenomes paths. | `string` | [] |  | True |
| `igenomes_ignore` | Do not load the iGenomes reference config. <details><summary>Help</summary><small>Do not load `igenomes.config` when running the pipeline. You may choose this option if you observe clashes between custom parameters and those supplied in `igenomes.config`.</small></details>| `boolean` | True |  | True |

## Institutional config options

Parameters used to describe centralised config profiles. These should not be edited.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `queue` | queue for Wash U RIS compute | `string` | dspencer | True |  |
| `dragen_container` | Dragen docker container (for Wash U RIS compute) | `string` | etycksen/dragen4:4.2.4 |  |  |
| `user_group` | user group for Wash U RIS compute | `string` | compute-dspencer | True |  |
| `custom_config_version` | Git commit id for Institutional configs. | `string` |  |  | True |
| `custom_config_base` | Base directory for Institutional configs. <details><summary>Help</summary><small>If you're running offline, Nextflow will not be able to fetch the institutional config files from the internet. If you don't need them, then this is not a problem. If you do need them, you should download the files from the repo and tell Nextflow where to find them with this parameter.</small></details>| `string` |  |  | True |
| `config_profile_name` | Institutional config name. | `string` |  |  | True |
| `config_profile_description` | Institutional config description. | `string` |  |  | True |
| `config_profile_contact` | Institutional config contact information. | `string` |  |  | True |
| `config_profile_url` | Institutional config URL link. | `string` |  |  | True |

## Max job request options

Set the top limit for requested resources for any single job.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `max_cpus` | Maximum number of CPUs that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the CPU requirement for each process. Should be an integer e.g. `--max_cpus 1`</small></details>| `integer` | 16 |  | True |
| `max_memory` | Maximum amount of memory that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the memory requirement for each process. Should be a string in the format integer-unit e.g. `--max_memory '8.GB'`</small></details>| `string` | 128.GB |  | True |
| `max_time` | Maximum amount of time that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the time requirement for each process. Should be a string in the format integer-unit e.g. `--max_time '2.h'`</small></details>| `string` | 240.h |  | True |
