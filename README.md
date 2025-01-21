# Ensemble_calling_cnv 
Here's an overview of how to get started with a Nextflow workflow that implements an ensemble calling, annotation, and visualization pipeline. This pipeline will process genomic data, including the steps of variant calling, merging multiple callsets, annotating the variants, and visualizing the results. 
## 1.1 Install Nextflow and Prerequisites 
Before starting, ensure you have **Nextflow** and all the necessary tools installed. You will also need to set up a computational environment that supports Nextflow execution. 
## 1.2 Configuring the nextflow.config File 

The **nextflow.config** file is crucial for setting up parameters, resource allocations, and execution profiles. Below is an outline of the essential configurations and their meanings. 
# 2. Configuration of nextflow.config 
This file allows you to define paths, memory allocations, environment variables, and execution profiles for running the workflow. 
## 2.1 Configuration Parameters 

In the `params` section of the **nextflow.config** file, you will define the following: 

  

params { 

      // Path to the reference genome file (FASTA format) 
    genome = "/path/to/reference/genome.fa"   
    // Path to the directory containing BAM files 
    input = "/path/to/bam_files"   
    // List of directories containing .VCF.GZ and .VCF.GZ.TBI files 
    directories = ["/path/to/vcf_directory1", "/path/to/vcf_directory2"]   
    // Pattern for selecting BAM files from the input directory 
    samples = "*.bam*"   
    // Pattern for selecting VCF files from the directories 
    ensamples = "*.vcf*"   

    // Memory configuration for each tool 

    clu_mem = '16GB'  // Memory for SurVClusterer 
    sur_mem = '24GB'  // Memory for SURVIVOR 
    tru_mem = '32GB'  // Memory for Truvari 
    annotsv_mem = '20GB'  // Memory for AnnotSV annotation 
    knotsv_mem = '20GB'  // Memory for KnotSV annotation (HTML and XLS) 
    vep_mem = '20GB'  // Memory for VEP annotation 

} 

2.1.1 Memory Parameters: 
- The memory (sur_mem, tru_mem, etc.) for each tool is adjusted to ensure the processes have enough resources based on your input files and the complexity of the tasks. 
Added annotsv_mem, knotsv_mem, and vep_mem for memory requirements for the AnnotSV, KnotSV, and VEP annotation processes, respectively. 

  

2.1.2 Input and Output File Patterns: 
- The input and samples parameters are set for the BAM files to match your workflow setup. 
- The directories and ensamples parameters are for selecting VCF files from different directories, ensuring flexibility in file selection. 

  

## 2.2 Environment Variables and Execution Profiles 
You will define environment variables, such as paths to libraries and executables, and specify execution profiles for different environments (SLURM or Local). 

  
  env { 

    LIBRARY_PATH = "/path/to/library"   
    PATH = "/path/to/bin:$PATH"   
    JAVA_HOME = "/path/to/java"   

} 

  

profiles { 

    // SLURM profile (for cluster execution) 

    slurm { 

        executor.name = 'slurm' 
        queue = 'batch' 
        memory = '32GB' 
        cpus = 4 
    } 

  // Local profile (for local execution) 

    local { 

        executor.name = 'local'
        memory = '16GB' 
        cpus = 2 

    } 

} 

  

# 3. Workflow Execution 
Once you’ve configured your **nextflow.config**, you can begin executing your workflow. 

## 3.1 Create Directories and Set Up Environment 
Before running the workflow, set up the necessary directories: 
```
mkdir nextflow_running 
cd nextflow_running 
```
## 3.2 Running the Workflow on SLURM Cluster 
To run the workflow on a **SLURM**-based cluster, use the following command: 
```
nextflow run../ Script.nf -profile slurm 
```

## 3.3 Running the Workflow Locally 
For local execution: 
```
nextflow run ../Script.nf -profile local 
```
## 3.4 Resuming the Workflow 
To resume the workflow from a previous point in case of failure: 
```
nextflow run ../Script.nf -profile slurm -resume 
```
# 4. Workflow Processes and Parameters 
Each process in the workflow is associated with specific tools and configurations for variant calling, merging, annotation, and visualization. 

## 4.1 survclusterer_ensemble (SurVClusterer) 
This process clusters structural variants using **SurVClusterer**, which groups similar variants into meaningful clusters. 

```
clusterer -d 1000 -t 4 file.txt ${params.genome} -o output_file.vcf
```
Parameters: 
- -d 1000    Clusters variants within 1000 base pairs. 
- -t 4       Uses 4 threads for parallel processing. 

## 4.2 survivor_ensemble (SURVIVOR) 
**SURVIVOR** merges multiple VCF files into a single output and filters variants based on various parameters. 
```
SURVIVOR merge file.txt 1000 1 1 -1 -1 -1 output_file.vcf 
```
 
Parameters:
- file.txt: Paths to the input VCF files for merging.
- 1000: Variants within 1000 base pairs are merged.
- 1: Minimum number of supporting samples required.
- 1: Includes structural variant type in merging (1 = Yes, 0 = No).
- -1: Excludes strand information from merging.
- -1: Excludes sequence information from merging.
- -1: Uses default distance for merging breakpoints.
- output_file.vcf: Output file containing merged variants.

Explanation: This command ensures that only variants with sufficient support are retained and combines multiple VCF files. 

## 4.3 truvari_ensemble (Truvari) 
**Truvari** collapses and merges structural variants into a concise VCF representation. 
```
truvari collapse -i bcftools.merge.vcf.gz -o truvari.vcf.gz -c output_file.vcf.gz 
```
Parameters:
- -i bcftools.merge.vcf.gz: Input VCF file containing structural variants to be collapsed.
- -o truvari.vcf.gz: Output VCF file with collapsed and filtered variants.
- -c output_file.vcf.gz: Comparison file for merging and filtering variants.

Explanation: **Truvari** reduces redundancy and improves structural variant analysis by merging similar variants. 

  
## 4.4 annotsv_annotation (AnnotSV) 
**AnnotSV** annotates structural variants with functional information, helping researchers understand their biological relevance. 
```
${ANNOTSV_DIR} -genomeBuild GRCh38 -SVinputFile truvari.vcf -outputFile annotsv -outputDir AnnotSV_output -SVinputInfo 1 
```
Parameters:
- ${ANNOTSV_DIR}: Path to the AnnotSV executable.
- genomeBuild GRCh38: Genome build version (e.g., GRCh38).
- SVinputFile truvari.vcf: Input VCF file containing SVs for annotation.
- outputFile annotsv: Base name for the annotated output file.
- outputDir AnnotSV_output: Directory for storing output files.
- SVinputInfo 1: Includes additional SV info from the input file.

Explanation: **AnnotSV** adds biological context, such as gene impact and disease associations, to the structural variants. 

## 4.5 knotsv_annotation_html and knotsv_annotation_xslm (KnotSV) 
**KnotSV** generates HTML and Excel reports with visualizations of annotated structural variants. 
```
perl /home/ale_sab/knotAnnotSV/knotAnnotSV2XL.pl --configFile /home/ale_sab/knotAnnotSV/config_AnnotSV.yaml --annotSVfile ${annot_out} --outDir KnotSV_dir --outPrefix knotsv_output --genomeBuild GRCh38 --LOEUFcolorRange 1 --geneCountThreshold 40 
```
Parameters: 

- -configFile /home/ale_sab/knotAnnotSV/config_AnnotSV.yaml: Path to the configuration file for knotAnnotSV.
- -annotSVfile ${annot_out}: Input AnnotSV file to process.
- -outDir KnotSV_dir: Directory for storing the output Excel file.
- -outPrefix knotsv_output: Prefix for the output file name.
- -genomeBuild GRCh38: Specifies the genome build (e.g., GRCh38).
- -LOEUFcolorRange 1: Enables LOEUF color-coding for enhanced visualization.
- -geneCountThreshold 40: Filters and highlights genes with a count above the threshold.

Explanation: **KnotSV** outputs detailed visual reports that aid in the interpretation of structural variants. 

## 4.6 vep_annotation (VEP) 
**VEP** annotates VCF files with detailed gene information, including predicted functional consequences of variants. 
```
singularity exec /home/ale_sab/project_cnv_vr/ensembl-vep_latest.sif vep -i cleaned_input.vcf --format vcf --output_file ${vep_annotated_vcf} --vcf --everything --assembly GRCh38 --symbol --canonical --vcf_info_field ANN --cache --offline --dir_cache ~/.vep --force_overwrite --cache_version 113 --sift b --polyphen b --biotype --hgvs --fasta 
```
 Parameters: 
- --singularity exec: Runs VEP within a Singularity container.
- --/home/ale_sab/project_cnv_vr/ensembl-vep_latest.sif: Path to the Singularity container image for VEP.
- --i cleaned_input.vcf: Input VCF file containing variants to annotate.
- --format vcf: Specifies that the input file is in VCF format.
- --output_file ${vep_annotated_vcf}: Name of the annotated output VCF file.
- --vcf: Ensures the output remains in VCF format.
- --everything: Enables all available annotations.
- --assembly GRCh38: Genome assembly version (GRCh38).
- --symbol: Includes gene symbols in the annotations.
- --canonical: Flags canonical transcripts.
- --vcf_info_field ANN: Adds annotations to the VCF INFO field using the ANN tag.
- --cache --offline: Uses local annotation cache for offline operation.
- --dir_cache ~/.vep: Path to the local VEP cache directory.
- --force_overwrite: Overwrites existing output files.
- --cache_version 113: Specifies the Ensembl cache version (e.g., v113).
- --sift b: Includes SIFT (Sorting Intolerant From Tolerant) predictions with scores.
- --polyphen b: Includes PolyPhen (Polymorphism Phenotyping) predictions with scores.
- --biotype: Adds biotype information for variants.
- --hgvs: Includes HGVS (Human Genome Variation Society) notations.
- --fasta: Uses a local FASTA file for sequence lookup
  
Explanation: VEP enriches input variants with functional annotations, pathogenicity predictions, and transcript details for downstream interpretation.
  

