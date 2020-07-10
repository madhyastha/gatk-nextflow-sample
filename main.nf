params.input = "s3://aws-batch-genomics-shared/secondary-analysis/example-files/fastq"
params.reference = "s3://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
params.dict = "s3://broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
params.sample_id = "NIST7035"
params.chromosomes = "chr21"


params.dbsnp = "s3://broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf"
params.golden_indel = "s3://broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz"
params.golden_indel_index = "s3://broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz_tbi"

params.hapmap = "s3://broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz"

// this is used as the publishDir in a couple processes.
// users need to specify a bucket that they have write access to for outputs
// otherwise you will get Access Denied errors that will end up terminating the workflow
params.output = "NONE"


// nextflow script is based on Groovy, so all language constructs therein
// can be used in workflow definitions.
// here is a mapping to define the specific container versions for the tools
// in the pipeline.
def containers = [
  bwa: "biocontainers/bwa:v0.7.15_cv4",
  samtools: "biocontainers/samtools:v1.7.0_cv4",
  bcftools: "biocontainers/bcftools:v1.5_cv3",
  gatk: "broadinstitute/gatk:4.1.3.0"
]

sample_id = params.sample_id
output_dir = "${params.output}/${params.sample_id}"


ref_name = file(params.reference).name

golden_indel = file(params.golden_indel).name
dbsnp = file(params.dbsnp).name


ref_dbsnp = Channel
  .fromPath("${params.dbsnp}*")
  .toList()

ref_golden_indel = Channel
  .fromPath("${params.golden_indel}*")
  .toList()


ref_indices = Channel
  .fromPath("${params.reference}*")
  .toList()

ref_dict = Channel
  .fromPath("${params.dict}")
  .toList()

reads = Channel
  .fromPath("${params.input}/${sample_id}_*{1,2}*{fastq.gz}")
  .toList()


chromosomes = Channel.fromList( params.chromosomes.tokenize(',') )
chromosomes.into { chromosomes_mpileup; chromosomes_call }

log.info """
script: ${workflow.scriptId}
session: ${workflow.sessionId}
sample-id: ${sample_id}
"""

process bwa_mem {
    container "${containers.bwa}"
    cpus 8
    memory "64 GB"

  input:
    file '*' from ref_indices
    file '*' from reads
  
  output:
    file "${sample_id}.sam" into sam_file
  
  script:
  """
  bwa mem -M -t 16 -p -R '@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:Illumina' \\
        ${ref_name} \
        ${sample_id}_*1*.fastq.gz \
        > ${sample_id}.sam
  """
}


process samtools_sort {
    container "${containers.samtools}"
    cpus 8
    memory "32 GB"

    if (params.output != 'NONE') {
      publishDir "${output_dir}", enabled: params.output != 'NONE'
    }

  input:
    file "${sample_id}.sam" from sam_file
  
  output:
    file "${sample_id}.bam" into bam_file
  
  script:
  """
  samtools sort \
        -@ 16 -O BAM \
        -o ${sample_id}.bam \
        ${sample_id}.sam
  """
}

process samtools_index {
    container "${containers.samtools}"
    cpus 8
    memory "32 GB"

    if (params.output != 'NONE') {
      publishDir "${output_dir}", enabled: params.output != 'NONE'
    }

  input:
    file "${sample_id}.bam" from bam_file
  
  output:
    file "${sample_id}.bam.bai" into bai_file
  
  script:
  """
  samtools index \
        ${sample_id}.bam
  """
}

process MarkDuplicates {
    container "${containers.gatk}"
    cpus 4
    memory "32 GB"
    
  input:
    file "*" from ref_indices
    file "${sample_id}.bam" from bam_file
    file "${sample_id}.bai" from bai_file

    output:
    file "${sample_id}_MarkDup.bam" into bam_markdup

    script:
    """
    gatk MarkDuplicates -I ${sample_id}.bam -M metrics.txt -O ${sample_id}_MarkDup.bam  
    """
}


process BaseRecalibrator {
    container "${containers.gatk}"
    cpus 4
    memory "32 GB"
    
    input:
    file "${sample_id}_MarkDup.bam" from bam_markdup
    file '*' from ref_indices
    file '*' from ref_dict
    file '*' from ref_golden_indel
    file '*' from ref_dbsnp

    output:
    file "${sample_id}_recal_data.table" into BaseRecalibrator_table

    script:
    """
    gatk BaseRecalibrator \
    -I ${sample_id}_MarkDup.bam \
    --known-sites $dbsnp \
    --known-sites $golden_indel \
    -O ${sample_id}_recal_data.table \
    -R ${ref_name}
    """
}


process ApplyBQSR {
    container "${containers.gatk}"
    cpus 4
    memory "32 GB"
    
    input:
    file "${sample_id}_MarkDup.bam" from bam_markdup
    file "${sample_id}_recal_data.table" from BaseRecalibrator_table

    output:
    file "${sample_id}_aln-pe_bqsr.bam" into bam_bqsr

    script:
    """
    gatk ApplyBQSR -I ${sample_id}_MarkDup.bam -bqsr ${sample_id}_recal_data.table -O ${sample_id}_aln-pe_bqsr.bam
    """
}



process HaplotypeCaller {
    container "${containers.gatk}"
    cpus 4
    memory "32 GB"
    
    input:
    file "${sample_id}_aln-pe_bqsr.bam" from bam_bqsr
    file '*' from ref_indices
    file '*' from ref_dict

    output:
    file "${sample_id}_haplotypecaller.g.vcf" into haplotypecaller_gvcf
    
    script:
    """
    gatk HaplotypeCaller -I ${sample_id}_aln-pe_bqsr.bam -O ${sample_id}_haplotypecaller.g.vcf --emit-ref-confidence GVCF -R ${ref_name}
    """
}


process GenotypeGVCFs {
    container "${containers.gatk}"
    cpus 4
    memory "32 GB"
        
    if (params.output != 'NONE') {
      publishDir "${output_dir}", enabled: params.output != 'NONE'
    }

    input:
    file "${sample_id}_haplotypecaller.g.vcf" from haplotypecaller_gvcf
    file '*' from ref_indices
    file '*' from ref_dict
    file '*' from ref_dbsnp

    output:
    file "${sample_id}.vcf" 
    

    script:
    """
    gatk GenotypeGVCFs  --dbsnp $dbsnp --variant ${sample_id}_haplotypecaller.g.vcf -R ${ref_name} -O ${sample_id}.vcf
    """
}


workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
