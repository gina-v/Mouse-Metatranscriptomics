SAMPLES_STM = ["278-1", "278-2", "278-3", "278-4", "278-5", "279-1", "279-2", "279-3", "279-4", "279-5"]
SAMPLES_BTH = ["278-1", "278-2", "278-3", "278-4", "278-5", "275-1", "275-2", "275-3", "275-4", "275-5"]
MOUSE_INDEX = "reference_genomes/GRCm39/GRCm39"
MICROBES = ["stm", "bth"]
ALL_SAMPLES = ["275-1", "275-2", "275-3", "275-4", "275-5", "278-1", "278-2", "278-3", "278-4", "278-5", "279-1", "279-2", "279-3", "279-4", "279-5"]

rule all:
  input:
    # expand("raw_rna_seq/trim/{sample}_bt_trim.fastq", sample=SAMPLES),
    # expand("mouse_filtered/{sample}_{microbes}_filtered.fq.gz", sample = SAMPLES_STM, microbes = "stm"),
    # expand("mouse_filtered/{sample}_{microbes}_filtered.fq.gz", sample = SAMPLES_BTH, microbes = "bth"),
    # expand("sam/{sample}_{microbes}.sam", sample = SAMPLES_STM, microbes = "stm"),
    # expand("sam/{sample}_{microbes}.sam", sample = SAMPLES_BTH, microbes = "bth"),
    expand("sam/sam_11_09/{sample}_{microbes}.sorted.bam", sample = SAMPLES_BTH, microbes = "bth"),
    expand("sam/sam_11_09/{sample}_{microbes}.sorted.bam", sample = SAMPLES_STM, microbes = "stm"),
   # expand("feature_counts/{microbes}_counts.txt", microbes = MICROBES),
   # expand("feature_counts/{microbes}_counts.txt", microbes = MICROBES),
rule fastp:
  input:
    fq = "raw_rna_seq/{microbes}_raw_seq/{sample}.fastq.gz"
  output: 
    trim = "raw_rna_seq/trim_11_09/{sample}_{microbes}_trim.fastq",
    json = "raw_rna_seq/reports_11_09/{sample}_{microbes}_trim.json", 
    html = "raw_rna_seq/reports_11_09/{sample}_{microbes}_trim.html",
  shell:
    """
    fastp --in1 {input.fq} \
    --out1 {output.trim} \
    --detect_adapter_for_pe \
    --qualified_quality_phred 4 \
    --length_required 31 --correction \
    --json {output.json} \
    --html {output.html}
    """
rule filter_mouse_reads:
  input:
    trim = "raw_rna_seq/trim_11_09/{sample}_{microbes}_trim.fastq"
  output:
    filt = "mouse_filtered/filt_11_09/{sample}_{microbes}_filtered.fq.gz"
  shell:
    """
    bowtie2 \
    --threads 32 \
    -x {MOUSE_INDEX} \
    -U {input.trim} \
    --un-gz {output.filt} > /dev/null
    """

rule align_to_microbe:
  input:
    filt = "mouse_filtered/filt_11_09/{sample}_{microbes}_filtered.fq.gz"
  output:
    sam = "sam/sam_11_09/{sample}_{microbes}.sam"
  threads: 32
  shell:
    """
    bowtie2 --threads {threads} \
    -x reference_genomes/{wildcards.microbes}/{wildcards.microbes} \
    -U {input.filt} -S {output.sam}
    """    

rule convert_sam_to_sorted_bam:
  input:
    sam = "sam/sam_11_09/{sample}_{microbes}.sam"
  output:
    bam = "sam/sam_11_09/{sample}_{microbes}.sorted.bam"
  threads: 16
  shell:
    """
    samtools sort -o {output.bam} -O BAM -@ {threads} {input.sam}
    """

# featureCounts -a reference_genomes/bth/bth_genome.gtf -o feature_counts/bth_1110.txt sam/sam_11_09/*bth.sorted.bam -t transcript --extraAttributes gene_name
 
#featureCounts -a reference_genomes/stm/stm_genome_temp.gtf -o feature_counts/stm_1113.txt sam/sam_11_09/*stm.sorted.bam -t transcript --extraAttributres gene_name

