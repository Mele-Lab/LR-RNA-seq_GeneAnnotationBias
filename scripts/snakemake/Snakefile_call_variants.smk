import pandas as pd

p = os.getcwd()
sys.path.append(p)

from utils import *

include: 'smkmod/blast.smk'
include: 'smkmod/samtools.smk'
include: 'smkmod/formatting.smk'
include: 'smkmod/minimap.smk'
include: 'smkmod/splitting.smk'
include: 'smkmod/qc.smk'
include: 'smkmod/clean.smk'
include: 'smkmod/variant_calling.smk'

configfile: 'snakemake/config.yml'
#config_tsv = 'snakemake/test_data_2.tsv'
config_tsv = '/gpfs/projects/bsc83/Projects/pantranscriptome/pclavell/01_basecalling/data/data_array.tsv'

# TODO map reads to the genome (NOT transcriptome)
use rule minimap_fq as vc_map_genome with:
    input:
        junc_bed = config['ref']['junc_bed'],
        fq = # TODO fill in fastq,
        ref_fa = config['ref']['fa_full']
    output:
        sam = temporary(config['data']['vc']['minimap']['sam'])

# filter reads for primary mappings
use rule sam_filt_for_primary as vc_filt_mappings with:
    input:
        align = rules.vc_map_genome.output.sam
    output:
        align = temporary(config['data']['vc']['minimap']['sam_filt'])

# add rg tag because gatk throws an absolute fit if not
use rule add_rg as vc_add_rg with:
    input:
        align = rules.vc_filt_mappings.output.align
    output:
        align = temporary(config['data']['vc']['minimap']['bam_rg'])

# TODO actually write this rule, there's nothing underneath ATM
# split reads based on where introns are because we have spliced data
use rule split_spliced_reads as vc_split_spliced_reads with:
    input:
        align = rules.vc_add_rg.output.align
    output:
        align = temporary(config['data']['vc']['minimap']['bam_split'])

use rule sort_bam as vc_sort_bam with:
    input:
        bam = rules.vc_split_spliced_reads.output.align
    output:
        bam = config['data']['vc']['minimap']['bam_sort']

use rule index_bam as vc_ind_bam with:
    input:
        bam = rules.vc_sort_bam.output.bam
    output:
        ind = config['data']['vc']['minimap']['bam_sort_ind']

###########################################
########### GATK
###########################################
use rule ref_make_gatk_dict as vc_make_dict with:
    input:
        fa = config['ref']['fa_full']
    output:
        fa_dict = config['ref']['fa_dict']

use rule ref_make_fa_ind as vc_fa_ind with:
    input:
        fa = rules.vc_make_dict.input.fa
    output:
        fai = config['ref']['fa_full_ind']

# actual variant calling
use rule call_variants_gatk as vc_call_variants with:
    input:
        bam = rules.vc_ind_bam.input.bam,
        bai = rules.vc_ind_bam.output.ind,
        fa = rules.vc_make_dict.input.fa,
        ind = rules.vc_fa_ind.output.fai,
        dict = rules.vc_make_dict.output.fa_dict
    output:
        vcf = temporary(config['data']['vc']['gatk']['vcf_gz']),
        bam = temporary(config['data']['vc']['gatk']['vcf_bam'])

use rule gunzip as vc_gunzip_vcf with:
    input:
        ifile = rules.vc_call_variants.output.vcf
    output:
        ofile = config['data']['vc']['gatk']['vcf']

# normalize variants so they're comparable across samples
use rule vcf_norm as vc_vcf_norm with:
    input:
        vcf = rules.vc_gunzip_vcf.output.ofile,
        fa = rules.vc_call_variants.input.fa
    output:
        vcf = config['data']['vc']['gatk']['vcf_norm']
