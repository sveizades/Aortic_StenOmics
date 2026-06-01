# utils.R — helper functions for generating HPC job scripts (Grid Engine SGE) and
# downstream analysis utilities used across the AorticStenOmics pipeline.

# Variables
supplemental_table_output_path <- "../03_output/03_supplemental_tables/"
main_figure_output_path <- "../03_output/01_main_figures/"
supplemental_figure_output_path <- "../03_output/02_supplemental_figures/"
source_data_output_path <- "../03_output/04_source_data/"
UBC_donor_colours <- c(AKY = "#91cd2b", C4J = "#9f40c7", C5H = "#02dd83", CBN = "#b30090", DA5 = "#9ea800", GW1 = "#8984ff", HFU = "#e6a600", I4Q = "#0058b4", INV = "#f17f0b", JI1 = "#0087cc", KS4 = "#b53600", LTO = "#27dec3", M6Y = "#e71b95", MK0 = "#7d7600", O2N = "#ff6ed7", PGS = "#8a4a00", PJN = "#d8a4ff", PUM = "#734a0f", Q8H = "#ff529d", QZJ = "#8a7745", R9L = "#d50035", RW7 = "#75406b", SFA = "#fcb787", UKG = "#a31a1a", V5R = "#ff9ccd", XZD = "#913040", Y84 = "#e89da0", YDA = "#ff6f82")
UBC_disease_colours <- c(control="#4D8FCB", mildmoderate="#E7BF47", severe="#D95E5B")
UBC_disease_colours_lighter <- c(control="#96C7FF", mildmoderate="#FFDD8D", severe="#FFAAA9")
UBC_atlas_colours <- c("VICs"= '#332288', "Myeloid"='#CC6677', 'B cells'='#DDCC77', 'T cells'='#117733', 'NK cells'='#88CCEE', 'Endothelial cells'='#882255', 'Adipocytes'='#44AA99', 'Mast cells'='#999933', "Proliferating" = "#AA4499")
UBC_atlas_colours_vics <- c(Transitional= '#1B9E77', Quiescent='#D95F02', Spongiosa='#7570B3', "FAP+ osteogenic"='#E7298A', Contractile='#0F52BA', "Neural crest-like"='#E6AB02', "IFN-stimulated"='black')
SALTIRE3_donor_colours <- c(ESJ = "#5c3989", I1Z = "#56ae6a", MIR = "#6973d8", OJC = "#ac9c3d", UCC = "#c169b9", UVN = "#ba543d", WQJ = "#6b93d9", XRU = "#b84873")
A1_celltype_cols <- c(
  "Valvular interstitial cell"        = "#332288",
  "Endothelial cell"        = "#882255",
  "Granulocyte"     = "#44AA99",
  "Myeloid"      = "#CC6677",
  "NK cell"   = "#88CCEE",
  "Proliferating"     = "#AA4499",
  "T cell"              = "#117733",
  "B cell"      = "#DDCC77"
)
SALTIRE_disease_colours <- c(
  "Control"    = "#4D8FCB",  
  "Sclerosis"  = "#7FB2D6",  
  "Mild"       = "#F4D06F",  
  "Moderate"   = "#F29E4C", 
  "Severe"     = "#D95E5B" 
)
SALTIRE_disease_colours_lighter <- c(
  "Control"    = "#A1C1E2",
  "Sclerosis"  = "#CFE3F3",
  "Mild"       = "#F1DDA1",
  "Moderate"   = "#F3C6A0",
  "Severe"     = "#EAA8A6"
)
vic_cols_black_background <- c(
  "VIC-IFN" = "#7A6E2F",
  "VIC-spongiosa"  = "#958BC2",
  "VIC-FAP-osteogenic" = "#E7298A",
  "VIC-quiescent" = "#F18D2D",
  "VIC-contractile" = "#0D52BA",
  "VIC-transitional" = "#1C9E76",
  "VIC-neural-crest"  = "#8A5A44",
  "other" = "#2D3132")
spatial_niche_cols <- c(
  "Cap" = "#4477AA",
  "Spongiosa" = "#CCBB44",
  "Peri-nodular" = "#66CC55",
  "Immune" = "#AA3377"
)
A2_celltype_cols <- c(
  # VICs
  "VIC-spongiosa"        = "#958BC2",
  "VIC-quiescent"        = "#F18D2D",
  "VIC-transitional"     = "#1C9E76",
  "VIC-contractile"      = "#0D52BA",
  "VIC-FAP-osteogenic"   = "#E7298A",
  "VIC-neural-crest"     = "#8A5A44",
  "VIC-IFN"              = "#111111",
  
  # Endothelium
  "VEC-endocardial"      = "#FF7F00",
  "VEC-vascular"         = "#00A6D6",
  
  # T / NK
  "T-CD4"                = "#0072B2",
  "T-CD8"                = "#E69F00",
  "T-gd"                 = "#56B4E9",
  "NK"                   = "#CC79A7",
  
  # Cycling
  "Proliferating"        = "#AA4499",
  
  # Myeloid
  "Macrophage"           = "#66C2A5",
  "Dendritic cell"       = "#984EA3",
  "Neutrophil"           = "#A6761D",
  "Monocyte"             = "#E41A1C",
  
  # B / plasma / mast
  "B cell"               = "#FB8072",
  "Plasma cell"          = "#80B1D3",
  "Mast cell"            = "#FDB462"
)
A3_celltype_cols <- c(
  
  # ---------------- VICs ----------------
  "VIC-spongiosa"                  = "#958BC2",
  "VIC-quiescent"                  = "#F18D2D",
  "VIC-transitional"               = "#1C9E76",
  "VIC-contractile"                = "#0D52BA",
  "VIC-FAP-osteogenic"             = "#E7298A",
  "VIC-neural-crest"               = "#8A5A44",
  "VIC-IFN"                        = "#111111",
  "VIC-cycling"                    = "#AA4499",
  
  # ---------------- VEC ----------------
  "VEC-arterial"                   = "#1F78B4",
  "VEC-capillary"                  = "#A6CEE3",
  "VEC-venous"                     = "#33A02C",
  "VEC-endocardial"                = "#FF7F00",
  "VEC-lymphatic"                  = "#6A3D9A",
  "VEC-activated"                  = "#B15928",
  
  # ---------------- T cells ----------------
  "T-CD4-naive"                    = "#66C2A5",
  "T-CD4-EM"                       = "#FC8D62",
  "T-CD4-Th1"                      = "#8DA0CB",
  "T-CD4-Th17"                     = "#E78AC3",
  "T-CD4-Treg"                     = "#A6D854",
  "T-CD4-cyto"                     = "#FFD92F",
  "T-CD8-EM"                       = "#E5C494",
  "T-CD8-cyto"                     = "#B3B3B3",
  "T-CD8-MAIT"                     = "#1B9E77",
  "T-gd-V9"                        = "#D95F02",
  "T-gd-cyto"                      = "#7570B3",
  "T-cycling"                      = "#F781BF",
  
  # ---------------- NK ----------------
  "NK-CD56-dim"                    = "#00BFC4",
  "NK-CD56-bright"                 = "#C77CFF",
  "NK-cycling"                     = "#FF61C3",
  
  # ---------------- Myeloid ----------------
  "Classical monocyte"             = "#E41A1C",
  "Non-classical monocyte"         = "#377EB8",
  "Resident macrophage"            = "#4DAF4A",
  "Resident activated macrophage"  = "#984EA3",
  "IFN-stimulated macrophage"      = "#FF7F00",
  "Scar-associated macrophage"     = "#FFFF33",
  "Osteoclast"                     = "#A65628",
  "cDC1"                           = "#F781BF",
  "cDC2"                           = "#999999",
  "moDC"                           = "#66A61E",
  "pDC"                            = "#E6AB02",
  "Neutrophil"                     = "#A6761D",
  "Myeloid-cycling"                = "#1F78B4",
  
  # ---------------- B / plasma / mast ----------------
  "B-naive"                        = "#8DD3C7",
  "B-memory"                       = "#FB8072",
  "Plasma cell"                    = "#80B1D3",
  "Mast cell"                      = "#FDB462"
)

A3_celltype_cols <- c(
  
  # ---------------- VICs ----------------
  "VIC-spongiosa"                  = "#958BC2",
  "VIC-quiescent"                  = "#F18D2D",
  "VIC-transitional"               = "#1C9E76",
  "VIC-contractile"                = "#0D52BA",
  "VIC-FAP-osteogenic"             = "#E7298A",
  "VIC-neural-crest"               = "#8A5A44",
  "VIC-IFN"                        = "#111111",
  "VIC-cycling"                    = "#AA4499",
  
  # ---------------- VEC ----------------
  "VEC-arterial"                   = "#1F78B4",
  "VEC-capillary"                  = "#A6CEE3",
  "VEC-venous"                     = "#33A02C",
  "VEC-endocardial"                = "#FF7F00",
  "VEC-lymphatic"                  = "#6A3D9A",
  "VEC-activated"                  = "#B15928",
  
  # ---------------- T cells ----------------
  "T-CD4-naive"                    = "#66C2A5",
  "T-CD4-EM"                       = "#FC8D62",
  "T-CD4-Th1"                      = "#8DA0CB",
  "T-CD4-Th17"                     = "#E78AC3",
  "T-CD4-Treg"                     = "#A6D854",
  "T-CD4-cyto"                     = "#FFD92F",
  "T-CD8-EM"                       = "#E5C494",
  "T-CD8-cyto"                     = "#B3B3B3",
  "T-CD8-MAIT"                     = "#1B9E77",
  "T-gd-V9"                        = "#D95F02",
  "T-gd-cyto"                      = "#7570B3",
  "T-cycling"                      = "#F781BF",
  
  # ---------------- NK ----------------
  "NK-CD56-dim"                    = "#00BFC4",
  "NK-CD56-bright"                 = "#C77CFF",
  "NK-cycling"                     = "#FF61C3",
  
  # ---------------- Myeloid ----------------
  "Classical monocyte"             = "#E41A1C",
  "Non-classical monocyte"         = "#377EB8",
  "Resident macrophage"            = "#4DAF4A",
  "Resident activated macrophage"  = "#984EA3",
  "IFN-stimulated macrophage"      = "#FF7F00",
  "Scar-associated macrophage"     = "#FFFF33",
  "Osteoclast"                     = "#A65628",
  "cDC1"                           = "#F781BF",
  "cDC2"                           = "#999999",
  "moDC"                           = "#66A61E",
  "pDC"                            = "#E6AB02",
  "Neutrophil"                     = "#A6761D",
  "Myeloid-cycling"                = "#1F78B4",
  
  # ---------------- B / plasma / mast ----------------
  "B-naive"                        = "#8DD3C7",
  "B-memory"                       = "#FB8072",
  "Plasma cell"                    = "#80B1D3",
  "Mast cell"                      = "#FDB462"
)

A2_celltype_cols_dark_background <- c(
  # VICs
  "VIC-spongiosa"        = "#958BC2",
  "VIC-quiescent"        = "#F18D2D",
  "VIC-transitional"     = "#1C9E76",
  "VIC-contractile"      = "#0D52BA",
  "VIC-FAP-osteogenic"   = "#E7298A",
  "VIC-neural-crest"     = "#8A5A44",
  "VIC-IFN"              = "#7A6E2F",
  
  # Endothelium
  "VEC-endocardial"      = "#FF7F00",
  "VEC-vascular"         = "#00A6D6",
  
  # T / NK
  "T-CD4"                = "#0072B2",
  "T-CD8"                = "#E69F00",
  "T-gd"                 = "#56B4E9",
  "NK"                   = "#CC79A7",
  
  # Cycling
  "Proliferating"        = "#AA4499",
  
  # Myeloid
  "Macrophage"           = "#ffd700",
  "Dendritic cell"       = "#984EA3",
  "Neutrophil"           = "#A6761D",
  "Monocyte"             = "#E41A1C",
  
  # B / plasma / mast
  "B cell"               = "#FB8072",
  "Plasma cell"          = "#80B1D3",
  "Mast cell"            = "#FDB462"
)
VIC_celltype_names <- c("VIC-transitional" = "Transitional", 
                        "VIC-quiescent" = "Quiescent", 
                        "VIC-spongiosa" = "Spongiosa", 
                        "VIC-FAP-osteogenic" = "FAP+ osteogenic", 
                        "VIC-contractile" = "Contractile", 
                        "VIC-neural-crest" = "Neural crest-like", 
                        "VIC-IFN" = "IFN-stimulated")
VIC_celltype_levels <- c("VIC-transitional", 
                        "VIC-quiescent", 
                        "VIC-spongiosa", 
                        "VIC-FAP-osteogenic", 
                        "VIC-contractile", 
                        "VIC-neural-crest", 
                        "VIC-IFN-stim")
VIC_solo_celltype_levels <- c("Transitional", 
                         "Quiescent", 
                         "Spongiosa", 
                         "FAP+ osteogenic", 
                         "Contractile", 
                         "Neural crest-like", 
                         "IFN-stimulated")
# Generates a Grid Engine bash script for QuantSeq alignment:
#   FastQC (pre- and post-trimming) → Cutadapt → STAR → featureCounts.
# Writes <sample_id>_alignment_script.sh to output_path.
# Args:
#   sample_id: unique sample identifier used for job naming and output filenames
#   sample_id_fastq: fastq filename prefix as stored on the NAS
#   output_path: directory in which the generated .sh script will be written
quantseq_alignment <- function(email = Sys.getenv("EDDIE_EMAIL"),
                             star_reference = Sys.getenv("STAR_REFERENCE"), 
                             sample_id, 
                             sample_id_fastq,
                             output_path) {
  quantseq_alignment_script <- paste0("#!/bin/bash

#$ -N ", sample_id, "-alignment
#$ -cwd
#$ -l h_rt=2:00:00
#$ -l h_rss=4G
#$ -pe sharedmem 8
#$ -e ./script_errors/", sample_id, "-alignment_error.txt
#$ -o ./script_output/", sample_id,"-alignment_output.txt
#$ -m beas
#$ -M ", email,"

# Initialising the environment modules
. /etc/profile.d/modules.sh
source ~/.bashrc
module load roslin/fastqc/0.12.1
module load igmm/apps/cutadapt/4.6
module load igmm/apps/STAR/2.7.11b
module load igmm/apps/subread/2.1.0
module load roslin/samtools/1.9

# Download sequencing data from NAS
STARTTIME=$(date +%s)
rsync -az coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "/Data/Sequencing/QuantSeq/Formatted/", sample_id_fastq, ".fastq.gz $TMPDIR
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo NAS download completed in $RUNTIME s

# Run fastq on untrimmed fastq
STARTTIME=$(date +%s)
fastqc_untrimmed_output_folder=$TMPDIR/", sample_id,"_fastqc_untrimmed
mkdir $fastqc_untrimmed_output_folder
fastqc -o $fastqc_untrimmed_output_folder -t 8 $TMPDIR/", sample_id_fastq, ".fastq.gz
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo FastQC completed in $RUNTIME s

# Upload untrimmed FastQC output to NAS
STARTTIME=$(date +%s)
rsync -az $fastqc_untrimmed_output_folder/ coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "Results/01-bulk-seq_Processing/01-FastQC/", sample_id,"
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo FastQC upload completed in $RUNTIME s

# Run cutadapt
cutadapt -m 20 -O 20 -a \"polyA=A{20}\" -a \"QUALITY=G{20}\" -n 2 $TMPDIR/", sample_id_fastq, ".fastq.gz | \\
cutadapt -m 20 -O 3 --nextseq-trim=10 -a \"r1adapter=A{18}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000\" - | \\
cutadapt -m 20 -O 20 -g \"r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20\" --discard-trimmed -o $TMPDIR/", sample_id_fastq, "_trimmed.fastq.gz -

# Upload cutadapt output to NAS
STARTTIME=$(date +%s)
rsync -az $TMPDIR/", sample_id_fastq, "_trimmed.fastq.gz coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "Results/01-bulk-seq_Processing/02-cutadapt/
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo cutadapt upload completed in $RUNTIME s

# Run fastq on trimmed fastq
STARTTIME=$(date +%s)
fastqc_trimmed_output_folder=$TMPDIR/", sample_id,"_fastqc_trimmed
mkdir $fastqc_trimmed_output_folder
fastqc -o $fastqc_trimmed_output_folder -t 8 $TMPDIR/", sample_id_fastq, "_trimmed.fastq.gz
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo FastQC completed in $RUNTIME s

# Upload trimmed FastQC output to NAS
STARTTIME=$(date +%s)
rsync -az $fastqc_trimmed_output_folder/ coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "Results/01-bulk-seq_Processing/01-FastQC/", sample_id,"/
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo FastQC upload completed in $RUNTIME s

# Run STAR
STARTTIME=$(date +%s)
STAR_output_folder=$TMPDIR/", sample_id,"_STAR
mkdir $STAR_output_folder
STAR --runThreadN 8 \\
  --genomeDir ", star_reference, " \\
  --readFilesIn $TMPDIR/", sample_id_fastq, "_trimmed.fastq.gz \\
  --readFilesCommand zcat \\
  --outFilterType BySJout \\
  --outFilterMultimapNmax 20 \\
  --alignSJoverhangMin 8 \\
  --alignSJDBoverhangMin 1 \\
  --outFilterMismatchNmax 999 \\
  --outFilterMismatchNoverLmax 0.1 \\
  --alignIntronMin 20 \\
  --alignIntronMax 1000000 \\
  --alignMatesGapMax 1000000 \\
  --outSAMattributes NH HI NM MD \\
  --outFileNamePrefix $STAR_output_folder/", sample_id,"_ \\
  --outSAMtype BAM SortedByCoordinate \\
  --quantMode GeneCounts
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo STAR completed in $RUNTIME s

# Index bam
samtools index -@ 8 $STAR_output_folder/", sample_id,"_Aligned.sortedByCoord.out.bam

# Upload STAR output to NAS
STARTTIME=$(date +%s)
rsync -az $STAR_output_folder/ coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "Results/01-bulk-seq_Processing/03-STAR/", sample_id,"/
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo STAR upload completed in $RUNTIME s

# Run featureCounts
STARTTIME=$(date +%s)
featureCounts_output_folder=$TMPDIR/", sample_id,"_featureCounts
mkdir $featureCounts_output_folder
$SCRATCH/subread-2.0.2-Linux-x86_64/bin/featureCounts -T 8 -s 1 -t exon -g gene_id -a $SCRATCH/genome/gencode.v44.annotation.gtf -o $featureCounts_output_folder/", sample_id,"_gene_featureCounts_output.txt $STAR_output_folder/", sample_id,"_Aligned.sortedByCoord.out.bam
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo featureCounts completed in $RUNTIME s

# Upload featureCounts output to NAS
STARTTIME=$(date +%s)
rsync -az $featureCounts_output_folder/ coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "Results/01-bulk-seq_Processing/04-featureCounts/", sample_id,"/
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo featureCounts upload completed in $RUNTIME s")
  writeLines(quantseq_alignment_script, paste0(output_path, sample_id, "_alignment_script.sh"))
}

# Generates a Grid Engine bash script for Cell Ranger count (GEX-only libraries):
#   FastQC → cellranger count.
# Writes <pool_id>_1-cellranger_script.sh to output_path.
# Args:
#   fastq_line: comma-separated fastq sample name prefix(es), formatted as
#               <project>_<library_id>_<library_type>_<flow_cell_id>
#               (e.g. CureAS_AO_GEX_BHNJKHDMXY)
#   pool_id: pool identifier used for job naming and output filenames
#   output_path: directory in which the generated .sh script will be written
cellranger_count <- function(email = Sys.getenv("EDDIE_EMAIL"),
                             transcriptome_path = Sys.getenv("CELLRANGER_REFERENCE"), 
                             fastq_line, 
                             pool_id,
                             output_path) {
  cellranger_script <- paste0("#!/bin/bash

#$ -N ", pool_id, "_1-cellranger
#$ -cwd
#$ -l h_rt=32:00:00
#$ -l h_rss=8G
#$ -pe sharedmem 16
#$ -e ./script_errors/", pool_id, "_1-cellranger_error.txt
#$ -o ./script_output/", pool_id,"_1-cellranger_output.txt
#$ -m beas
#$ -M ", email,"

# Initialising the environment modules
. /etc/profile.d/modules.sh
source ~/.bashrc
module load roslin/fastqc/0.12.1
module load igmm/apps/cellranger/10.0.0

# Download sequencing data from NAS
STARTTIME=$(date +%s)
rsync -az coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "/Data/Sequencing/scnRNA-seq/Formatted/", pool_id, " $TMPDIR
input_folder=$TMPDIR/", pool_id,"/", pool_id,"_GEX
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo NAS download completed in $RUNTIME s

# Run fastq on all files of pool ID
STARTTIME=$(date +%s)
fastqc_output_folder=$TMPDIR/", pool_id,"_fastqc
mkdir $fastqc_output_folder
fastqc -o $fastqc_output_folder -t 16 $input_folder/*fastq.gz
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo FastQC completed in $RUNTIME s

# Upload FastQC output to NAS
STARTTIME=$(date +%s)
rsync -az $fastqc_output_folder/ coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "Results/02-scnRNA-seq_Processing/01-FastQC/", pool_id,"
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo FastQC upload completed in $RUNTIME s

# Run cellranger count
STARTTIME=$(date +%s)
cellranger_output_folder=$TMPDIR/", pool_id,"_cellranger
mkdir $cellranger_output_folder
cd $cellranger_output_folder
cellranger count --id=", pool_id," \\
--transcriptome=", transcriptome_path, " \\
--fastqs=$input_folder \\
--sample=", fastq_line, " \\
--create-bam=true \\
--cell-annotation-model=auto \\
--jobmode=local \\
--localcores 16 \\
--localmem 115
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo Cell Ranger completed in $RUNTIME s

# Upload Cell Ranger output to NAS
STARTTIME=$(date +%s)
rsync -az $cellranger_output_folder/ coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "Results/02-scnRNA-seq_Processing/02-CellRanger
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo Cell Ranger upload completed in $RUNTIME s")
  writeLines(cellranger_script, paste0(output_path, pool_id, "_1-cellranger_script.sh"))
}
# Generates a Grid Engine bash script for Cell Ranger multi (GEX + VDJ + hashtag/CSP libraries):
#   FastQC → cellranger multi.
# Writes <pool_id>_1-cellranger_config.csv and <pool_id>_1-cellranger_script.sh to output_path.
# Args:
#   fastq_table: data frame with columns pool_id, library_type (GEX/TCR/BCR/CSP), flow_cell
#   sample_table: data frame with columns sample_id and hashtag_ids for cell hashing demux;
#                 pass "empty" (default) to skip the [samples] section
#   feature_ref: filename (without extension) of the hashtag/antibody feature reference CSV
#   pool_id: pool identifier used for job naming and output filenames
#   output_path: directory in which the generated files will be written
cellranger_multi <- function(email = Sys.getenv("EDDIE_EMAIL"),
                             transcriptome_path = "$GROUP_HOME/SV/references/refdata-gex-GRCh38-2024-A", 
                             vdj_reference_path = "$GROUP_HOME/SV/references/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0",
                             feature_ref_dir = "$GROUP_HOME/SV/references/",
                             feature_ref,
                             fastq_table,
                             sample_table = "empty",
                             pool_id,
                             output_path){
  fastq_table$fastq_id <- paste0(fastq_table$pool_id, "_", fastq_table$library_type, "_", fastq_table$flow_cell)
  fastq_table$fastqs <- paste0("$TMPDIR/", fastq_table$pool_id, "/", fastq_table$pool_id, "_", fastq_table$library_type)
  fastq_table$feature_types <- mapvalues(fastq_table$library_type, from=c("GEX", "CSP", "BCR", "TCR"), to=c("Gene Expression", "Antibody Capture", "VDJ-B", "VDJ-T"))
  cellranger_csv <- paste0("[gene-expression]
reference,", transcriptome_path, "
create-bam,true
cell-annotation-model,auto

[vdj]
reference,", vdj_reference_path, "

[feature]
reference,", feature_ref_dir, feature_ref, ".csv

[libraries]
fastq_id,fastqs,feature_types
", paste0(paste0(fastq_table$fastq_id, ",",fastq_table$fastqs, ",",fastq_table$feature_types), collapse = "\n"))
  if(!is.character(sample_table)){
    cellranger_csv <- paste0(cellranger_csv,"

[samples]
sample_id,hashtag_ids
", paste0(paste0(sample_table$sample_id, ",",sample_table$hashtag_ids), collapse = "\n"))
  }
  writeLines(cellranger_csv, paste0(output_path, pool_id, "_1-cellranger_config.csv"))
  
  cellranger_script <- paste0("#!/bin/bash

#$ -N ", pool_id, "_1-cellranger
#$ -cwd
#$ -l h_rt=36:00:00
#$ -l h_rss=8G
#$ -pe sharedmem 16
#$ -e ./script_errors/", pool_id, "_1-cellranger_error.txt
#$ -o ./script_output/", pool_id,"_1-cellranger_output.txt
#$ -m beas
#$ -M ", email,"

# Initialising the environment modules
. /etc/profile.d/modules.sh
source ~/.bashrc
module load roslin/fastqc/0.12.1

# Download sequencing data from NAS
STARTTIME=$(date +%s)
rsync -az coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "/Data/Sequencing/scnRNA-seq/Formatted/", pool_id, " $TMPDIR
input_folder=$TMPDIR/", pool_id,"/
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo NAS download completed in $RUNTIME s

# Run fastq on all files of pool ID
STARTTIME=$(date +%s)
fastqc_output_folder=$TMPDIR/", pool_id,"_fastqc
mkdir $fastqc_output_folder
fastqc -o $fastqc_output_folder -t 16 $input_folder/", pool_id,"_GEX/*fastq.gz $input_folder/", pool_id,"_TCR/*fastq.gz $input_folder/", pool_id,"_BCR/*fastq.gz $input_folder/", pool_id,"_CSP/*fastq.gz
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo FastQC completed in $RUNTIME s

# Upload FastQC output to NAS
STARTTIME=$(date +%s)
rsync -az $fastqc_output_folder/ coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "Results/02-scnRNA-seq_Processing/01-FastQC/", pool_id,"
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo FastQC upload completed in $RUNTIME s

# Run cellranger multi
STARTTIME=$(date +%s)
sed \"s|\\$TMPDIR|${TMPDIR}|g\" $SCRATCH/scripts/", pool_id, "_1-cellranger_config.csv > $TMPDIR/", pool_id, "_1-cellranger_config.csv
cellranger_output_folder=$TMPDIR/", pool_id,"_cellranger
mkdir $cellranger_output_folder
cd $cellranger_output_folder
", Sys.getenv("CELLRANGER_PATH"), " multi --id=", pool_id," \\
--csv=$TMPDIR/", pool_id, "_1-cellranger_config.csv \\
--jobmode=local \\
--localcores 16 \\
--localmem 115
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo Cell Ranger completed in $RUNTIME s

# Upload Cell Ranger output to NAS
STARTTIME=$(date +%s)
rsync -az $cellranger_output_folder/ coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "Results/02-scnRNA-seq_Processing/02-CellRanger
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo Cell Ranger upload completed in $RUNTIME s")
  writeLines(cellranger_script, paste0(output_path, pool_id, "_1-cellranger_script.sh"))
  
}
# Generates a Grid Engine bash script for Cell Ranger multi for OCM (On-Chip Multiplexing)
#   libraries (GEX + VDJ only, no hashtag/CSP):
#   FastQC → cellranger multi.
# Writes <pool_id>_1-cellranger_config.csv and <pool_id>_1-cellranger_script.sh to output_path.
# Args:
#   fastq_table: data frame with columns pool_id, library_type (GEX/TCR/BCR), flow_cell
#   sample_table: data frame with columns sample_id and hashtag_ids for OCM demux;
#                 pass "empty" (default) to skip the [samples] section
#   pool_id: pool identifier used for job naming and output filenames
#   output_path: directory in which the generated files will be written
cellranger_multi_OCM <- function(email = Sys.getenv("EDDIE_EMAIL"),
                             transcriptome_path = "$GROUP_HOME/SV/references/refdata-gex-GRCh38-2024-A", 
                             vdj_reference_path = "$GROUP_HOME/SV/references/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0",
                             fastq_table,
                             sample_table = "empty",
                             pool_id,
                             output_path){
  fastq_table$fastq_id <- paste0(fastq_table$pool_id, "_", fastq_table$library_type, "_", fastq_table$flow_cell)
  fastq_table$fastqs <- paste0("$TMPDIR/", fastq_table$pool_id, "/", fastq_table$pool_id, "_", fastq_table$library_type)
  fastq_table$feature_types <- mapvalues(fastq_table$library_type, from=c("GEX","BCR", "TCR"), to=c("Gene Expression", "VDJ-B", "VDJ-T"))
  cellranger_csv <- paste0("[gene-expression]
reference,", transcriptome_path, "
create-bam,true
cell-annotation-model,auto

[vdj]
reference,", vdj_reference_path, "

[libraries]
fastq_id,fastqs,feature_types
", paste0(paste0(fastq_table$fastq_id, ",",fastq_table$fastqs, ",",fastq_table$feature_types), collapse = "\n"))
  if(!is.character(sample_table)){
    cellranger_csv <- paste0(cellranger_csv,"

[samples]
sample_id,ocm_barcode_ids
", paste0(paste0(sample_table$sample_id, ",",sample_table$hashtag_ids), collapse = "\n"))
  }
  writeLines(cellranger_csv, paste0(output_path, pool_id, "_1-cellranger_config.csv"))
  
  cellranger_script <- paste0("#!/bin/bash

#$ -N ", pool_id, "_1-cellranger
#$ -cwd
#$ -l h_rt=36:00:00
#$ -l h_rss=8G
#$ -pe sharedmem 16
#$ -e ./script_errors/", pool_id, "_1-cellranger_error.txt
#$ -o ./script_output/", pool_id,"_1-cellranger_output.txt
#$ -m beas
#$ -M ", email,"

# Initialising the environment modules
. /etc/profile.d/modules.sh
source ~/.bashrc
module load roslin/fastqc/0.12.1

# Download sequencing data from NAS
STARTTIME=$(date +%s)
rsync -az coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "/Data/Sequencing/scnRNA-seq/Formatted/", pool_id, " $TMPDIR
input_folder=$TMPDIR/", pool_id,"/
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo NAS download completed in $RUNTIME s

# Run fastq on all files of pool ID
STARTTIME=$(date +%s)
fastqc_output_folder=$TMPDIR/", pool_id,"_fastqc
mkdir $fastqc_output_folder
fastqc -o $fastqc_output_folder -t 16 $input_folder/", pool_id,"_GEX/*fastq.gz $input_folder/", pool_id,"_TCR/*fastq.gz $input_folder/", pool_id,"_BCR/*fastq.gz
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo FastQC completed in $RUNTIME s

# Upload FastQC output to NAS
STARTTIME=$(date +%s)
rsync -az $fastqc_output_folder/ coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "Results/02-scnRNA-seq_Processing/01-FastQC/", pool_id,"
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo FastQC upload completed in $RUNTIME s

# Run cellranger multi
STARTTIME=$(date +%s)
sed \"s|\\$TMPDIR|${TMPDIR}|g\" $SCRATCH/scripts/", pool_id, "_1-cellranger_config.csv > $TMPDIR/", pool_id, "_1-cellranger_config.csv
cellranger_output_folder=$TMPDIR/", pool_id,"_cellranger
mkdir $cellranger_output_folder
cd $cellranger_output_folder
", Sys.getenv("CELLRANGER_PATH"), " multi --id=", pool_id," \\
--csv=$TMPDIR/", pool_id, "_1-cellranger_config.csv \\
--jobmode=local \\
--localcores 16 \\
--localmem 115
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo Cell Ranger completed in $RUNTIME s

# Upload Cell Ranger output to NAS
STARTTIME=$(date +%s)
rsync -az $cellranger_output_folder/ coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "Results/02-scnRNA-seq_Processing/02-CellRanger
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo Cell Ranger upload completed in $RUNTIME s")
  writeLines(cellranger_script, paste0(output_path, pool_id, "_1-cellranger_script.sh"))
  
}
# Generates a Grid Engine bash script for CellBender ambient RNA removal.
# Writes <pool_id>_2-cellbender_script.sh to output_path.
# The job holds until the corresponding cellranger job completes.
# Args:
#   assay: one of "snC5pV2", "scC5pV2", or "scC5pV3OCM" — controls which
#          Cell Ranger output path and cell count extraction are used
#   sample_table: for scC5p assays, data frame with a sample_id column used to
#                 locate the per-sample metrics_summary.csv; pass "empty" for snC5pV2
#   pool_id: pool identifier used for job naming and output filenames
#   output_path: directory in which the generated .sh script will be written
cellbender <- function(email = Sys.getenv("EDDIE_EMAIL"),
                       assay = "snC5pV2", 
                       sample_table = "empty",
                       pool_id,
                       output_path) {
  if(assay == "snC5pV2" | pool_id == "CureAS_MG"){
    metrics_summary_path <- paste0("coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "Results/02-scnRNA-seq_Processing/02-CellRanger/", pool_id,"/outs/metrics_summary.csv")
    cellcount <- paste0("$(awk -F',' 'NR==2 {match($0, /\"[^\"]+\"/); val=substr($0, RSTART+1, RLENGTH-2); gsub(\",\", \"\", val); print val}' $TMPDIR/metrics_summary.csv)")
    h5_path <- paste0("coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "Results/02-scnRNA-seq_Processing/02-CellRanger/", pool_id,"/outs/raw_feature_bc_matrix.h5")
  } else if(assay == "scC5pV2" | assay == "scC5pV3OCM"){
    cellcount <- paste0('$(grep "Gene Expression" $TMPDIR/metrics_summary.csv | grep "GEX_1,Cells,\\"" | awk -F"\\"" \'{gsub(",", "", $2); print $2}\')')
    if(!is.character(sample_table)){
      metrics_summary_path <- paste0("coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "Results/02-scnRNA-seq_Processing/02-CellRanger/", pool_id,"/outs/per_sample_outs/", sample_table$sample_id[1],"/metrics_summary.csv")
    } else{
      metrics_summary_path <- paste0("coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "Results/02-scnRNA-seq_Processing/02-CellRanger/", pool_id,"/outs/per_sample_outs/", pool_id,"/metrics_summary.csv")
    }
    h5_path <- paste0("coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "Results/02-scnRNA-seq_Processing/02-CellRanger/", pool_id,"/outs/multi/count/raw_feature_bc_matrix.h5")
  }

  cellbender_script <- paste0("#!/bin/bash

#$ -N ", pool_id, "_2-cellbender
#$ -cwd
#$ -l h_rt=6:00:00
#$ -q gpu
#$ -l gpu-mig=1
#$ -e ./script_errors/", pool_id, "_2-cellbender_error.txt
#$ -o ./script_output/", pool_id,"_2-cellbender_output.txt
#$ -m beas
#$ -M ", email,"
#$ -hold_jid ", pool_id,"_1-cellranger

# Initialising the environment modules
. /etc/profile.d/modules.sh
source ~/.bashrc
module load cuda/12.1.1
mamba activate cellbender-env

# Download data from NAS
STARTTIME=$(date +%s)
rsync -az ", metrics_summary_path, " $TMPDIR/metrics_summary.csv
rsync -az ", h5_path, " $TMPDIR/feature_bc_matrix.h5
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo NAS download completed in $RUNTIME s

# Run cellbender
STARTTIME=$(date +%s)
cellbender_output_folder=$TMPDIR/", pool_id,"_cellbender
mkdir $cellbender_output_folder
cd $TMPDIR
cellranger_count=", cellcount, "
cellbender remove-background \\
                 --input $TMPDIR/feature_bc_matrix.h5 \\
                 --output $cellbender_output_folder/cellbender_output.h5 \\
                 --expected-cells $cellranger_count	 \\
                 --total-droplets-included 40000 \\
                 --cuda \\
                 --fpr 0.01 \\
                 --epochs 150
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo CellBender completed in $RUNTIME s

# Upload CellBender output to NAS
STARTTIME=$(date +%s)
rsync -az $cellbender_output_folder/ coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "Results/02-scnRNA-seq_Processing/03-CellBender/", pool_id,"
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo CellBender upload completed in $RUNTIME s")
  writeLines(cellbender_script, paste0(output_path, pool_id, "_2-cellbender_script.sh"))
}

# Generates a Grid Engine bash script for Souporcell SNP-based donor demultiplexing.
# Writes <pool_id>_3-souporcell_script.sh to output_path.
# The job holds until the corresponding cellbender job completes.
# Args:
#   assay: one of "snC5pV2", "scC5pV2", or "scC5pV3OCM" — controls BAM path
#   pool_id: pool identifier used for job naming and output filenames
#   n_pool: number of donors in the pool (passed as -k to souporcell_pipeline.py)
#   output_path: directory in which the generated .sh script will be written
souporcell <- function(email = Sys.getenv("EDDIE_EMAIL"),
                       assay = "snC5pV2",
                       pool_id,
                       n_pool,
                       output_path) {
  if(assay == "snC5pV2"){
    bam_loc <- paste0("coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "Results/02-scnRNA-seq_Processing/02-CellRanger/", pool_id,"/outs/possorted_genome_bam.bam")
    bam_bai_loc <- paste0("coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "Results/02-scnRNA-seq_Processing/02-CellRanger/", pool_id,"/outs/possorted_genome_bam.bam.bai")
  } else if(assay == "scC5pV2" | assay == "scC5pV3OCM"){
    bam_loc <- paste0("coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "Results/02-scnRNA-seq_Processing/02-CellRanger/", pool_id,"/outs/per_sample_outs/", pool_id,"/count/sample_alignments.bam")
    bam_bai_loc <- paste0("coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "Results/02-scnRNA-seq_Processing/02-CellRanger/", pool_id,"/outs/per_sample_outs/", pool_id,"/count/sample_alignments.bam.bai")
  }
  souporcell_script <- paste0("#!/bin/bash

#$ -N ", pool_id, "_3-souporcell
#$ -cwd
#$ -l h_rt=12:00:00
#$ -l h_rss=8G
#$ -pe sharedmem 16
#$ -e ./script_errors/", pool_id, "_3-souporcell_error.txt
#$ -o ./script_output/", pool_id,"_3-souporcell_output.txt
#$ -m beas
#$ -M ", email,"
#$ -hold_jid ", pool_id,"_2-cellbender

# Initialising the environment modules
. /etc/profile.d/modules.sh
source ~/.bashrc
module load singularity

# $TMPDIR is automatically set by gridengine, otherwise
# singularity will use /tmp which will fill
export SINGULARITY_TMPDIR=$TMPDIR/SINGULARITY_TMPDIR
mkdir $SINGULARITY_TMPDIR
cachedir=$SCRATCH/singularity
if [ ! -d \"$cachedir\" ]; then
  mkdir \"$cachedir\"
  echo \"Folder '$cachedir' created.\"
else
  echo \"Folder '$cachedir' already exists.\"
fi
export SINGULARITY_CACHEDIR=$cachedir

# Download data from NAS
STARTTIME=$(date +%s)
rsync -az coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "Results/02-scnRNA-seq_Processing/03-CellBender/", pool_id,"/cellbender_output_cell_barcodes.csv $TMPDIR/cellbender_output_cell_barcodes.tsv
rsync -az ", bam_loc, " $TMPDIR/sample_bam.bam
rsync -az ", bam_bai_loc, " $TMPDIR/sample_bam.bam.bai
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo NAS download completed in $RUNTIME s

# Run souporcell
STARTTIME=$(date +%s)
souporcell_output_folder=$TMPDIR/", pool_id,"_souporcell
mkdir $souporcell_output_folder
export SINGULARITY_BINDPATH=$GROUP_HOME,$SCRATCH,$TMPDIR
singularity exec $GROUP_HOME/software/sv_singularity/235ff94a409a24d2c758b92ed00a06db.sif souporcell_pipeline.py \\
    -i $TMPDIR/sample_bam.bam \\
    -b $TMPDIR/cellbender_output_cell_barcodes.tsv \\
    -f $GROUP_HOME/SV/references/refdata-gex-GRCh38-2024-A/fasta/genome.fa \\
    -t 16 \\
    -o $souporcell_output_folder \\
    -k ", n_pool, "
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo souporcell completed in $RUNTIME s

# Upload souporcell output to NAS
STARTTIME=$(date +%s)
rsync -az $souporcell_output_folder/ coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "Results/02-scnRNA-seq_Processing/04-souporcell/", pool_id,"
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo souporcell upload completed in $RUNTIME s")
  writeLines(souporcell_script, paste0(output_path, pool_id, "_3-souporcell_script.sh"))
}

# Generates a Grid Engine bash script for Velocyto RNA velocity quantification.
# Writes <pool_id>_4-velocyto_script.sh to output_path.
# Expects pre-staged filtered matrix files under $SCRATCH/velocyto_input/<pool_id>/.
# Args:
#   assay: one of "snC5pV2", "scC5pV2", or "scC5pV3OCM" — controls BAM path
#   pool_id: pool identifier used for job naming and output filenames
#   output_path: directory in which the generated .sh script will be written
velocyto <- function(email = Sys.getenv("EDDIE_EMAIL"),
                       assay = "snC5pV2",
                       pool_id,
                       output_path) {
  if(assay == "snC5pV2"){
    bam_loc <- paste0("coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "Results/02-scnRNA-seq_Processing/02-CellRanger/", pool_id,"/outs/possorted_genome_bam.bam")
    bam_bai_loc <- paste0("coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "Results/02-scnRNA-seq_Processing/02-CellRanger/", pool_id,"/outs/possorted_genome_bam.bam.bai")
  } else if(assay == "scC5pV2" | assay == "scC5pV3OCM"){
    bam_loc <- paste0("coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "Results/02-scnRNA-seq_Processing/02-CellRanger/", pool_id,"/outs/per_sample_outs/", pool_id,"/count/sample_alignments.bam")
    bam_bai_loc <- paste0("coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "Results/02-scnRNA-seq_Processing/02-CellRanger/", pool_id,"/outs/per_sample_outs/", pool_id,"/count/sample_alignments.bam.bai")
  }
  velocyto_script <- paste0("#!/bin/bash

#$ -N ", pool_id, "_4-velocyto
#$ -cwd
#$ -l h_rt=12:00:00
#$ -l h_rss=8G
#$ -pe sharedmem 16
#$ -e ./script_errors/", pool_id, "_4-velocyto_error.txt
#$ -o ./script_output/", pool_id,"_4-velocyto_output.txt
#$ -m beas
#$ -M ", email,"

# Initialising the environment modules
. /etc/profile.d/modules.sh
source ~/.bashrc
mamba activate velocyto-env

mkdir $TMPDIR/", pool_id,"
mkdir $TMPDIR/", pool_id,"/outs/

# Download data from NAS
STARTTIME=$(date +%s)
rsync -az ", bam_loc, " $TMPDIR/", pool_id,"/outs/possorted_genome_bam.bam
rsync -az ", bam_bai_loc, " $TMPDIR/", pool_id,"/outs/possorted_genome_bam.bam.bai
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo NAS download completed in $RUNTIME s

# Run velocyto
STARTTIME=$(date +%s)
mkdir $TMPDIR/", pool_id,"/outs/filtered_feature_bc_matrix
cp $SCRATCH/velocyto_input/", pool_id,"/barcodes.tsv.gz $TMPDIR/", pool_id,"/outs/filtered_feature_bc_matrix/
cp $SCRATCH/velocyto_input/", pool_id,"/features.tsv.gz $TMPDIR/", pool_id,"/outs/filtered_feature_bc_matrix/
cp $SCRATCH/velocyto_input/", pool_id,"/matrix.mtx.gz $TMPDIR/", pool_id,"/outs/filtered_feature_bc_matrix/
velocyto run10x $TMPDIR/", pool_id," ", Sys.getenv("VELOCYTO_REFERENCE"), "
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo velocyto completed in $RUNTIME s
                            
# Upload velocyto output to NAS
STARTTIME=$(date +%s)
rsync -az $TMPDIR/", pool_id,"/velocyto/ coconut:", Sys.getenv("COCONUT_CUREAS_PATH"), "Results/02-scnRNA-seq_Processing/05-velocyto/", pool_id,"
ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)
echo velocyto upload completed in $RUNTIME s")
  writeLines(velocyto_script, paste0(output_path, pool_id, "_4-velocyto_script.sh"))
}

# Scripts for data checking
checksum_featurecounts <- function(counts_df) {
  content <- counts_df[, c(1, 7:ncol(counts_df))]
  digest::digest(content, algo = "md5")
}

validate_featurecounts <- function(counts_df, sample_id, reference_checksums) {
  # Checksum only gene IDs and count columns (columns 7+), matching the bash approach
  content <- counts_df[, c(1, 7:ncol(counts_df))]
  observed <- digest::digest(content, algo = "md5")
  
  if (!sample_id %in% rownames(reference_checksums)) {
    warning(paste0("No reference checksum for ", sample_id, 
                   ". Storing: ", observed))
    return(observed)
  }
  
  expected <- reference_checksums[sample_id, "md5"]
  if (observed != expected) {
    stop(paste0("Checksum mismatch for ", sample_id,
                "\n  Expected: ", expected,
                "\n  Observed: ", observed))
  }
  message(paste0(sample_id, ": checksum OK"))
  return(observed)
}

# Dispatches to the right checksum method by file type
compute_checksum <- function(f) {
  if (grepl("\\.csv$", f)) {
    digest::digest(file = f, algo = "md5")
  } else if (grepl("\\.rds$", f)) {
    checksum_seurat(f)
  } else {
    digest::digest(file = f, algo = "md5")
  }
}

checksum_seurat <- function(path) {
  obj <- readRDS(path)
  
  if (inherits(obj, "Seurat")) {
    
    assay_name <- intersect(c("RNA", "Xenium"), names(obj@assays))
    
    if (length(assay_name) == 0) {
      stop("No RNA or Xenium assay found in object.")
    }
    
    if (length(assay_name) > 1) {
      stop("Both RNA and Xenium assays found. Please specify which one to use.")
    }
    
    assay <- obj@assays[[assay_name]]
    
    if (inherits(assay, "Assay5")) {
      counts <- SeuratObject::LayerData(
        obj,
        assay = assay_name,
        layer = "counts"
      )
    } else {
      counts <- assay@counts
    }
    
    content <- list(
      assay = assay_name,
      counts_digest = digest::digest(counts, algo = "md5"),
      metadata_digest = digest::digest(obj@meta.data, algo = "md5"),
      cell_names = Cells(obj),
      gene_names = rownames(counts)
    )
    
  } else {
    content <- obj
  }
  
  rm(obj)
  gc()
  
  digest::digest(content, algo = "md5")
}

generate_output_checksums <- function(output_path, checksum_file) {
  files <- list.files(output_path, 
                      pattern = "\\.(csv|rds)$", 
                      full.names = TRUE)
  checksums <- sapply(files, compute_checksum)
  df <- data.frame(
    file = basename(names(checksums)),
    md5 = unname(checksums)
  )
  write.csv(df, checksum_file, row.names = FALSE)
  cat(paste0("Saved ", nrow(df), " checksums to ", checksum_file, "\n"))
  return(df)
}

validate_output_checksums <- function(output_path, checksum_file) {
  reference <- read.csv(checksum_file)
  files <- list.files(output_path, 
                      pattern = "\\.(csv|rds)$", 
                      full.names = TRUE)
  
  results <- data.frame(file = character(), 
                        status = character(), 
                        stringsAsFactors = FALSE)
  
  for (f in files) {
    fname <- basename(f)
    observed <- compute_checksum(f)
    
    ref_row <- reference[reference$file == fname, ]
    
    if (nrow(ref_row) == 0) {
      status <- "NO REFERENCE"
      cat(paste0("WARNING  ", fname, ": no reference checksum\n"))
    } else if (observed == ref_row$md5) {
      status <- "OK"
      cat(paste0("OK       ", fname, "\n"))
    } else {
      status <- "MISMATCH"
      cat(paste0("MISMATCH ", fname, "\n",
                 "           Expected: ", ref_row$md5, "\n",
                 "           Observed: ", observed, "\n"))
    }
    results <- rbind(results, data.frame(file = fname, status = status))
  }
  
  missing <- setdiff(reference$file, basename(files))
  for (m in missing) {
    cat(paste0("MISSING  ", m, ": in reference but not found\n"))
    results <- rbind(results, data.frame(file = m, status = "MISSING"))
  }
  
  n_fail <- sum(results$status != "OK")
  if (n_fail > 0) {
    warning(paste0(n_fail, " file(s) failed validation"))
  } else {
    cat(paste0("\nAll ", nrow(results), " files validated OK\n"))
  }
  return(results)
}

generate_input_checksums <- function(paths, checksum_file_path, base_path = "../01_data/") {
  checksums <- sapply(paths, function(path) {digest::digest(file = path, algo = "md5")})
  relative <- sub(paste0("^", gsub("([.])", "\\\\\\1", base_path)), "", names(checksums))
  df <- data.frame(file = relative, md5 = unname(checksums))
  write.csv(df, checksum_file_path, row.names = FALSE)
  cat(paste0("Saved ", nrow(df), " input checksums to ", checksum_file_path, "\n"))
  return(df)
}

validate_input_checksums <- function(paths, checksum_file_path, base_path = "../01_data/") {
  reference <- read.csv(checksum_file_path)
  for (p in paths) {
    relative <- sub(paste0("^", gsub("([.])", "\\\\\\1", base_path)), "", p)
    observed <- digest::digest(file = p, algo = "md5")
    ref_row <- reference[reference$file == relative, ]
    if (nrow(ref_row) == 0) {
      cat(paste0("WARNING  ", relative, ": no reference checksum\n"))
    } else if (observed == ref_row$md5) {
      cat(paste0("OK       ", relative, "\n"))
    } else {
      stop(paste0("MISMATCH ", relative,
                  "\n  Expected: ", ref_row$md5,
                  "\n  Observed: ", observed))
    }
  }
}

# Draws a pheatmap with a subset of row labels, connected by line segments to
# their original positions. Useful when the heatmap has too many rows to label all.
# Adapted from https://stackoverflow.com/questions/52599180/partial-row-labels-heatmap-r
# Args:
#   pheatmap: a pheatmap object
#   kept.labels: character vector of row labels to display
#   repel.degree: value in [0, 1] controlling label spread (0 = minimal, 1 = full y-axis)
add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {
  library(grid)
  # repel.degree = number within [0, 1], which controls how much 
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis
  
  heatmap <- pheatmap$gtable
  
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  
  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, "")
  
  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant
    
    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }
      
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }
    
    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    
    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)
  
  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions
  
  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                     grobs = new.flag,
                                     t = 4, 
                                     l = 4
  )
  
  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  
  # plot result
  grid.newpage()
  grid.draw(heatmap)
  
  # return a copy of the heatmap invisibly
  invisible(heatmap)
}
# Loads a VCF file and returns a matrix of per-variant genotype calls.
# Row names are formatted as <chrom>:<pos>; column names are sample IDs from the VCF.
# Args:
#   vcf_file: path to a VCF file (hg38)
load_vcf_genotypes <- function(vcf_file) {
  # Read VCF file
  vcf <- readVcf(vcf_file, "hg38")  # Adjust genome version as needed
  
  # Extract genotype data
  genotypes <- geno(vcf)$GT
  rownames(genotypes) <- paste0(seqnames(rowRanges(vcf)), ":", start(rowRanges(vcf)))
  
  return(genotypes)
}

# Matches donors between two genotype matrices by computing pairwise genotype concordance.
# Returns a list with:
#   $SimilarityMatrix: matrix of concordance scores (rows = donors in genotypes1, cols = genotypes2)
#   $DonorMatches: data frame of best-match pairs and their concordance scores
# Args:
#   genotypes1, genotypes2: genotype matrices as returned by load_vcf_genotypes()
match_donors <- function(genotypes1, genotypes2) {
  # Find common variants between the two genotype sets
  common_variants <- intersect(rownames(genotypes1), rownames(genotypes2))
  if (length(common_variants) == 0) {
    stop("No common variants found between the two VCF files.")
  }
  
  # Subset the genotypes to the common variants
  genotypes1_common <- genotypes1[common_variants, ]
  genotypes2_common <- genotypes2[common_variants, ]
  
  # Compute similarity score between each donor pair
  similarity_matrix <- matrix(0, nrow = ncol(genotypes1_common), ncol = ncol(genotypes2_common))
  rownames(similarity_matrix) <- colnames(genotypes1_common)
  colnames(similarity_matrix) <- colnames(genotypes2_common)
  
  for (i in 1:ncol(genotypes1_common)) {
    for (j in 1:ncol(genotypes2_common)) {
      # Calculate the proportion of matching genotypes
      matching_genotypes <- sum(genotypes1_common[, i] == genotypes2_common[, j], na.rm = TRUE)
      total_variants <- length(common_variants)
      similarity_matrix[i, j] <- matching_genotypes / total_variants
    }
  }
  
  # Determine the best matches
  matching_results <- apply(similarity_matrix, 1, which.max)
  donor_matches <- data.frame(
    Donor1 = rownames(similarity_matrix),
    Donor2 = colnames(similarity_matrix)[matching_results],
    SimilarityScore = apply(similarity_matrix, 1, max)
  )
  
  return(list(SimilarityMatrix = similarity_matrix, DonorMatches = donor_matches))
}

# Runs the standard Harmony integration and clustering pipeline on a Seurat object:
#   (optionally) NormalizeData → variable feature selection → ScaleData → PCA →
#   RunHarmony (by donor_id) → FindNeighbors → FindClusters → RunUMAP.
# Returns the updated Seurat object.
# Args:
#   obj: Seurat object
#   dims: PCA/UMAP dimensions to use (default 1:30)
#   resolutions: vector of Louvain resolutions to test simultaneously
#   normalize: run NormalizeData before variable feature selection (default FALSE;
#              set TRUE when the object has been freshly subsetted/merged)
#   use_integration_features: use SelectIntegrationFeatures across per-donor splits
#                             instead of FindVariableFeatures (default FALSE; set TRUE
#                             for cell types with strong donor batch effects)
#   return.model: store the UMAP model for future projection (default TRUE)
#   vars.to.regress: character vector of metadata columns to regress out in ScaleData
#                    (e.g. "percent_ieg" to remove IEG artefacts; default NULL)
harmony_cluster <- function(obj,
                            dims = 1:30,
                            resolutions = c(0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4),
                            normalize = FALSE,
                            use_integration_features = FALSE,
                            return.model = TRUE,
                            vars.to.regress = NULL) {
  if (normalize) obj <- NormalizeData(obj)
  if (use_integration_features) {
    VariableFeatures(obj) <- SelectIntegrationFeatures(
      lapply(X = SplitObject(obj, split.by = "donor_id"),
             FUN = function(y) { FindVariableFeatures(y, verbose = FALSE) })
    )
  } else {
    obj <- FindVariableFeatures(obj)
  }
  obj <- ScaleData(obj, vars.to.regress = vars.to.regress) %>%
    RunPCA() %>%
    RunHarmony(group.by.vars = "donor_id") %>%
    FindNeighbors(reduction = "harmony", dims = dims) %>%
    FindClusters(resolution = resolutions) %>%
    RunUMAP(dims = dims, reduction = "harmony", return.model = return.model)
  return(obj)
}

# Plots a UMAP DimPlot for each clustering resolution stored in the Seurat object,
# to facilitate visual comparison across resolutions.
# Args:
#   obj: Seurat object with RNA_snn_res.<resolution> metadata columns
#   resolutions: vector of resolutions to plot (must match those used in FindClusters)
plot_cluster_umaps <- function(obj,
                               resolutions = c(0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3, 0.4)) {
  for (res in resolutions) {
    print(
      DimPlot(obj, reduction = "umap", group.by = paste0("RNA_snn_res.", res),
              cols = "polychrome", raster = FALSE) + coord_equal()
    )
  }
}

# Plots per-cluster QC violin plots and per-cell UMAP feature plots to aid
# identification of doublet/low-quality clusters before filtering.
# Args:
#   obj: Seurat object
#   cluster_col: metadata column to group violin plots by (e.g. "RNA_snn_res.0.3")
#   show_doublets: if TRUE, also plots scDblFinder.score and known_doublet overlay
#                  with a cluster × doublet contingency table (use TRUE for pre-filter
#                  QC review, FALSE for a lighter post-filter quality check)
plot_qc_by_cluster <- function(obj, cluster_col, show_doublets = TRUE) {
  vln_features <- c("nCount_RNA", "nFeature_RNA", "percent_mito")
  if (show_doublets) vln_features <- c(vln_features, "scDblFinder.score")
  for (f in vln_features) {
    print(VlnPlot(obj, features = f, group.by = cluster_col, pt.size = 0))
  }
  for (f in c("nCount_RNA", "nFeature_RNA", "percent_mito")) {
    print(FeaturePlot(obj, features = f) + coord_equal())
  }
  if (show_doublets) {
    print(FeaturePlot(obj, features = "scDblFinder.score", order = TRUE) + coord_equal())
    print(DimPlot(obj, group.by = "known_doublet") + coord_equal())
    print(table(obj$known_doublet, obj[[cluster_col, drop = TRUE]]))
  }
}

# Runs FindAllMarkers on a Seurat object, filters significant markers, and saves
# a dot plot and DEG table to output_path.
# Args:
#   x: Seurat object
#   name: label used in output filenames (<name>_deg_dot_plot.pdf, <name>_degs.csv)
#   output_path: directory in which output files will be written
FAM <- function(x, name, output_path){
  library(future)
  plan(multisession)
  DEG <- FindAllMarkers(x, only.pos = T, densify = T, max.cells.per.ident = 2000)
  DEG %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
  DEG$pct_diff <- abs(DEG$pct.2 - DEG$pct.1)
  DEG$sig <- "no"
  DEG$sig[DEG$pct_diff > 0.1 & DEG$p_val_adj < 0.05 & abs(DEG$avg_log2FC > 0.5)] <- "yes"
  plot <- DotPlot(x, features = unique(top10$gene), col.min = -1, col.max = 1) +
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
    scale_colour_gradientn(colours = scico(10, palette = "roma", direction = -1)) +
    guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+coord_equal()+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text=element_text(size=14, colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.title=element_blank(),
          panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
  ggsave(filename = paste0(output_path, name, "_deg_dot_plot.pdf"), plot = plot, width = 21, units = "in")
  write.csv(DEG, file = paste0(output_path, name, "_degs.csv"))
  plan(sequential)
}



# Figure functions
theme_minimal_fig <- function(
    axis_text_x = FALSE,
    axis_text_y = TRUE,
    legend = FALSE,
    margin = c(1, 1, 1, 1)
) {
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background  = element_blank(),
    axis.title = element_blank(),
    axis.ticks.x = if (axis_text_x) element_line(colour = "black", linewidth = 0.25) else element_blank(),
    axis.ticks.y = element_line(colour = "black", linewidth = 0.25),
    axis.text.x = if (axis_text_x) element_text(colour = "black", size = 5) else element_blank(),
    axis.text.y = if (axis_text_y) element_text(colour = "black", size = 5) else element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.25),
    plot.margin = unit(margin, "mm"),
    legend.position = if (legend) "right" else "none"
  )
}


