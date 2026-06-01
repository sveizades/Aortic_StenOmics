from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
import os
os.makedirs("../03_output/11_S3_scRNAseq_CellPhoneDB/", exist_ok=True)
cpdb_results2 = cpdb_statistical_analysis_method.call(
    cpdb_file_path = '../01_data/v5.0.0/cellphonedb.zip',
    meta_file_path = '../03_output/10_S3_scRNAseq_annotations_combining/saltire3sc_atlas_celltypes.tsv',
    counts_file_path = '../03_output/10_S3_scRNAseq_annotations_combining/saltire3sc_atlas_counts_mtx.h5ad',
    counts_data = 'hgnc_symbol',             
    score_interactions = True,
    threshold = 0.1,
    threads = 15,
    debug_seed = 42,
    result_precision = 3,
    separator = '|',
    debug = False,
    output_path = '../03_output/11_S3_scRNAseq_CellPhoneDB/results',
    output_suffix = None,
    microenvs_file_path = "../01_data/microenvironments.tsv"
)