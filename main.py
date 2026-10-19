# imports
# ----------------------------------------------------------------------------------------------------------------------
import pandas as pd
from pathlib import Path
# local imports
from utils import (save_as_json_dict,
                   bulk_sc_cnv,
                   mult_algin_bins,
                   bin_size_norm)
# ----------------------------------------------------------------------------------------------------------------------


##############
# INITIALIZE #
##############

bin_size_set = 10_000
p_file_methyl_scRNA_connect = Path("/run/media/feilerwe/WF_SUPER/hannah_data_analysis/PFSE_Methylierungs_match/KK__PFSE__scRNA__Methyl__link.csv")
p_dir_scRNA_csv = Path("/run/media/feilerwe/WF_SUPER/hannah_data_analysis/PFSE_Methylierungs_match/PF-SE/KK_mtx_out")
out_path = Path("/run/media/feilerwe/WF_SUPER/hannah_data_analysis/PFSE_Methylierungs_match/results_corr")

if __name__ == "__main__":
    if not out_path.exists():
        out_path.mkdir(exist_ok=True, parents=True)

    # import central df
    df_conn__meth_scRNA = pd.read_csv(p_file_methyl_scRNA_connect, index_col = "ID").dropna(subset=["Sentrix ID path", "KK-Nummer - fixed RNA"])
    dict_json_out = {}  # add all different comparison options with name for visualization
    for idx, row in df_conn__meth_scRNA.iterrows():
        dict_rename_cols = {
            'Chromosome': 'CHR',
            'Start': 'START',
            'End': 'END',
            Path(row["Sentrix ID path"]).stem.replace("_CNV", ""): "bulk_profile"
        }
        # dataframes import
        methyl_df = pd.read_csv(row["Sentrix ID path"], sep="\t").drop("Feature", axis=1).rename(dict_rename_cols, axis=1)
        scRNA_df = pd.read_csv(p_dir_scRNA_csv / f"{row['KK-Nummer - fixed RNA']}__cnv_mtx.csv", index_col="Gene")
        scRNA_df.index.name = "index"
        scRNA_df_bulk = bulk_sc_cnv(scRNA_df.copy())
        # alignment
        dict_chr_bin_algined = mult_algin_bins(scRNA_df_bulk, methyl_df)
        print(dict_chr_bin_algined)
        df_rebin__rna = bin_size_norm(dict_chr_bin_algined["rna_cnv_reshaped"], bin_size_set)
        df_rebin__methyl = bin_size_norm(dict_chr_bin_algined["methyl_cnv_reshaped"], bin_size_set)


