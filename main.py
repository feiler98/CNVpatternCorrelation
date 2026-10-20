# imports
# ----------------------------------------------------------------------------------------------------------------------
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr
# local imports
from utils import (save_as_json_dict,
                   dict_q_arm_start__hg38,
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
    # add all different comparison options with name for visualization
    genome_dict = {"total_genome": []}
    chr_dict = {}
    chr_arm_dict = {}

    for idx, row in df_conn__meth_scRNA.iterrows():
        dict_json_out = {"sample_name": row['KK-Nummer - fixed RNA'], "chr_seg": {}, "chr_seg_arm": {}}
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
        df_bin_align__rna = dict_chr_bin_algined["rna_cnv_reshaped"].drop("index", axis=1).reset_index()
        df_bin_align__methyl = dict_chr_bin_algined["methyl_cnv_reshaped"].drop("index", axis=1).reset_index()
        # concat and reassign START and END positions based on rule: take bigger number START, take lower number STOP per index
        df_concat_signatures = pd.concat([df_bin_align__rna[["CHR", "bulk_profile"]].rename({"bulk_profile": "bulk_profile_rna"}, axis=1),
                                          df_bin_align__methyl[["bulk_profile"]].rename({"bulk_profile": "bulk_profile_methyl"}, axis=1)],
                                         axis=1)
        df_concat_signatures["START"] = [pos_rna if pos_rna >= pos_methyl else pos_methyl for pos_rna, pos_methyl in
                                         zip(df_bin_align__rna["START"], df_bin_align__methyl["START"])]
        df_concat_signatures["END"] = [pos_rna if pos_rna <= pos_methyl else pos_methyl for pos_rna, pos_methyl in
                                         zip(df_bin_align__rna["END"], df_bin_align__methyl["END"])]
        df_normbin_signatures = bin_size_norm(df_concat_signatures, bin_size_set)

        # Spearman correlation
        # --------------------------------------------------------------------------------------------------------------
        # full genome
        spearman_score = float("%.3f" % round(spearmanr(df_normbin_signatures["bulk_profile_rna"].to_numpy(), df_normbin_signatures["bulk_profile_methyl"].to_numpy(), alternative="two-sided").statistic, 3))
        dict_json_out["total_genome"] = spearman_score
        list_append = genome_dict["total_genome"].copy()
        list_append.append(spearman_score)
        genome_dict["total_genome"] = list_append

        # chr segmented
        for chr_id in list(df_concat_signatures["CHR"].unique()):
            chr_df_slice = df_concat_signatures.where(df_concat_signatures["CHR"]==chr_id).dropna(subset="CHR")
            spearman_score = float("%.3f" % round(spearmanr(chr_df_slice["bulk_profile_rna"], chr_df_slice["bulk_profile_methyl"], alternative="two-sided").statistic, 3))
            dict_extend = dict_json_out["chr_seg"].copy()
            dict_extend[chr_id] = spearman_score
            dict_json_out["chr_seg"] = dict_extend
            if chr_id not in chr_dict.keys():
                chr_dict[chr_id] = [spearman_score]
            else:
                list_append = chr_dict[chr_id].copy()
                list_append.append(spearman_score)
                chr_dict[chr_id] = list_append
            # chr arm segmented --> ref: dict_q_arm_start__hg38
            chr_df_slice_p = df_concat_signatures.where((df_concat_signatures["CHR"] == chr_id) &
                                                        (df_concat_signatures["START"] <= dict_q_arm_start__hg38[chr_id])).dropna(subset="CHR")
            chr_df_slice_q = df_concat_signatures.where((df_concat_signatures["CHR"] == chr_id) &
                                                        (df_concat_signatures["START"] > dict_q_arm_start__hg38[chr_id])).dropna(subset="CHR")
            for tag_arm, df_slice_chr_arm in zip(["p", "q"], [chr_df_slice_p, chr_df_slice_q]):
                full_tag = f"{chr_id}_{tag_arm}"
                spearman_score = float("%.3f" % round(
                    spearmanr(df_slice_chr_arm["bulk_profile_rna"], df_slice_chr_arm["bulk_profile_methyl"],
                              alternative="two-sided").statistic, 3))
                dict_extend = dict_json_out["chr_seg_arm"].copy()
                dict_extend[full_tag] = spearman_score
                dict_json_out["chr_seg_arm"] = dict_extend
                if full_tag not in chr_dict.keys():
                    chr_arm_dict[full_tag] = [spearman_score]
                else:
                    list_append = chr_arm_dict[full_tag].copy()
                    list_append.append(spearman_score)
                    chr_arm_dict[full_tag] = list_append
        # --------------------------------------------------------------------------------------------------------------
        print(dict_json_out)
        print(genome_dict)
        print(chr_dict)
        print(chr_arm_dict)

