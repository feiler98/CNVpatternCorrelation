# imports
# ----------------------------------------------------------------------------------------------------------------------
import multiprocessing as mp
import os
import pandas as pd
import numpy as np
# ----------------------------------------------------------------------------------------------------------------------

# https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
dict_q_arm_start__hg38 = {
    'chr1': 123400000,
    'chr10': 39800000,
    'chr11': 53400000,
    'chr12': 35500000,
    'chr13': 17700000,
    'chr14': 17200000,
    'chr15': 19000000,
    'chr16': 36800000,
    'chr17': 25100000,
    'chr18': 18500000,
    'chr19': 26200000,
    'chr2': 93900000,
    'chr20': 28100000,
    'chr21': 12000000,
    'chr22': 15000000,
    'chr3': 90900000,
    'chr4': 50000000,
    'chr5': 48800000,
    'chr6': 59800000,
    'chr7': 60100000,
    'chr8': 45200000,
    'chr9': 43000000,
    'chrX': 61000000,
    'chrY': 10400000
}

def mult_algin_bins(df_rnaCNVbulk: pd.DataFrame,
                    df_methylCNV: pd.DataFrame) -> dict:
    """
    Parameters
    ----------
    df_rnaCNVbulk: pd.DataFrame
        DataFrame must contain: ["CHR", "START", "END", "bulk_profile"]
    df_methylCNV: pd.DataFrame
        DataFrame must contain: ["CHR", "START", "END", "bulk_profile"]

    Returns
    -------
    dict
        Dictionary with the curated / aligned genomic bins of Methylation-CNV and bulkified-scRNA-CNV.
    """

    list_chr_keys_order = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
                           "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                           "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
                           "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

    # create "index" column for both dataframes
    df_methylCNV = df_methylCNV.reset_index().drop("index", axis=1).reset_index()
    df_rnaCNV = df_rnaCNVbulk.reset_index().drop("index", axis=1).reset_index()
    df_methyl_idx, df_methyl_bulk = (df_methylCNV[["index", "CHR", "START", "END"]],
                                df_methylCNV[["bulk_profile"]])
    df_rna_idx, df_rna_bulk = (df_rnaCNV[["index", "CHR", "START", "END"]],
                                df_rnaCNV[["bulk_profile"]])

    # available chromosomes
    chr_available_list = [chr_tag for chr_tag in list(df_methylCNV["CHR"].unique()) if
                          chr_tag in list(df_rnaCNV["CHR"].unique())]
    # multiprocessing
    with mp.Pool(os.cpu_count() if len(chr_available_list) >= os.cpu_count() else len(chr_available_list)) as pool:
        df_slice_list = pool.starmap(mult_chr_bin_match, [
            (chr_key,
             df_methyl_idx.where(df_methyl_idx["CHR"] == chr_key).dropna(subset="CHR"),
             df_rna_idx.where(df_rna_idx["CHR"] == chr_key).dropna(subset="CHR"))
             for chr_key in chr_available_list])
    # concat list of dictionaries
    dict_chr_df_slice = {k: v for dict_slice in df_slice_list for k, v in dict_slice.items()}
    # concat and sort by list_chr_keys_order
    df_concat = pd.concat(
        [dict_chr_df_slice[chr_key] for chr_key in list_chr_keys_order if chr_key in chr_available_list], axis=0)
    df_concat = df_concat.astype({"START": int,
                                  "END": int,
                                  "indexMethyl": int,
                                  "indexRNA": int})

    df_rnaCNVbulk_reshape = df_rnaCNV.copy().loc[list(df_concat["indexRNA"]), :]
    df_methylCNV_reshape = df_methylCNV.copy().loc[list(df_concat["indexMethyl"]), :]

    return {"rna_cnv_reshaped": df_rnaCNVbulk_reshape,
            "methyl_cnv_reshaped": df_methylCNV_reshape}


def mult_chr_bin_match(chr_tag: str,
                       df_CNV_slice1: pd.DataFrame,
                       df_CNV_slice2: pd.DataFrame) -> dict:
    """
    Parameters
    ----------
    chr_tag: str
    df_CNV_slice1: pd.DataFrame
    df_CNV_slice2: pd.DataFrame

    Returns
    -------
    dict
        Contains the chr_tag key and the DataFrame slice
    """

    df_CNV_slice1 = df_CNV_slice1.sort_values("START")
    df_CNV_slice2 = df_CNV_slice2.sort_values("START")
    dict_bin_s1CNV = {col: np.array(df_CNV_slice1[col]) for col in df_CNV_slice1.columns}
    dict_bin_s2CNV = {col: np.array(df_CNV_slice2[col]) for col in df_CNV_slice2.columns}
    chr_list = []
    start_list = []
    end_list = []
    index1 = []
    index2 = []
    len_dict_END = len(dict_bin_s2CNV["END"])
    start_pos_rna = dict_bin_s2CNV["START"][0]

    for idx, chr_tag, start_pos, end_pos in zip(dict_bin_s1CNV["index"], dict_bin_s1CNV["CHR"],
                                                dict_bin_s1CNV["START"], dict_bin_s1CNV["END"]):
        start_idx = 0 if start_pos < start_pos_rna else np.where(dict_bin_s2CNV["START"] <= start_pos)[0][-1]
        end_idx = np.argmax(dict_bin_s2CNV["END"] > end_pos) + 1
        end_idx = end_idx if end_idx <= len_dict_END else len_dict_END

        # small AI silliness since end pos is 0 at the end due to overflow
        end_idx = len(dict_bin_s2CNV["END"]) if start_idx == len(dict_bin_s2CNV["START"]) - 1 else end_idx

        if start_idx >= end_idx:
            continue
        else:
            start_slice = dict_bin_s2CNV["START"].copy()[start_idx:end_idx]
            end_slice = dict_bin_s2CNV["END"].copy()[start_idx:end_idx]

            start_slice_idx = start_slice[0]
            end_slice_idx = end_slice[-1]
            start_replace = start_pos if start_slice_idx < start_pos else start_slice_idx
            end_replace = end_pos if end_slice_idx > end_pos else end_slice_idx

            if start_replace >= end_slice[0]:
                start_idx += 1
                start_slice = np.delete(start_slice, 0)
                end_slice = np.delete(end_slice, 0)
            else:
                start_slice[0] = start_replace
            if len(start_slice) == 0:
                continue
            if end_replace <= start_slice[-1]:
                end_idx -= 1
                start_slice = np.delete(start_slice, -1)
                end_slice = np.delete(end_slice, -1)
            else:
                end_slice[-1] = end_replace
            if len(start_slice) == 0:
                continue

            # extend lists
            # ------------
            chr_list.extend(dict_bin_s2CNV["CHR"][start_idx:end_idx])
            start_list.extend(start_slice)
            end_list.extend(end_slice)
            index1.extend([idx] * (end_idx - start_idx))
            index2.extend(dict_bin_s2CNV["index"][start_idx:end_idx])

    return {chr_tag: pd.DataFrame({"CHR": chr_list,
                                   "START": start_list,
                                   "END": end_list,
                                   "indexMethyl": index1,
                                   "indexRNA": index2})}


def bin_size_norm(df_reshape: pd.DataFrame, bin_size: int = 10_000) -> pd.DataFrame:
    """
    Parameters
    ----------
    df_reshape: pd.DataFrame
        DataFrame must contain: ["CHR", "START", "END", "bulk_profile"].
        Should be applied post alignment of dataframe genomic bins!
    bin_size: int
        Size of new bins, values rounded

    Returns
    -------
    pd.DataFrame
    """

    df_reshape.reset_index(inplace=True)
    df_reshape.index.name = "std_bin_idx"
    list_difference = list(df_reshape["END"] - df_reshape["START"])
    list_rebin_idx = []
    list_rebin_idx.extend(idx for idx, diff in zip(list(df_reshape.index), list_difference) for _ in range(round(diff/bin_size)))
    return df_reshape.copy().loc[list_rebin_idx, :]


