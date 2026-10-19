# imports
# ----------------------------------------------------------------------------------------------------------------------
import pandas as pd
# ----------------------------------------------------------------------------------------------------------------------


# Pseudobulk generation based on scRNA data
# here described for counts, not CNV patterns!
# --------------------------------------------------
# Crowell et al. (2020, Nature Communications) --> pseudobulk (mean / sum of counts)
# Squair et al. (2021, Nature Communications) --> aggregating counts across cells
# --------------------------------------------------
# it seems generally normal aggregation seems to be standard procedure
# Based on literature: use median (less sensitive to outliers) instead of mean, for the CNV patterns

def norm_bulk_sc_cnv(df_sc_cnv: pd.DataFrame) -> pd.DataFrame:
    """
    Parameters
    ----------
    df_sc_cnv: pd.DataFrame
        DataFrame must contain: ["CHR", "START", "END", cell_1,... cell_n] with n ∈ N (total cell pool)
    Returns
    -------
    pd.DataFrame
        DataFrame with shape: ["CHR", "START", "END", "bulk_profile"]
        with 'bulk_profile' calculated as the median of each genomic bin
    """
    df_idx_slice, df_cells = df_sc_cnv.copy()[["CHR", "START", "END"]], df_sc_cnv.copy().drop(["CHR", "START", "END"])
    # function across columns (cells), for every row (genomic bin)
    median_df = df_cells.median(axis=1).to_frame(name="bulk_profile")
    return pd.concat([df_idx_slice, median_df], axis=1)
