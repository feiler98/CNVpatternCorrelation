# imports
# ----------------------------------------------------------------------------------------------------------------------
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from pathlib import Path
# ----------------------------------------------------------------------------------------------------------------------


def vis_boxplot_cnv_spearman(dict_boxplots: dict,
                             title: str | None = None,
                             out_dir: Path | str | None = Path.cwd(),
                             show: bool = True,
                             for_boring_paper: bool = False):
    if out_dir is not None:
        out_dir = Path(out_dir).resolve()
        if not out_dir.exists() or not out_dir.is_dir():
            raise ValueError("Invalid 'out_dir' value provided!")

    if not for_boring_paper:  # only used by true adams
        plt.style.use('seaborn-v0_8-darkgrid')
        mpl.rcParams["boxplot.patchartist"] = True
        # rcParams for line & bar plots
        for param in ['figure.facecolor',
                      'axes.facecolor',
                      'savefig.facecolor']:
            plt.rcParams[param] = '#212946'  # bluish dark grey
        for param in ['text.color',
                      'axes.labelcolor']:
            plt.rcParams[param] = '0.9'
        for param in ['grid.color',
                      'xtick.color',
                      'ytick.color']:
            plt.rcParams[param] = "#8084d1"  # very light grey
        plt.rcParams['axes.prop_cycle'] = plt.cycler(color=["#fca138", "0.7"])  # orange

        # rcParams for boxplots
        for param in ['boxplot.whiskerprops.color',
                      'boxplot.flierprops.markeredgecolor',
                      'boxplot.capprops.color',
                      'boxplot.meanprops.markerfacecolor',
                      'boxplot.meanprops.markeredgecolor',
                      'boxplot.boxprops.color',
                      'boxplot.medianprops.color']:
            plt.rcParams[param] = "#fca138"  # orange

        for param in ['boxplot.flierprops.markerfacecolor',
                      'patch.facecolor']:
            plt.rcParams[param] = "#fca13870"
        plt.rcParams['boxplot.medianprops.linewidth'] = 2.5
        sns.despine()
    title = str(title) if title is not None else 'Spearman Correlation of Methylation x Pseudobulk-scRNA CNV Patterns'

    fig, ax = plt.subplots(figsize=(7+0.4*len(dict_boxplots), 7))
    ax.boxplot(dict_boxplots.values())
    ax.set_title(
        f"Boxplot | {title}")
    plt.xticks(rotation=45, ha="right")
    ax.set_xlabel("genomic sections")
    ax.set_ylabel("CNV spearman correlation [methyl x pseudobulk-scRNA]")
    ax.set_xticklabels(dict_boxplots.keys())
    if show:
        plt.show()
    if out_dir is not None:
        plt.savefig(out_dir / f"boxplot__{title.replace(" ", "_")}", dpi=300, bbox_inches="tight")
    plt.rcParams.update(plt.rcParamsDefault)