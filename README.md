# ECA_2020_Sed

The files ECA2_FTICR_BetaDisp.csv and VGC_texture.csv need to be in the repository for a piece of the code to run, but the contents are not used for any analysis associated with the manuscript.

The file merged_weights.csv contains moisture content data in the format used by the scripts, which deviates from the format used on the ESS-DIVE data package. It is read in by the scripts below.

For analyses in the manuscript, once FTICR data are processed, you can run the following scripts in the following order:
- 01_fticr_blanks_github.r "This file cleans the FTICR data so it is free of contamination. The outputs are stored in the repository, so it does not need to be rerun."
- 02_Generate_MCD_github.R "This file generates MCDs for each field site. The outputs are stored in the repository, so it does not need to be rerun."
- 03_FTICR_bNTI_create-nulls_github.R "This file generates null model runs. The outputs from this file are not stored in the repository due to a large number of files."
- 04_FTICR_bNTI_merge-nulls_github.R "This file compiles null model runs. For this script to run, it is necessary to first run 03_FTICR_bNTI_create-nulls_github.R. The outputs are stored in the repository, so it does not need to be rerun."
- 05_ECA_2020_ICR_Null_Figures_Stats_github.r "This file generates figures and does analyses. All figures in the paper should be generated by this script."
