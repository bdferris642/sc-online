library(fitdistrplus)
library(harmony)
library(Seurat)
library(dplyr)
library(qs)

N_BOOTS = 5

MEAN_PROPS = c(
    "SOX6_AGTR1" = 0.1,
    "SOX6_NXPH1" = 0.1,
    "SOX6_GFRA2" = 0.1,
    "SOX6_GADPH" = 0.1,
    "MT" = 0.1,
    "CALB1_CALCR_NTNG1" = 0.1,
    "CALB1_CRYM" = 0.1,
    "CALB1_TRHR" = 0.1,  
    "CALB1_PPPR1R7" = 0.1,
    "CALB1_GAD2_CHRM2" = 0.1,
    "CALB1_GEM" = 0.1,
    "CALB1_RBP4_SST" = 0.1
)

VAR_PROPS = c(
    "SOX6_AGTR1" = 0.1,
    "SOX6_NXPH1" = 0.1,
    "SOX6_GFRA2" = 0.1,
    "SOX6_GADPH" = 0.1,
    "MT" = 0.1,
    "CALB1_CALCR_NTNG1" = 0.1,
    "CALB1_CRYM" = 0.1,
    "CALB1_TRHR" = 0.1,  
    "CALB1_PPPR1R7" = 0.1,
    "CALB1_GAD2_CHRM2" = 0.1,
    "CALB1_GEM" = 0.1,
    "CALB1_RBP4_SST" = 0.1
)

BASE_PATH = "/mnt/accessory/seq_data/pd_all/240514"

# R
# Define the list of vectors
# alist = list(
#     c(0.5, 0.3, 0.2),
#     c(0.4, 0.4, 0.2),
#     c(0, 0.4, 0.6)
#     # ...
# )

# # Initialize empty vectors to store the shape parameters
# shape1 <- numeric(length(alist))
# shape2 <- numeric(length(alist))

# # Loop over each vector in the list
# for (i in seq_along(alist)) {
#     # Calculate the mean and variance of the vector
#     mean_i <- mean(alist[[i]])
#     var_i <- var(alist[[i]])
    
#     # Estimate the shape parameters for a Beta distribution
#     shape1[i] <- ((1 - mean_i) / var_i - 1 / mean_i) * mean_i ^ 2
#     shape2[i] <- shape1[i] * (1 / mean_i - 1)
# }

# print(shape1)
# print(shape2)


# Load the data
# (1) Get cluster centers
# Bootstrap:
# (2) Split controls by participant
    # In the mock-cases, subsample the DANs 
# (3) Integrate everything with harmony and cluster at a variety of resolutions
    # Pick the lowest resolution that gives you the “right” number of clusters
    # If you skip, default to the lower end and fewer clusters
# (4) Label new clusters by their proximity to the original cluster centers
# (5) Make pseudocells
    # TODO: Run Peer?
# (6) Run DE, make a distribution of logFCs

nurr_da = qread(file.path(BASE_PATH, "nurr_seurat_da_clean_clustered.qs"))
nurr_da_ctr = qread(file.path(BASE_PATH, "nurr_seurat_da_ctr_clean_clustered.qs"))

nurr_da_ctr_counts = as.data.frame(nurr_da_ctr@assays$RNA@counts)
nurr_da_ctr_counts$da_subtype = nurr_da_ctr$da_subtype
nurr_da_ctr_counts$da_subtype

# is this valid? Should I be grouping by individuals first?
nurr_da_ctr_counts_cluster_centers = (
    nurr_da_ctr_counts %>%
    group_by(da_subtype) %>%
    summarize_all(mean)
)

# in our own data, 41% of participants are cases

for (i in 1:N_BOOTS){
    
    # randomly assign 41% of the participants to the mock_case group
    nurr_da_ctr_participants = unique(nurr_da_ctr$participant_id)
    mock_case_control = sample(
        c("pd", "ctr"), 
        size = length(nurr_da_ctr_participants), 
        replace = TRUE, 
        prob = c(0.41, 0.59))
    names(mock_case_control) = nurr_da_ctr_participants

    nurr_da_ctr$mock_case_control = mock_case_control[nurr_da_ctr$participant_id]

    # subsample the mock cases
    nurr_da_ctr__mock_ctr = nurr_da_ctr[,nurr_da_ctr$mock_case_control == "ctr"]
    nurr_da_ctr__mock_pd = nurr_da_ctr[,nurr_da_ctr$mock_case_control == "pd"]

    for participant in unique(nurr_da_ctr__mock_pd$participant_id){
        nurr_da_ctr__mock_pd_participant = nurr_da_ctr__mock_pd[,nurr_da_ctr__mock_pd$participant_id == participant]
        
        for (da_subtype in unique(nurr_da_ctr__mock_pd_participant$da_subtype)){
            nurr_da_ctr__mock_pd_participant = subsample_category(
                obj=nurr_da_ctr__mock_pd_participant, 
                sample_col = "da_subtype", 
                sample_val = da_subtype, 
                #sample_prop = 
            )
        }


        nurr_da_ctr__mock_pd_participant_subsampled = nurr_da_ctr__mock_pd_participant[, sample(colnames(nurr_da_ctr__mock_pd_participant), size = 2999, replace = FALSE)]
    }


    pbmc.subsampled <- pbmc[, sample(colnames(pbmc), size =2999, replace=F)]

    
}



