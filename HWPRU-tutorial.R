##############################################################################
# R-Tutorial: Sampling algorithm
##############################################################################

# We were asked to provide a statistical approach through which to "match"
# schools based upon defined variables 
# (obesity prevalence, deprivation and ethnicity).

##### Step 1: If necessary install relevant libraries

# you can also install via the packages tab (right panel) 
install.packages("dplyr", "tidyr", "readr")

# Step 2: Load libraries
library(dplyr) # supports data manipulation
library(tidyr) # allows you to organise and reshape data
library(readr) # faster read in of CSV documents

# Step 3: Import the csv. I have set this up as an accessible GITHUB
# repository. You would usually have to "direct" R to the relevant folder
# and file.
# 
schls_data <- read_csv("https://raw.githubusercontent.com/UCLJRAHILLY/HWPRU-Tutorial060725/refs/heads/main/HWPRU_tutorial2.csv")

# Step 4: Real‑world CSVs often embed non‑numeric clutter in numeric column
# (e.g “5,000”, “23%). We therefore need to clean the data.
# The code below strips everything except digits 
# and a decimal point.

# We can create a simple function that will convert the data for us
# and be reusable/reproducible

clean_numeric <- function(x) as.numeric(gsub("[^0-9.]", "", x))

# Then we can apply it to calculate derived fields, 
# and drop incomplete rows:

schls_data_amend <- schls_data %>%
 mutate(
  sampleSize = clean_numeric(Year.6.count.NCMP), # converts the school headcount
  obesityRate = clean_numeric(Very.Overweight.Prevalance....) #removes % from the obesity prevalence value
  ) %>% filter(
    if_all(c("ethnicity1", "ethnicity2", 
             "Pupil.Index.Of.Multiple.Deprivation.Score", 
             "Very.Overweight.Prevalance...."), ~ !is.na(.)) # Keeps only rows with values
  )

# Which school do you think might be the best match for school_1?
schls_data_amend[schls_data_amend$School.Name == "school_1",
                 c("sampleSize","obesityRate",
                   "ethnicity1","ethnicity2",
                   "Pupil.Index.Of.Multiple.Deprivation.Score")]

# Look at the raw dataframe

View(schls_data_amend)
###########################################################################
# Question: How many records had NAs?
###########################################################################

# As we are using variables that are on different scales
# we need to normalise - using z-score
# This transforms all variables to a common scale with a mean of 0 and 
# a standard deviation of 1, ensuring that each variable contributes 
# equally to the matching. 
# This normalisation makes the comparison fair and 
# improves the accuracy of the matching process.

schls_data_amend <- schls_data_amend %>%
  mutate(
    obesityPrev_z = scale(obesityRate)[,1], #scale is an inbuilt function that can do this for you
    IMD_z = scale(Pupil.Index.Of.Multiple.Deprivation.Score)[,1],
    ethnicity1_z = scale(ethnicity1)[,1],
    ethnicity2_z = scale(ethnicity2)[,1]
  )

# separate out the variables we want to match on by name
# By explicitly listing these variable names, we can easily select and work
# with only these specific variables later in the analysis.

match_vars <- c("ethnicity1_z",
                "ethnicity2_z",
                "IMD_z",
                "obesityPrev_z")

# we can utilise existing functions to calculate the Euclidean distance
# between our variables.
# Euclidean distance quantifies how different two observations 
# (e.g., schools) are across the selected variables we're matching on — 
# such as ethnicity proportions, deprivation (IMD), and obesity prevalence.

# Since all variables have been standardised (using z-scores), the Euclidean distance is essentially the overall difference between two schools, taking into account all matching variables at once.

# A smaller distance means the schools are more similar across those variables.

# A larger distance means they differ more significantly.

# So in this case, Euclidean distance is a single number that summarizes the combined difference across all selected variables — allowing us to identify which schools are most alike (and therefore good matches).
dist_matrix <- dist(schls_data_amend[, match_vars], method = "euclidean")
dist_mat    <- as.matrix(dist_matrix)     

# Adds names to our matrix so that we can inspect it afterwards
rownames(dist_mat) <- schls_data_amend$School.Name
colnames(dist_mat) <- schls_data_amend$School.Name

###########################################################################
 # Inspect the dist_mat dataframe
###########################################################################

# Get distances from school_1 to all others
dist_to_school1 <- dist_mat[which(schls_data_amend$School.Name == "school_1"), ]

# Exclude self-distance (zero)
dist_to_school1[which(schls_data_amend$School.Name == "school_1")] <- NA

# Find index of closest school
best_match_idx <- which.min(dist_to_school1)

# Get school name of best match
if(length(best_match_idx) > 0) {
  best_match_school <- names(dist_to_school1)[best_match_idx]
} else {
  best_match_school <- NA
}
cat("On Euclidean distance the 
    best matching school for school_1 is:", best_match_school, "\n")

###########################################################################
# While Euclidean distance treats all variables equally and assumes 
# they are uncorrelated, it can be misleading if the variables you're 
# matching on are correlated or have different variances
# 
# In our case, some of the matching variables (like deprivation and obesity
# prevalence, or ethnicity proportions) may be correlated with each other.
# This means that differences along those dimensions might be redundant or
# overemphasized in a simple Euclidean distance calculation
# 
# Mahalanobis distance solves this by:

# Taking into account the correlation structure of the variables.

# Adjusting the distance calculation to avoid overweighting related or redundant variables.

# Providing a more accurate measure of true multivariate difference between observations
###########################################################################

# Define the variables to use for matching
# We're keeping only one ethnicity variable to avoid redundancy
match_vars <- c("ethnicity1_z",   # standardised ethnicity variable
                #"ethnicity2_z",  # this one is excluded as it is redundant
                "IMD_z",          # standardised deprivation score
                "obesityPrev_z")  # standardised obesity prevalence

# Subset the main data to include only the selected matching variables
# This creates a matrix-like data frame used in distance calculations
match_data <- schls_data_amend[, match_vars]

# Calculate the covariance matrix of the matching variables
# This is needed for Mahalanobis distance, which adjusts for correlations
cov_mat <- cov(match_data)

# Find the row number (index) for the school we want to match from
idx_school1 <- which(schls_data_amend$School.Name == "school_1")

# Extract the matching variable values for "school_1"
# This becomes the reference point for calculating distances
school1_vec <- match_data[idx_school1, ]

# Calculate Mahalanobis distance from "school_1" to every other school
# This accounts for correlation and scale among the variables
mahal_dist <- apply(match_data, 1, function(x) 
  mahalanobis(as.numeric(x), center = as.numeric(school1_vec), cov = cov_mat)
)

# Remove the self-distance (distance from "school_1" to itself), which is 0
# We set it to NA so it won't be picked as the closest match
mahal_dist[idx_school1] <- NA

# Identify the index of the school with the smallest Mahalanobis distance
# This is the best match for "school_1" based on the selected variables
best_match_idx_mahal <- which.min(mahal_dist)

#Get school name of best Mahalanobis match
best_match_school_mahal <- schls_data_amend$School.Name[best_match_idx_mahal]

cat("On Mahalanobis distance the 
    best matching school for school_1 is:", 
    best_match_school_mahal, "\n")

############################################################################



