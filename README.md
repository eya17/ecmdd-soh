# Evidential c-Medoids with DTW Distance

Time series clustering implementation using Evidential c-Medoids algorithm with Dynamic Time Warping distance metric.

## Requirements
```R
install.packages(c("evclust", "cluster", "dtw"))
```

## Usage
Run the script directly using:
```bash
Rscript ecmdd_soh.R
```

This will:
1. Load time series data from SOH.csv
2. Compute DTW distance matrix
3. Perform clustering with 3 clusters
4. Output clustering results and silhouette analysis

## Author
Eya Laffet
