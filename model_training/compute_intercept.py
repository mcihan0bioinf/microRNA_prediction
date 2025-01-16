import pandas as pd
import argparse

# Parse arguments
parser = argparse.ArgumentParser(description="Compute intercepts for miRNA data.")
parser.add_argument("--feature-matrix", required=True, help="Path to the feature matrix CSV file.")
parser.add_argument("--label-matrix", required=True, help="Path to the label matrix CSV file.")
parser.add_argument("--coefficients", required=True, help="Path to the coefficients CSV file.")
parser.add_argument("--output-intercepts", required=True, help="Path to save the output intercepts CSV file.")
args = parser.parse_args()

# Paths to input files
feature_matrix_file = args.feature_matrix  # Gene expression matrix (raw features)
label_matrix_file = args.label_matrix      # miRNA expression matrix (raw targets)
coefficients_file = args.coefficients      # Coefficients for each miRNA
output_intercepts_file = args.output_intercepts  # File to store computed intercepts

# Load data
feature_matrix = pd.read_csv(feature_matrix_file)  # Features: Genes in columns, samples as rows
label_matrix = pd.read_csv(label_matrix_file)      # Labels: miRNAs in columns, samples as rows
coefficients = pd.read_csv(coefficients_file, index_col=0)  # Coefficients: Genes in rows, miRNAs in columns

# Transpose feature and label matrices
feature_matrix = feature_matrix.T  # Genes become rows, samples become columns
label_matrix = label_matrix.T      # miRNAs become rows, samples become columns

# Add headers to feature and label matrices if necessary
feature_matrix.columns = [f"Sample_{i+1}" for i in range(feature_matrix.shape[1])]
label_matrix.columns = [f"Sample_{i+1}" for i in range(label_matrix.shape[1])]

# Ensure alignment of indices
coefficients.index = coefficients.index.str.strip().str.lower()  # Gene names
feature_matrix.index = feature_matrix.index.str.strip().str.lower()
label_matrix.index = label_matrix.index.str.strip().str.lower()

# Compute raw means and standard deviations for features (genes)
raw_feature_means = feature_matrix.mean(axis=1)  # Mean for each gene
raw_feature_stds = feature_matrix.std(axis=1)    # Standard deviation for each gene

# Compute raw means for targets (miRNAs)
target_means = label_matrix.mean(axis=1)         # Mean for each miRNA
target_stds = label_matrix.std(axis=1)           # Standard deviation for each miRNA (new)

# Adjust coefficients to raw scale
adjusted_weights = coefficients.div(raw_feature_stds, axis=0)

# Compute weighted sum of raw feature means
weighted_sum = adjusted_weights.T.dot(raw_feature_means)

# Compute intercepts
intercepts = target_means - weighted_sum

# Create the intercept DataFrame
intercepts_df = pd.DataFrame({
    "mirna_name": intercepts.index,
    "Intercept": intercepts.values,
    "miRNA_STD": target_stds.values,  # Add the standard deviation of miRNA
    "z_scale_factor": (target_means.values / target_stds.values)  # Add the z-scaling factor
})

# Save intercepts to file
intercepts_df.to_csv(output_intercepts_file, index=False)
print(f"Intercepts with miRNA STD and z-scaling factors saved to {output_intercepts_file}")
