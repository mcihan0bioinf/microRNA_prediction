import pandas as pd
import argparse

# Function to predict miRNA expressions
def predict_mirnas(gene_expression_matrix, weights, intercepts, benchmarking):
    # Normalize column names to avoid mismatch
    gene_expression_matrix.columns = gene_expression_matrix.columns.str.strip().str.lower()
    gene_expression_matrix["gene_name"] = gene_expression_matrix["gene_name"].str.strip().str.lower()
    weights.columns = weights.columns.str.strip()
    weights["gene_name"] = weights["gene_name"].str.strip().str.lower()
    intercepts["mirna_name"] = intercepts["mirna_name"].str.strip().str.lower()

    # Ensure all required genes are present
    required_genes = set(weights["gene_name"])
    available_genes = set(gene_expression_matrix["gene_name"])
    missing_genes = required_genes - available_genes

    if missing_genes:
        print(f"Warning: The following required genes are missing: {', '.join(missing_genes)}")

    # Filter gene expression matrix to include only required genes
    filtered_gene_expression = gene_expression_matrix.loc[gene_expression_matrix["gene_name"].isin(required_genes)]
    filtered_gene_expression = filtered_gene_expression.set_index("gene_name")["expression_value"]

    # Ensure alignment between intercepts and filtered_gene_expression
    intercepts = intercepts.set_index("mirna_name")
    intercept_genes = intercepts.index.intersection(filtered_gene_expression.index)
    if len(intercept_genes) < len(filtered_gene_expression.index):
        missing_genes = set(filtered_gene_expression.index) - set(intercept_genes)
        print(f"Warning: The following genes are missing in the intercept file: {', '.join(missing_genes)}")

    # Align means and stds
    means = intercepts.loc[intercept_genes, "Intercept"]
    stds = intercepts.loc[intercept_genes, "miRNA_STD"]

    # Scale the gene expression values
    scaled_gene_expression = (filtered_gene_expression.loc[intercept_genes] - means) / stds

    # Align weights with the scaled genes
    weights = weights.set_index("gene_name")
    aligned_weights = weights.loc[scaled_gene_expression.index]

    # Predict miRNA expression using the weights
    raw_predictions = aligned_weights.T.dot(scaled_gene_expression)

    # Add intercepts to the predictions
    intercept_values = intercepts["Intercept"]
    final_predictions = raw_predictions + intercept_values

    # Convert predictions to a DataFrame
    predictions_df = final_predictions.reset_index()
    predictions_df.columns = ["miRNA", "Predicted_Expression"]

    # Merge predictions with R-squared values from benchmarking
    merged_results = predictions_df.merge(benchmarking[["microRNA", "R-squared"]], left_on="miRNA", right_on="microRNA", how="left")
    merged_results.drop(columns=["microRNA"], inplace=True)

    # Set R-squared values less than 0 to 0
    merged_results["R-squared"] = merged_results["R-squared"].apply(lambda x: max(0, x))

    # Sort by R-squared in descending order
    merged_results = merged_results.sort_values(by="R-squared", ascending=False)

    return merged_results

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Predict miRNA expressions.")
    parser.add_argument("-i", "--input", required=True, help="Path to the gene expression matrix file.")
    parser.add_argument("-o", "--output", required=True, help="Path to save the output file.")
    args = parser.parse_args()

    # Define file paths
    gene_expression_matrix_file = args.input
    output_file = args.output
    weights_file = "data/weight_coefficients/weight_coefficients_all.csv"
    intercepts_file = "data/weight_coefficients/intercepts.csv"
    benchmarking_file = "data/benchmarking_metrics/benchmark_results_all.csv"

    # Load the input files
    gene_expression_matrix = pd.read_csv(gene_expression_matrix_file)
    weights = pd.read_csv(weights_file)
    intercepts = pd.read_csv(intercepts_file)
    benchmarking = pd.read_csv(benchmarking_file)

    # Predict miRNA expressions and save the results
    results = predict_mirnas(gene_expression_matrix, weights, intercepts, benchmarking)
    results.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    main()