import pandas as pd
import numpy as np
import argparse

def filter_matrices(feature_matrix_path, label_matrix_path, threshold_feature, threshold_label, output_dir):
    # Load raw matrices
    feature_matrix = pd.read_csv(feature_matrix_path)
    label_matrix = pd.read_csv(label_matrix_path)

    # Preserve the first row with sample names
    feature_matrix_header = feature_matrix.iloc[0, :]
    label_matrix_header = label_matrix.iloc[0, :]

    # Filter rows with too many zeros or NaNs in feature matrix
    rows_to_drop_feature = feature_matrix.iloc[1:].index[
        ((feature_matrix.iloc[1:] == 0) | feature_matrix.iloc[1:].isna()).sum(axis=1) > threshold_feature
    ]
    feature_matrix.drop(rows_to_drop_feature, inplace=True)

    # Filter rows with too many zeros or NaNs in label matrix
    rows_to_drop_label = label_matrix.iloc[1:].index[
        ((label_matrix.iloc[1:] == 0) | label_matrix.iloc[1:].isna()).sum(axis=1) > threshold_label
    ]
    label_matrix.drop(rows_to_drop_label, inplace=True)

    # Reattach the header row
    feature_matrix = pd.concat([feature_matrix_header.to_frame().T, feature_matrix], ignore_index=True)
    label_matrix = pd.concat([label_matrix_header.to_frame().T, label_matrix], ignore_index=True)

    # Save the filtered matrices
    feature_matrix.to_csv(f"{output_dir}/filtered_feature_matrix.csv", index=False)
    label_matrix.to_csv(f"{output_dir}/filtered_label_matrix.csv", index=False)

    print(f"Filtered matrices saved to {output_dir}")

def main():
    parser = argparse.ArgumentParser(description="Filter raw TCGA matrices based on thresholds for zeros and NaNs.")
    parser.add_argument("--feature-matrix", required=True, help="Path to the raw feature matrix CSV file.")
    parser.add_argument("--label-matrix", required=True, help="Path to the raw label matrix CSV file.")
    parser.add_argument("--threshold-feature", type=int, default=1050, help="Maximum allowed zeros or NaNs in feature matrix rows.")
    parser.add_argument("--threshold-label", type=int, default=10000, help="Maximum allowed zeros or NaNs in label matrix rows.")
    parser.add_argument("--output-dir", required=True, help="Directory to save the filtered matrices.")

    args = parser.parse_args()

    filter_matrices(
        feature_matrix_path=args.feature_matrix,
        label_matrix_path=args.label_matrix,
        threshold_feature=args.threshold_feature,
        threshold_label=args.threshold_label,
        output_dir=args.output_dir,
    )

if __name__ == "__main__":
    main()
