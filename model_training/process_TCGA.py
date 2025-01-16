import os
import pandas as pd

# Paths
manifest_file = "data/training_data/TCGA_barcodes.csv"
features_dir = "data/TCGA_data/features/"
labels_dir = "data/TCGA_data/labels/"
output_feature_matrix = "data/raw_feature_matrix.csv"
output_label_matrix = "data/raw_label_matrix.csv"

# Load the manifest file
manifest = pd.read_csv(manifest_file)

# Function to find the first .tsv or .txt file in a given folder
def find_first_file(base_dir, uuid, file_extension):
    uuid_path = os.path.join(base_dir, uuid)
    if os.path.isdir(uuid_path):
        files = [f for f in os.listdir(uuid_path) if f.endswith(file_extension)]
        if files:
            return os.path.join(uuid_path, files[0])  # Return the first matching file
    print(f"Warning: No {file_extension} file found for UUID {uuid}")
    return None

# Function to create the feature matrix
def create_feature_matrix(features_dir, manifest):
    matrix = None
    for _, row in manifest.iterrows():
        uuid = row["feature"]
        barcode = row["barcode"]
        file_path = find_first_file(features_dir, uuid, ".tsv")
        if not file_path:
            continue

        # Read the second row as the header
        raw_header = pd.read_csv(file_path, sep="\t", skiprows=1, nrows=1)
        new_header = raw_header.columns  # Extract column names

        # Read the data starting from row 6 (after skipping 4 rows of data)
        data = pd.read_csv(file_path, sep="\t", skiprows=6, names=new_header, usecols=["gene_name", "tpm_unstranded"])
        data = data.rename(columns={"tpm_unstranded": barcode})

        # Combine columns
        if matrix is None:
            matrix = data
        else:
            matrix = pd.concat([matrix, data[barcode]], axis=1)

    # Transpose the matrix
    matrix = matrix.set_index("gene_name").T.reset_index()
    matrix.rename(columns={"index": "sample"}, inplace=True)
    return matrix

# Function to create the label matrix
def create_label_matrix(labels_dir, manifest):
    matrix = None
    for _, row in manifest.iterrows():
        uuid = row["label"]
        barcode = row["barcode"]
        file_path = find_first_file(labels_dir, uuid, ".txt")
        if not file_path:
            continue

        # Read and extract relevant columns
        data = pd.read_csv(file_path, sep="\t", usecols=["miRNA_ID", "reads_per_million_miRNA_mapped"])
        data = data.rename(columns={"reads_per_million_miRNA_mapped": barcode})

        # Combine columns
        if matrix is None:
            matrix = data
        else:
            matrix = pd.concat([matrix, data[barcode]], axis=1)

    # Transpose the matrix
    matrix = matrix.set_index("miRNA_ID").T.reset_index()
    matrix.rename(columns={"index": "sample"}, inplace=True)
    return matrix

# Create the feature matrix
print("Creating feature matrix...")
feature_matrix = create_feature_matrix(features_dir, manifest)

# Create the label matrix
print("Creating label matrix...")
label_matrix = create_label_matrix(labels_dir, manifest)

# Save matrices to CSV files
print("Saving matrices...")
feature_matrix.to_csv(output_feature_matrix, index=False)
label_matrix.to_csv(output_label_matrix, index=False)

print(f"Feature matrix saved to {output_feature_matrix}")
print(f"Label matrix saved to {output_label_matrix}")
