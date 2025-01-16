import requests
import csv
import json
import re
import os
import argparse

# GDC API endpoint
data_endpt = "https://api.gdc.cancer.gov/data"

# Function to download data
def download_data(uuid_list, output_dir, data_type):
    params = {"ids": uuid_list}
    response = requests.post(
        data_endpt,
        data=json.dumps(params),
        headers={"Content-Type": "application/json"},
    )
    
    # Extract filename from headers
    response_head_cd = response.headers.get("Content-Disposition", "")
    file_name_match = re.findall(r"filename=(.+)", response_head_cd)
    file_name = file_name_match[0] if file_name_match else f"{data_type}_data.tar.gz"

    # Save the file
    file_path = os.path.join(output_dir, file_name)
    with open(file_path, "wb") as output_file:
        output_file.write(response.content)
    
    print(f"Downloaded {file_name} to {output_dir}")

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Download data using UUIDs from a CSV file.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input CSV file containing UUIDs.")
    args = parser.parse_args()

    # Paths
    input_file = args.input
    base_dir = "/data/TCGA_data"  # Base folder for storing downloaded data
    output_dir_features = os.path.join(base_dir, "features")
    output_dir_labels = os.path.join(base_dir, "labels")

    # Create folder structure if not exists
    os.makedirs(output_dir_features, exist_ok=True)
    os.makedirs(output_dir_labels, exist_ok=True)

    # Read the CSV and collect UUIDs
    with open(input_file, "r") as file:
        csv_reader = csv.DictReader(file)
        feature_uuids = []
        label_uuids = []
        
        for row in csv_reader:
            feature_uuids.append(row["feature"])
            label_uuids.append(row["label"])

    # Download features in batches (GDC API supports limited UUIDs per request)
    batch_size = 100
    for i in range(0, len(feature_uuids), batch_size):
        batch = feature_uuids[i:i+batch_size]
        download_data(batch, output_dir_features, "features")

    # Download labels in batches
    for i in range(0, len(label_uuids), batch_size):
        batch = label_uuids[i:i+batch_size]
        download_data(batch, output_dir_labels, "labels")

if __name__ == "__main__":
    main()
