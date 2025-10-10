import os
import pandas as pd
import requests

# Input and output paths
input_csv = "../Data/mouse_projection_data_sets.csv"
output_dir = "../Data/meso_projection"
os.makedirs(output_dir, exist_ok=True)

# Read dataset IDs
df = pd.read_csv(input_csv)
if "data_set_id" not in df.columns:
    raise ValueError("CSV must contain a column named 'data_set_id'")
if "projection_structure_unionizes_file_url" not in df.columns:
    raise ValueError("CSV must contain a column named 'projection_structure_unionizes_file_url'")

for _, row in df.iterrows():
    exp_id = int(row["data_set_id"])
    url = row["projection_structure_unionizes_file_url"]
    output_file = os.path.join(output_dir, f"experiment_{exp_id}.csv")

    try:
        print(f"Downloading {url} -> {output_file}")
        response = requests.get(url, stream=True)
        response.raise_for_status()

        with open(output_file, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)

    except Exception as e:
        print(f"⚠️ Failed to download {exp_id}: {e}")
