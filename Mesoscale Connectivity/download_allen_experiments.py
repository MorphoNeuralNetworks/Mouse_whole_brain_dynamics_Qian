import os
import pandas as pd
import requests

# Input and output paths
input_csv = "/mnt/data/mouse_projection_data_sets.csv"
output_dir = "../Data/meso_projection"
os.makedirs(output_dir, exist_ok=True)

# Read dataset IDs
df = pd.read_csv(input_csv)
if "data_set_id" not in df.columns:
    raise ValueError("CSV must contain a column named 'data_set_id'")

# Base URL template (Allen Brain Atlas projection_structure)
url_template = "http://connectivity.brain-map.org/projection/experiment/{exp_id}.csv"

for exp_id in df["data_set_id"].dropna().astype(int).unique():
    url = url_template.format(exp_id=exp_id)
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
