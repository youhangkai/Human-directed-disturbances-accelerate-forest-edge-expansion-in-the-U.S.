import rasterio
import numpy as np
import os
import pandas as pd
from tqdm import trange

def calculate_forest_edge_statistics_by_year(directory_path, start_year=1985, end_year=2021):
    results = []

    for year in trange(start_year, end_year + 1):
        file_path = os.path.join(directory_path, f'LCMAP_{year}_edges.tif')
        if os.path.exists(file_path):
            with rasterio.open(file_path) as src:
                data = src.read(1)
                # Calculate pixel counts for each category
                pixel_counts = {i: np.sum(data == i) for i in range(1, 16)}

                results.append({'year': year, **pixel_counts})

    # Convert the list of dictionaries to a DataFrame
    df_results = pd.DataFrame(results)
    return df_results

# Usage
directory_path = '/Users/hyou/mnt/Hangkai/CONUS_Forest_Edge_LCMAP/LCMAP_edges'
forest_edge_stats_by_year = calculate_forest_edge_statistics_by_year(directory_path)

# Save the results to CSV
csv_path = os.path.join(directory_path, 'forest_edge_statistics.csv')
forest_edge_stats_by_year.to_csv(csv_path, index=False)
print(f"Results saved to {csv_path}")

