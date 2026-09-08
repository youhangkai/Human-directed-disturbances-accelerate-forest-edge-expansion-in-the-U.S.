import rasterio
import numpy as np
import gc  # Import the garbage collector
from tqdm import tqdm

def save_tif(output_path, data, transform, crs):
    """Saves a multi-band TIFF file where each band represents edge dynamics."""
    with rasterio.open(
        output_path, 'w',
        driver='GTiff',
        height=data.shape[0],
        width=data.shape[1],
        count=data.shape[2],  # Number of bands
        dtype=data.dtype,
        crs=crs,
        transform=transform
    ) as dst:
        for i in range(data.shape[2]):
            dst.write(data[:, :, i], i + 1)  # Write each band separately

def read_tif(path):
    with rasterio.open(path) as src:
        return src.read(1)  # Read the first band

def convert_to_binary(lc_map, forest_val=4):
    """Convert land cover map to binary map where forest pixels are True."""
    return lc_map == forest_val

def classify_edge_dynamics(edge_indices,prev_edge_map, curr_edge_map, prev_forest_map, curr_forest_map):
    # Define the dynamics map with dimensions for each edge of the pixel (north, south, east, west)
    dynamics_map = np.zeros((*prev_edge_map.shape, 4), dtype=np.uint8)

    # Bit positions for north, south, east, west
    positions = {
        'north': 0b1000,
        'south': 0b0100,
        'east': 0b0001,
        'west': 0b0010
    }


    # Process only pixels that have edges in either year
    for i, j in tqdm(edge_indices):
        if 1 <= i < prev_edge_map.shape[0] - 1 and 1 <= j < prev_edge_map.shape[1] - 1:
            for k, pos in positions.items():
                prev_edge = prev_edge_map[i, j] & pos
                curr_edge = curr_edge_map[i, j] & pos
                prev_forest = prev_forest_map[i, j]
                curr_forest = curr_forest_map[i, j]

                # Classify dynamics based on edge changes and forest changes
                if prev_edge and not curr_edge:
                    if not curr_forest:
                        dynamics_map[i, j, list(positions).index(k)] = 1  # Edge decrease due to forest decrease
                    elif curr_forest:
                        dynamics_map[i, j, list(positions).index(k)] = 2  # Edge decrease due to forest increase
                elif not prev_edge and curr_edge:
                    if not prev_forest:
                        dynamics_map[i, j, list(positions).index(k)] = 4  # Edge increase due to forest increase
                    elif prev_forest:
                        dynamics_map[i, j, list(positions).index(k)] = 3  # Edge increase due to forest decrease
                elif prev_edge and curr_edge:
                    dynamics_map[i, j, list(positions).index(k)] = 5  # Stable edge

    return dynamics_map

# Paths and processing as previously defined
years = range(1985, 2021)
forest_val = 4

for year in years:
    print(f"Processing year {year}")
    prev_edge_path = f"/Users/hyou/mnt/Hangkai/CONUS_Forest_Edge_LCMAP/LCMAP_edges/LCMAP_{year}_edges.tif"
    curr_edge_path = f"/Users/hyou/mnt/Hangkai/CONUS_Forest_Edge_LCMAP/LCMAP_edges/LCMAP_{year + 1}_edges.tif"
    prev_lc_path = f"/Users/hyou/mnt/Public/LCMAP/LCMAP_CU_{year}_V13_LCPRI.tif"
    curr_lc_path = f"/Users/hyou/mnt/Public/LCMAP/LCMAP_CU_{year + 1}_V13_LCPRI.tif"
    output_path = f"/Users/hyou/mnt/Hangkai/CONUS_Forest_Edge_LCMAP/LCMAP_edge_dynamics/LCMAP_{year}_{year + 1}_edge_dynamics.tif"

    prev_lc_map = read_tif(prev_lc_path)
    prev_forest_map = convert_to_binary(prev_lc_map, forest_val)
    del prev_lc_map
    gc.collect()
    curr_lc_map = read_tif(curr_lc_path)
    curr_forest_map = convert_to_binary(curr_lc_map, forest_val)
    del curr_lc_map
    gc.collect()

    prev_edge_map = read_tif(prev_edge_path)
    curr_edge_map = read_tif(curr_edge_path)

    edges_present = (prev_edge_map > 0) | (curr_edge_map > 0)
    edge_indices = np.argwhere(edges_present)
    del edges_present
    gc.collect()

    dynamics_map = classify_edge_dynamics(edge_indices,prev_edge_map, curr_edge_map, prev_forest_map, curr_forest_map)

    with rasterio.open(prev_edge_path) as src:
        save_tif(output_path, dynamics_map, src.transform, src.crs)


