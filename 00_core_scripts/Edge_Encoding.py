import numpy as np
import rasterio
from tqdm import tqdm

def get_forest_edge_encoding(landcover_map):
    with rasterio.open(landcover_map) as src:
        lc_array = src.read(1)

        edge_encoding = np.zeros_like(lc_array, dtype=np.uint8)

        forest_val = 4

        rows, cols = lc_array.shape

        for i in (range(1, rows - 1)):
            for j in range(1, cols - 1):
                if lc_array[i, j] == forest_val:
                    code = 0b0000
                    if lc_array[i - 1, j] != forest_val:
                        code |= 0b1000
                    if lc_array[i + 1, j] != forest_val:
                        code |= 0b0100
                    if lc_array[i, j - 1] != forest_val:
                        code |= 0b0010
                    if lc_array[i, j + 1] != forest_val:
                        code |= 0b0001
                    edge_encoding[i, j] = code

    return edge_encoding

def save_encoding_map(encoding_map, input_tif, output_tif):
    with rasterio.open(input_tif) as src:
        profile = src.profile
        profile.update(
            dtype=rasterio.uint8,
            count=1,
            compress='deflate'
        )

        with rasterio.open(output_tif, 'w', **profile) as dst:
            dst.write(encoding_map, 1)

input_tif = "/mnt/cephfs/scratch/groups/chen_group/hangkai/CONUS/Landcover/LCMAP_CU_2020_V13_LCPRI.tif"
output_tif = "/mnt/cephfs/scratch/groups/chen_group/hangkai/CONUS/LCMAP_edges/LCMAP_2020_edges.tif"
print(input_tif)
print(output_tif)

encoding_map = get_forest_edge_encoding(input_tif)
save_encoding_map(encoding_map, input_tif, output_tif)
print("Edge encoding map saved successfully.")
