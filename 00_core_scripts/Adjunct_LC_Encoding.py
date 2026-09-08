import numpy as np
import rasterio
from tqdm import tqdm


def process_year(year):
    edge_file = f"Users\\hyou\\mnt\\Hangkai\\CONUS_Forest_Edge_LCMAP\\LCMAP_edges\\LCMAP_{year}_edges.tif"
    landcover_file = f"Users\\hyou\\mnt\\Public\\LCMAP\\LCMAP_CU_{year}_V13_LCPRI.tif"
    output_file = f"Users\\hyou\\mnt\\Hangkai\\CONUS_Forest_Edge_LCMAP\\Edge_adjunct_LC\\Adjunct_LC_{year}.tif"

    with rasterio.open(edge_file) as edge_src, rasterio.open(landcover_file) as lc_src:
        edge_array = edge_src.read(1)
        lc_array = lc_src.read(1)
        rows, cols = lc_array.shape
        encoded_final_map = np.zeros_like(lc_array, dtype=np.uint16)

        for i in tqdm(range(1, rows - 1)):
            for j in range(1, cols - 1):
                code = edge_array[i, j]
                if code > 0:  # It's an edge pixel
                    top = lc_array[i - 1, j] if code & 0b1000 else 0
                    bottom = lc_array[i + 1, j] if code & 0b0100 else 0
                    left = lc_array[i, j - 1] if code & 0b0010 else 0
                    right = lc_array[i, j + 1] if code & 0b0001 else 0
                    encoded_final_map[i, j] = 1 * 10000 + top * 1000 + bottom * 100 + left * 10 + right

        # Save the final encoded map
        with rasterio.open(output_file, 'w', driver='GTiff', height=encoded_final_map.shape[0],
                           width=encoded_final_map.shape[1], count=1, dtype='uint16',
                           crs=lc_src.crs, transform=lc_src.transform) as dst:
            dst.write(encoded_final_map, 1)


# Process each year from 1985 to 2021
for year in (range(1985, 1995)):
    process_year(year)
