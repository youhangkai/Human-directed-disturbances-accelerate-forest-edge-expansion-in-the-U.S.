import os
import numpy as np
import rasterio
from scipy.ndimage import binary_dilation
from tqdm import tqdm


def classify_forest_depth(input_filepath, output_folder):
    # Make sure the output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # Reading input forest cover data
    with rasterio.open(input_filepath) as src:
        data = src.read(1)
        meta = src.meta
        # Forest data: True for forest pixels, False for non-forests
        forest_data = data == 4
        del data  # Delete the data array and free up memory

    # Create an output array to store the depth markers for each pixel
    depth_layers = np.zeros_like(forest_data, dtype=np.uint8)

    # Set the structure element to define the neighborhood of the pixel
    structure = np.array([[0, 1, 0],
                          [1, 1, 1],
                          [0, 1, 0]])

    # Check the depth of each forest pixel
    for i in tqdm(range(forest_data.shape[0])):
        for j in range(forest_data.shape[1]):
            if forest_data[i, j]:  # Processing of forest pixels only
                found_non_forest = False
                for d in range(1, 5):  # Check 1 to 4 pixel distance
                    # Use binary_dilation to increase the scope of checking
                    mask = binary_dilation(np.array([[1]]), structure, iterations=d)
                    min_i, max_i = max(0, i - d), min(forest_data.shape[0], i + d + 1)
                    min_j, max_j = max(0, j - d), min(forest_data.shape[1], j + d + 1)
                    # Get the data in the surrounding d-grid
                    neighborhood = forest_data[min_i:max_i, min_j:max_j]
                    if not mask.shape == neighborhood.shape:
                        mask = mask[:neighborhood.shape[0], :neighborhood.shape[1]]
                    # Check for the presence of non-forest pixels
                    if np.any(~neighborhood & mask):
                        depth_layers[i, j] = d
                        found_non_forest = True
                        break
                if not found_non_forest:
                    depth_layers[i, j] = 5

    # Setting non-forested areas to 0
    depth_layers[~forest_data] = 0

    output_filepath = os.path.join(output_folder, os.path.basename(input_filepath))

    # write out data
    with rasterio.open(output_filepath, 'w', **meta) as dst:
        dst.write(depth_layers, 1)


input_folder = '/Users/hyou/mnt/Public/LCMAP'
output_folder = '/Users/hyou/mnt/Hangkai/CONUS_Forest_Edge_LCMAP/Forest_Depth_Classification'
for year in range(2005, 2007):
    print(year)
    input_filepath = os.path.join(input_folder, f'LCMAP_CU_{year}_V13_LCPRI.tif')
    classify_forest_depth(input_filepath, output_folder)
