import geopandas as gpd
import rasterio
from rasterio.mask import mask
import numpy as np
import pandas as pd
from tqdm import tqdm

# Load the Shapefile
shapefile_path = '/Users/hyou/mnt/Hangkai/CONUS Forest Edge Mapping/CONUS shapefile/CONUS_ECOSYSTEM.shp'
shapes = gpd.read_file(shapefile_path)


def clip_and_count_values(raster, geometry, name, na_l1name):
    """Clip raster with geometry and count unique pixel values greater than 0."""
    out_image, out_transform = mask(raster, [geometry], crop=True)
    unique, counts = np.unique(out_image[out_image > 0], return_counts=True)
    result = {'NAME': name, 'NA_L1NAME': na_l1name}  # Initialize with name fields
    counts_dict = dict(zip(unique, counts))
    result.update(counts_dict)  # Add count results
    return result


# Iterate over years
for year in range(1986, 2021):
    print(year)
    raster_path = f'/Users/hyou/mnt/Hangkai/CONUS_Forest_Edge_LCMAP/Edge_adjunct_LC/Adjunct_LC_{year}.tif'
    raster = rasterio.open(raster_path)

    # Ensure CRS consistency
    if shapes.crs != raster.crs:
        shapes = shapes.to_crs(raster.crs)

    results = []

    # Process each polygon in the shapefile
    for index, row in tqdm(shapes.iterrows(), total=shapes.shape[0], desc=f'Processing {year}'):
        result = clip_and_count_values(raster, row['geometry'], row['NAME'], row['NA_L1NAME'])
        result['shape_id'] = index  # Add shape identifier
        results.append(result)

    # Create DataFrame and save to CSV for the current year
    df = pd.DataFrame(results)
    df.to_csv(
        f'/Users/hyou/mnt/Hangkai/CONUS_Forest_Edge_LCMAP/Edge_adjunct_LC/Ecoregion_Classification/eco_{year}.csv',
        index=False)
    print(f"Data for {year} processed and saved.")

    raster.close()  # Close the raster file
