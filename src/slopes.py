import concurrent.futures
import multiprocessing
import threading

import numpy as np
import rasterio
import rasterio.warp

NODATA = -1

def slope_thresholds_at_full_resolution(
    path_to_eudem_slope_data, max_slope, outpath_to_tech_slope_limit, threads
):
    """
    EUDEM slope data is given in integer values, requiring a translation to degrees to
    then assign pixels either a 1 (slope is low enough to allow technology deployment)
    or 0 (slope is too steep to allow technology deployment).

    Capitalises on the possible concurrency of rasterio windows to improve memory and time efficiency,
    see https://rasterio.readthedocs.io/en/latest/topics/concurrency.html
    """
    _lim = _deg_to_int(max_slope)

    print(
        f"Max slope limit of {max_slope} translated from degrees to integer slope limit of "
        f"{_lim} corresponding to slope limits of {_int_to_deg(_lim):.2f} degrees"
    )

    with rasterio.open(path_to_eudem_slope_data, 'r') as src:
        profile = src.profile.copy()
        profile.update(
            compress='lzw',
            dtype=rasterio.float32,
            nodata=NODATA
        )
        with rasterio.open(outpath_to_tech_slope_limit, 'w', **profile) as dst:
            windows = [window for ij, window in dst.block_windows(1)]
            # We cannot write to the same file from multiple threads
            # without causing race conditions. To safely read/write
            # from multiple threads, we use a lock to protect the
            # DatasetReader/Writer
            read_lock = threading.Lock()
            write_lock = threading.Lock()
            def process(window):
                with read_lock:
                    src_array = src.read(1, window=window)
                result = _get_valid_pixels_from_tech_slope_limit(src_array, _lim)
                with write_lock:
                    dst.write(result, 1, window=window)
            # We map the process() function over the list of windows.
            with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
                executor.map(process, windows)


def _deg_to_int(deg):
    return int(round(np.cos(deg / (180 / np.pi)) * 250, 0))


def _int_to_deg(r):
    return np.arccos(r / 250) * 180 / np.pi


def _get_valid_pixels_from_tech_slope_limit(r, lim):
        _arr = np.where(r <= lim, 0, 1)  # steeper = lower integer value
        return np.where(r == 0, NODATA, _arr).astype(np.float32)  # integer value of zero = 90 degrees (i.e. NaN)


if __name__ == "__main__":
    slope_thresholds_at_full_resolution(
        path_to_eudem_slope_data=snakemake.input.slopes_in_europe,
        max_slope=snakemake.params.max_slope,
        threads=snakemake.params.max_threads,
        outpath_to_tech_slope_limit=snakemake.output[0]
    )
