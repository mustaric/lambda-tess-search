"""Modify the amount of sub-cube you wish to extract in the
``fits_mmap()`` function. Then run this::

    python measure_memory.py

My results below. This shows that memory mapping should only
consume the memory used in accessed sub-cube, not the entire cube.
Trying to calculate the mean of the entire cube did not work,
as expected.

Line #    Mem usage    Increment   Line Contents
================================================
    40     41.8 MiB     41.8 MiB   @profile
    41                             def fits_mmap():
    42     41.8 MiB      0.0 MiB       with fits.open(FILENAME, memmap=True) as pf:
    43     47.6 MiB      5.8 MiB           data = pf[1].data[0:10, 0:10, 0, 0]
    44     48.1 MiB      0.4 MiB           print(data.mean())
    45     48.1 MiB      0.0 MiB       return data

Line #    Mem usage    Increment   Line Contents
================================================
    40     41.8 MiB     41.8 MiB   @profile
    41                             def fits_mmap():
    42     41.8 MiB      0.0 MiB       with fits.open(FILENAME, memmap=True) as pf:
    43     47.6 MiB      5.8 MiB           data = pf[1].data[0:100, 0:100, :, 0]
    44    145.8 MiB     98.2 MiB           print(data.mean())
    45    145.8 MiB      0.0 MiB       return data

"""  # noqa
from astropy.io import fits

from memory_profiler import profile

# Change this to point to your local cube file.
# shape = (2078, 2136, 1282, 2)
FILENAME = '../../../tess-s0001-1-1-cube.fits'


@profile
def fits_mmap():
    with fits.open(FILENAME, memmap=True) as pf:
        data = pf[1].data[0:100, 0:100, :, 0]
        print(data.mean())
    return data


if __name__ == '__main__':
    fits_mmap()
