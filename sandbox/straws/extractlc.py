import numpy as np

#from pdb import set_trace as debug
from astroquery.mast import Catalogs, Tesscut
from astropy.coordinates import SkyCoord



"""
Quickest, Dirtiest possible lightcurve extraction library.

The goal is to go end to end from Tess cutout to folded lightcurve
in as short a development time as possible

Improved algorithms will be subbed in as necessary.

Sample usage:
cube, hdr = getWasp126()
main(cube, hdr)

"""

__version__ = "$Id: extractlc.py 37 2019-03-14 01:46:27Z mustaric $"
__URL__ = "$URL: https://svn.code.sf.net/p/greataunttess/code/trunk/extractlc.py $"


def getWasp126():
    ticId = 25155310
    sector = 1
    size_pixels = 20

    time, cube, hdr = loadSingleSector(ticId, sector, size_pixels)
    return time, cube, hdr


def main(cube, hdr):
    sky = measureSky(cube)

    centroid = np.array(cube.shape[1:]) / 2.
    aper = chooseAperture(cube, 2 * sky, centroid)

    input()

    lc = performAperturePhotometry(cube, sky, aper)

    input()



def loadSingleSector(ticId, sector, size_pixels):
    starname = "TIC %08i" %(ticId)

    catalogData = Catalogs.query_object(starname, radius=1/60., catalog="TIC")
    ra = catalogData[0]['ra']
    dec = catalogData[0]['dec']
    coord = SkyCoord(ra, dec, unit='deg')

    data, hdr, wcshdr = getTessCutout(coord, size_pixels, sector)
    time = data['TIME']
    cube = getTargetPixelArrayFromFits(data, hdr)

    return time, cube, hdr, wcshdr


def getTessCutout(coord, size_pixels, sector):
    """
    Returns first returned for that sector.
    """

    hdulist = Tesscut.get_cutouts(coord, size_pixels, sector=sector)

    if len(hdulist)<1:
        raise IOError

    hdu = hdulist[0]
    return hdu[1].data, hdu[1].header, hdu[2].header



def getTargetPixelArrayFromFits(data, hdr, column='FLUX'):
    """Convert a FITS data structure to a cube

    cube.shape == (ncin, nc, nr)
    where ncin = "Number of cadences in sector, nc is number of columns
    in image, and nr is number of rows
    """

    nRows = int(hdr['NAXIS2'])
    #All imgs are the same size, so it doesn't matter which TDIM
    #I pick here.
    shape = eval(hdr['TDIM5'])
    #The reverse puts the columns on the x-axis and row
    #on the y-axis when you call mp.imshow()
    shape = tuple(reversed(shape))

    tpfArray = np.empty((nRows, shape[0], shape[1]))

    for i, cadence in enumerate(data[column]):
        tpfArray[i,:,:] = cadence.reshape(shape)

    return tpfArray


def chooseAperture(cube, skyThreshold, centroids, radius=1.5):
    """
    Pick an aperture using CCL algorithm and input sky threshold
    But if nothing is returned, then pick a circular aperture.
    Either way it expects the star to be at the center of the image.
    """
    try:
        aperture = chooseCclAperture(cube, skyThreshold, centroids)

    except ValueError:
        dim = np.shape(cube)[1:3]
        ap = chooseSimpleAperture(dim, radius=radius, starloc=centroids)
        aperture = ap == 1

    return aperture


def chooseCclAperture(cube, skyThreshold, centroids):
    """Choose an aperture using a connected component labelling
    algorithm
    """
    img = cube.sum(axis=0)
    assert np.all(np.isfinite(img))

    threshold = np.sum(skyThreshold)

    #@TODO Some sort of graph theory approach would be much
    #faster here
    labels = label_image(img, threshold)

    c0, r0 = centroids.astype(int)
    i0 = labels[c0, r0]
    if i0 == 0:
        raise ValueError("No star found in center of image in chooseCClAperture.")

    out = np.zeros_like(labels, dtype=bool)
    out[labels == i0] = True
    return out



def chooseSimpleAperture(dim, starloc=None, radius=1.5):
    """
    dim = (x,y) size of the image
    Choose a simple circular aperture of pixel radius .
    If starloc is None, assume star is in the center
    otherwise starloc is a tuple of the (x,y) indicies of where the star is located.

    Returns 2dimage where 1 indicates that it is in the aperture.
    """

    if starloc is None:
        x0 = np.floor(dim[0]/2)
        y0 = np.floor(dim[1]/2)
    else:
        x0 = starloc[0]
        y0 = starloc[1]
    assert x0 < dim[0]
    assert y0 < dim[1]

    aperture = np.ones(dim[::-1])


    for x in np.arange(0,dim[0]):
        for y in np.arange(0,dim[1]):

            dist = ( (x-x0)**2 + (y-y0)**2)**(1/2)

            if dist > radius:
                aperture[y,x] = 0

    return aperture


def performAperturePhotometry(cube, sky, aperture):
    cube = cube.copy()

    cube -= sky[:, None, None]
    cube *= aperture[None, :,:]
    tmp = np.sum(cube, axis=2)

    lc = np.sum(tmp, axis=1)
    assert np.all(np.isfinite(lc))
    return lc


def measureSky(cube):
    ncin = cube.shape[0]

    sky = np.zeros(ncin)
    for i in range(ncin):
        sky[i] = np.nanmedian(cube[i])

    return sky

def measureEmptySky(cube, skyFactor=1.9):
    """
     Remove bright pixels before final measure of the sky.
    Uses the same CCL algorithm to identify the bright pixels.
    """
    
    initsky = measureSky(cube)
    img = cube.sum(axis=0)
    threshold = skyFactor * np.sum(initsky)
    labels = label_image(img, threshold)
    want = labels == 0
    if len(want[want] > 1):
        sky = measureSky(cube[:,want])
    else:
        sky = initsky
        print("No sky pixels remain. Sky falling back to simple median of cutout.")
    return sky
    
    

def measureHistSky(cube, bins=np.arange(10,300,2)):
    ncin = cube.shape[0]
    sky = np.zeros(ncin)
    
    for i in range(ncin):
    
        h,b = np.histogram(cube[i,:,:].ravel(), bins=bins)
        bestbin = bins[np.argmax(h)]
        newbins = np.arange(bestbin,bestbin+1,(b[1]-b[0])/50)
        h,b = np.histogram(cube[i,:,:].ravel(), bins=newbins)
        sky[i] = newbins[np.argmax(h)]
    
    return sky

    
    
    
def generateFakeDataCube():
    cube = np.zeros( (1500, 20, 20))
    c0 = 10
    r0 = 9
    seeing_pix = float(.7)

    x, y = np.meshgrid( np.arange(20), np.arange(20) )

    zx  = (x - c0) / seeing_pix
    zy = (y - r0) / seeing_pix
    psf = 100 * np.exp(-.5 * zx**2) * np.exp(-.5 * zy**2)

    cube += psf

    cube += 10  #Some background
    sigma = np.sqrt(cube)

    poisson = np.random.randn(*cube.shape) * sigma
    cube += poisson
    return cube



def label_image(img, threshold):
    """
    Connecte component labeling (Rosenfeld 1966)
    """
    nr, nc = img.shape
    label = img * 0
    labelId = 1

    for i in range(nr):
        for j in range(nc):
            if img[i, j] > threshold:
                if j == 0 or label[i, j-1] == 0:
                    #Create a new labeled region
                    label[i, j] = labelId
                    labelId += 1
                else:
                    #Extend the current labeled region
                    label[i, j] = label[i, j-1]

    #At this point, we have a set of labeled regions, but each
    #is confined to its own column. Now we merged regions that touch
    #By scanning across the rows(instead of down the columns )
    for j in range(nc):
        for i in range(nr-1):
            v1 = label[i, j]
            v2 = label[i+1, j]

            if v1 == v2:
                continue

            if v1 > 0 and v2 > 0:
                idx = label == v2
                label[idx] = v1

    return label
