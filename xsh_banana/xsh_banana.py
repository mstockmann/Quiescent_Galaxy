__author__ = 'Johannes Zabl (jz)'
__copyright__ = '2016 jz'
__email__ = 'johannes@observingtheuniverse.com'

import numpy as np
import astropy.io.fits as fits
import numpy.ma as ma

def remove_median(data, mask,axis=0, minel=25):

    arr1 = ma.masked_array(data, mask=~mask) # Note that mask is currently inverted, as I had defined the mask in a positive sense


    n_elements_ok = np.sum(mask, axis=axis)
    to_few_elements_mask = n_elements_ok < minel

    inter = ma.median(arr1, axis=axis)
    inter2 = inter.data.copy()
    mask_not_ok = inter.mask | to_few_elements_mask
    inter2[mask_not_ok] = np.median(inter2[~mask_not_ok])
    inter = inter2

    if axis == 1:
        inter = inter[:, np.newaxis] * np.ones(data.shape[1])

    data_new = data - inter
    return data_new

def correct_dark_variations(fnames, fname_mask_good, fname_mask_ignore, extra_corr=False):
    """

    :param fnames: Must be a list of X-Shooter raw science frames
    :param fname_mask_good: Defines the region, based on which the median is determined
    :param fname_mask_ignore: Defines a region, for which no correction is done;
                             (important, as otherwise part of some order is corrected)
    :param extra_corr: default False (currently, no other option should be chosen);
                        plan is to get the correction also redwards of the currently corrected orders implement

    :return:
    Warning: Automatically files are created in the directory of the raw data, with the ending '.median_back.fits'
    The function returns a list with the updated file names; this list can be directly
    used in Martin's pipeline manager
    """
    
    mask_good = fits.getdata(fname_mask_good)
    mask_good = mask_good.astype(bool)

    mask_ignore = fits.getdata(fname_mask_ignore)
    mask_ignore = mask_ignore.astype(bool)

    fnames_out = []

    for fname in fnames:
        data = fits.getdata(fname)
        hdu = fits.open(fname)
        data_new = remove_median(data, mask_good, axis=1)
        fname_out = fname.replace('.fits', '.median_back.fits')
        data_new[mask_ignore] = data[mask_ignore]
        if extra_corr == True: # FIXME This is currently not correctly working
            median_row_extra = np.median(data_new[556, 4:49])
            median_surrounding = np.median(data_new[np.array([555,557]),4:49])
            median_diff = median_row_extra - median_surrounding
            data_new[556, :] -= median_diff

        hdu[0].data = data_new


        hdu.writeto(fname_out, output_verify='fix', clobber=True)
        fnames_out.append(fname_out)

    return fnames_out
