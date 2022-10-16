# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 23:15:44 2018
@author: user

https://nbviewer.org/github/ysbach/SNU_AOclass/blob/master/Notebooks/05-Differential_Phot.ipynb

"""
#%%
import os
from warnings import warn
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.time import Time
from astropy.visualization import ZScaleInterval, ImageNormalize
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats

from astroquery.jplhorizons import Horizons

from photutils.aperture import CircularAperture as CAp
from photutils.aperture import CircularAnnulus as CAn
from photutils.detection import DAOStarFinder
from photutils.psf.groupstars import DAOGroup
from astropy.nddata import CCDData, Cutout2D

from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib import gridspec

import Python_utilities

from astropy.table import Table, vstack
from astroquery.vizier import Vizier
from astroquery.mast import Catalogs
from scipy.interpolate import UnivariateSpline
#from xarray import Coordinate

#import ysfitsutilpy as yfu
#import ysphotutilpy as ypu
#import ysvisutilpy as yvu

import warnings
#%%

log_dir = "logs/"
log_file = "{}{}.log".format(log_dir, os.path.basename(__file__)[:-3])
err_log_file = "{}{}_err.log".format(log_dir, os.path.basename(__file__)[:-3])
print ("log_file: {}".format(log_file))
print ("err_log_file: {}".format(err_log_file))

base_dr = "../Post_processing/M35_Light_-_2018-10-31_-_TMB130ss_STF-8300M_-_1bin/"

### make all fits file list...
fullnames = Python_utilities.getFullnameListOfallFiles("{}/input".format(base_dr))
#print ("fullnames: {}".format(fullnames))
print ("len(fullnames): {}".format(len(fullnames)))

c_method = 'median'
master_dr = "master_files/"
reduced_dr = "readuced_files/"
result_dr = "result_files/"

if not os.path.exists('{0}'.format("{}{}".format(base_dr, result_dr))):
    os.makedirs("{}{}".format(base_dr, result_dr))
    print("{}{}is created".format(base_dr, result_dr))

fullnames_light = [w for w in fullnames \
            if ("_bias_" not in w.lower()) \
                and ("_dark_" not in w.lower()) \
                    and ("_flat_" not in w.lower())]
print ("len(fullnames_light): {}".format(len(fullnames_light)))

#%%
n = 0
for fullname in fullnames_light[:]:
    n += 1
    print('#'*40,
        "\n{2:.01f}%  ({0}/{1}) {3}".format(n, len(fullnames_light), 
                                            (n/len(fullnames_light))*100, os.path.basename(__file__)))
    print ("Starting...\nfullname: {}".format(fullname))

#%%
fullname = fullnames_light[0]
    
fullname_el = fullname.split("/")
hdul = fits.open(fullname)
hdr = hdul[0].header
data = hdul[0].data

img = hdul[0].data
print("img: {}".format(img))
print("img.shape: {}".format(img.shape))

# Set WCS and print for your information
w = WCS(hdr)
print("WCS: {}".format(w))


#%%
def panstarrs_query(ra_deg, dec_deg, rad_deg, columns=None, column_filters={},
                    maxsources=10000):
    """
    Query PanSTARRS @ VizieR using astroquery.vizier
    :param ra_deg: RA in degrees
    :param dec_deg: Declination in degrees
    :param rad_deg: field radius in degrees
    :param maxmag: upper limit G magnitude (optional)
    :param maxsources: maximum number of sources
    :return: astropy.table object
    Note
    ----
    All columns: http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=II/349
    """
    if columns is None:
        columns = ['objID', 'RAJ2000', 'DEJ2000', 'e_RAJ2000', 'e_DEJ2000',
                   'gmag', 'e_gmag', 'rmag', 'e_rmag', 'imag', 'e_imag', 
                   'zmag', 'e_zmag', 'ymag', 'e_ymag']
    vquery = Vizier(columns=columns,
                    column_filters=column_filters,
                    row_limit=maxsources)

    field = SkyCoord(ra=ra_deg, dec=dec_deg,
                     unit=(u.deg, u.deg),
                     frame='icrs')
    return vquery.query_region(field,
                               width=("{}d".format(rad_deg)),
                               catalog="II/349/ps1")[0]


image_center = np.array(hdul[0].shape) / 2 - 0.5
cent_coord = w.wcs_pix2world(image_center[0], image_center[1], 0)
# 0 means that we are using 0-indexing. You may use 1 for 1-indexing.

#result = panstarrs_query(cent_coord[0], cent_coord[1], rad_deg=0.1,
#                        column_filters={"gmag":"13.0..20.0", "e_gmag":"<0.10"})
#result
#%%
# Query object
objname = "Kleopatra"
observat = {"lon": 127.0, "lat": 37.5, "elevation": 101}
# Observatory code: see https://en.wikipedia.org/wiki/List_of_observatory_codes
# dict(lon=, lat=, elevation=)

t_obs = Time(hdr["DATE-OBS"]) + hdr["EXPTIME"] * u.s / 2  # middle of observation time

obj = Horizons(id=objname, location=observat, epochs=t_obs.jd)
obj_q = obj.ephemerides()

pos_sky = SkyCoord(obj_q["RA"][0], obj_q["DEC"][0], unit='deg')
pos_pix = pos_sky.to_pixel(wcs=w)

print("pos_sky, pos_pix: {}, {}".format(pos_sky, pos_pix))

#%%
def center_radec(
        ccd_or_header,
        center_of_image=True,
        ra_key="RA",
        dec_key="DEC",
        equinox=None,
        frame=None,
        equinox_key="EPOCH",
        frame_key="RADECSYS",
        ra_unit=u.hourangle,
        dec_unit=u.deg,
        mode="all",
        verbose=True,
        plain=False,
):
    """ Returns the central ra/dec from header or WCS.
    Notes
    -----
    Even though RA or DEC is in sexagesimal, e.g., "20 53 20", astropy
    correctly reads it in such a form, so no worries.
    Parameters
    ----------
    ccd_or_header : CCD-like, Header
        The ccd or header to extract the central RA/DEC from keywords or WCS.
    center_of_image : bool, optional
        If `True`, WCS information will be extracted from the ccd or header,
        rather than relying on the `ra_key` and `dec_key` keywords directly. If
        `False`, `ra_key` and `dec_key` from the header will be understood as
        the "center" and the RA, DEC of that location will be returned.
    equinox, frame : str, optional
        The `equinox` and `frame` for SkyCoord. Default (`None`) will use the
        default of SkyCoord. Important only if ``usewcs=False``.
    XX_key : str, optional
        The header key to find XX if ``XX`` is `None`. Important only if
        ``usewcs=False``.
    XX_unit : Quantity, optional
        The unit of ``XX``. Important only if ``usewcs=False``.
    mode : 'all' or 'wcs', optional
        Whether to do the transformation including distortions (``'all'``) or
        only including only the core WCS transformation (``'wcs'``). Important
        only if ``usewcs=True``.
    plain : bool
        If `True`, only the values of RA/DEC in degrees will be returned.
    """
    if isinstance(ccd_or_header, CCDData):
        header = ccd_or_header.header
        w = ccd_or_header.wcs
    elif isinstance(ccd_or_header, fits.Header):
        header = ccd_or_header
        w = WCS(header)

    if center_of_image:
        nx, ny = float(header["NAXIS1"]), float(header["NAXIS2"])
        centx = nx / 2 - 0.5
        centy = ny / 2 - 0.5
        coo = SkyCoord.from_pixel(centx, centy, wcs=w, origin=0, mode=mode)
    else:
        ra = get_from_header(header, ra_key, verbose=verbose)
        dec = get_from_header(header, dec_key, verbose=verbose)
        if equinox is None:
            equinox = get_from_header(header, equinox_key, verbose=verbose, default=None)
        if frame is None:
            frame = get_from_header(
                header, frame_key, verbose=verbose, default=None
            ).lower()
        coo = SkyCoord(
            ra=ra, dec=dec, unit=(ra_unit, dec_unit), frame=frame, equinox=equinox
        )

    if plain:
        return coo.ra.value, coo.dec.value
    return coo

# Position of the telescope FOV center 
# (RA/DEC of the pixel at the center)
cent_coord = center_radec(ccd_or_header=hdr, center_of_image=True)
print("cent_coord: {}".format(cent_coord))

#%%

PS1_DR1_DEL_FLAG = [
    0,   # FEW: Used within relphot; skip star.
    1,   # POOR: Used within relphot; skip star.
    2,   # ICRF_QSO: object IDed with known ICRF quasar
    3,   # HERN_QSO_P60: identified as likely QSO, P_QSO >= 0.60
    5,   # HERN_RRL_P60: identified as likely RR Lyra, P_RRLyra >= 0.60
    7,   # HERN_VARIABLE: identified as a variable based on ChiSq
    8,   # TRANSIENT: identified as a non-periodic (stationary) transient
    9,   # HAS_SOLSYS_DET: at least one detection identified with a known SSO
    10,  # MOST_SOLSYS_DET: most detections identified with a known SSO
    23,  # EXT: extended in our data (eg, PS)
    24   # EXT_ALT: extended in external data (eg, 2MASS)
]

def xyinFOV(
        table,
        header=None,
        wcs=None,
        shape=None,
        ra_key='ra',
        dec_key='dec',
        unit=None,
        bezel=0,
        bezel_x=None,
        bezel_y=None,
        origin=0,
        mode='all',
        return_mask=False,
        col_x="x",
        col_y="y",
        verbose=1
):
    ''' Convert RA/DEC to pixel with rejection at bezels
    Parameters
    ----------
    table : astropy.table.Table or pandas.DataFrame
        The queried result table.
    header : astropy.io.fits.Header, optional
        The header to extract WCS information. One and only one of `header` and
        `wcs` must be given.
    wcs : astropy.wcs.WCS, optional
        The WCS to convert the RA/DEC to XY. One and only one of `header` and
        `wcs` must be given.
    shape : 2-element tuple, optional
        The shape of the image (used to reject stars at bezels). If `header` or
        `wcs` is given, it is ignored and ``NAXISi`` keywords will be used.
    ra_key, dec_key : str, optional
        The column names containing RA/DEC.
    unit : `~astropy.Quantity`, optional
        The unit of the RA/DEC. If not given, it will be inferred from the
        table (only if astropy's QTable.)
    bezel : int, float, array-like, optional
        The bezel size. If array-like, it should be ``(lower, upper)``. If only
        this is given and `bezel_x` and/or `bezel_y` is/are `None`, it/both
        will be replaced by `bezel`. If you want to keep some stars outside the
        edges, put negative values (e.g., ``-5``).
    bezel_x, bezel_y : int, float, 2-array-like, optional
        The bezel (border width) for x and y axes. If array-like, it should be
        ``(lower, upper)``. Mathematically put, only objects with center
        ``(bezel_x[0] + 0.5 < center_x) & (center_x < nx - bezel_x[1] - 0.5)``
        (similar for y) will be selected. If you want to keep some stars
        outside the edges, put negative values (e.g., ``-5``).
    origin : int, optional
       Whether to return 0 or 1-based pixel coordinates.
    mode: 'all' or 'wcs', optional
        Whether to do the transformation including distortions (``'all'``) or
        only including only the core WCS transformation (``'wcs'``).
    col_x, col_y : str, optional
        The column names for x and y positions.
    '''
    bezel = np.atleast_1d(bezel)
    if bezel.shape == (1,):
        bezel = np.tile(bezel, 2)
    elif bezel.shape != (2,):
        raise ValueError(f"bezel must have shape (1,) or (2,). Now {bezel.shape}")

    _tab = table.copy()
    N_old = len(_tab)
    if isinstance(table, pd.DataFrame):
        _tab = Table.from_pandas(table)
        # Better to unify as astropy Table, because of the unit in RA/DEC may
        # be missing if converted to pandas..
    elif not isinstance(table, Table):
        raise TypeError("table must be either astropy Table or pandas DataFrame.")

    if wcs is None:
        wcs = WCS(header)

    if unit is None:
        coo = SkyCoord(_tab[ra_key], _tab[dec_key])
    else:
        coo = SkyCoord(_tab[ra_key], _tab[dec_key], unit=unit)

    x, y = coo.to_pixel(wcs=wcs, origin=origin, mode=mode)
    _tab[col_x] = x
    _tab[col_y] = y

    if shape is not None:
        nx, ny = shape[1], shape[0]
    else:
        if header is not None:
            nx, ny = header["NAXIS2"], header["NAXIS1"]
        elif wcs is not None:
            nx, ny = wcs._naxis
        else:  # If none of them is available, no way to reject stars at bezels.
            if verbose >= 1:
                print("No way to reject stars at bezels: provide `header` or `shape`.")
            if return_mask:
                return _tab, None
            return _tab

    mask = bezel_mask(x, y, nx, ny, bezel=bezel, bezel_x=bezel_x, bezel_y=bezel_y)
    _tab.remove_rows(mask)

    N_new = len(_tab)
    if verbose >= 1:
        mask_str(N_new, N_old, f"{bezel}-pixel bezel")

    if return_mask:
        return _tab, mask

    return _tab

def mask_str(n_new, n_old, msg):
    dn = n_old - n_new
    print(f"{n_new:3d} objects remaining: {dn:3d} masked out of {n_old:3d} based on {msg:s}.")


def bezel_mask(
        xvals,
        yvals,
        nx,
        ny,
        bezel=(0, 0),
        bezel_x=None,
        bezel_y=None
):
    '''
    Parameters
    ----------
    xvals, yvals : array-like
        The x and y position values.
    nx, ny : int or float
        The number of x and y pixels (``NAXIS2`` and ``NAXIS1``, respectively
        from header).
    bezel : int, float, array-like, optional
        The bezel size. If array-like, it should be ``(lower, upper)``. If only
        this is given and ``bezel_x`` and/or ``bezel_y`` is/are ``None``,
        it/both will be replaced by ``bezel``. If you want to keep some stars
        outside the edges, put negative values (e.g., ``-5``).
    bezel_x, bezel_y : int, float, 2-array-like, optional
        The bezel (border width) for x and y axes. If array-like, it should be
        ``(lower, upper)``. Mathematically put, only objects with center
        ``(bezel_x[0] + 0.5 < center_x) & (center_x < nx - bezel_x[1] - 0.5)``
        (similar for y) will be selected. If you want to keep some stars
        outside the edges, put negative values (e.g., ``-5``).
    '''
    bezel = np.array(bezel)
    if len(bezel) == 1:
        bezel = np.repeat(bezel, 2)

    if bezel_x is None:
        bezel_x = bezel.copy()
    else:
        bezel_x = np.atleast_1d(bezel_x)
        if len(bezel_x) == 1:
            bezel_x = np.repeat(bezel_x, 2)

    if bezel_y is None:
        bezel_y = bezel.copy()
    else:
        bezel_y = np.atleast_1d(bezel_y)
        if len(bezel_y) == 1:
            bezel_y = np.repeat(bezel_y, 2)

    mask = ((xvals < bezel_x[0] + 0.5)
            | (yvals < bezel_y[0] + 0.5)
            | (xvals > (nx - bezel_x[1]) - 0.5)
            | (yvals > (ny - bezel_y[1]) - 0.5)
            )
    return mask


def fov_radius(header, unit=u.deg):
    """ Calculates the rough radius (cone) of the (square) FOV using WCS.
    Parameters
    ----------
    header: Header
        The header to extract WCS information.
    Returns
    -------
    radius: `~astropy.Quantity`
        The radius in degrees
    """
    w = WCS(header)
    nx, ny = float(header["NAXIS1"]), float(header["NAXIS2"])
    # Rough calculation, so use mode='wcs'
    c1 = SkyCoord.from_pixel(0, 0, wcs=w, origin=0, mode="wcs")
    c2 = SkyCoord.from_pixel(nx, 0, wcs=w, origin=0, mode="wcs")
    c3 = SkyCoord.from_pixel(0, ny, wcs=w, origin=0, mode="wcs")
    c4 = SkyCoord.from_pixel(nx, ny, wcs=w, origin=0, mode="wcs")

    # TODO: Can't we just do ``return max(r1, r2).to(unit)``???
    #   Why did I do this? I can't remember...
    #   2020-11-09 14:29:29 (KST: GMT+09:00) ysBach
    r1 = c1.separation(c3).value / 2
    r2 = c2.separation(c4).value / 2
    r = max(r1, r2) * u.deg
    return r.to(unit)

# Get the radius of the smallest circle which encloses all the pixels
rad = fov_radius(header=hdr, unit=u.deg)

print("rad: {}".format(rad))

#%%
class PanSTARRS1:
    def __init__(self, ra, dec, frame="icrs", radius=None, inner_radius=None,
                 width=None, height=None, columns=["**", "+_r"],
                 column_filters={}):
        """ Query PanSTARRS @ VizieR using astroquery.vizier
        Parameters
        ----------
        ra, dec, radius : float or `~astropy.Quantity`
            The central RA, DEC and the cone search radius. If not
            `~astropy.Quantity`, assuming it is in degrees unit.
        frame : str
            The frame for the coordinates.
            Default: "icrs"
        inner_radius : cfloat or `~astropy.Quantity`
            When set in addition to `radius`, the queried region becomes
            annular, with outer radius `radius` and inner radius
            `inner_radius`. If not `~astropy.Quantity`, assuming it is in
            degrees unit.
        width : convertible to `~astropy.coordinates.Angle`
            The width of the square region to query. If not
            `~astropy.Quantity`, assuming it is in degrees unit.
        height : convertible to `~astropy.coordinates.Angle`
            When set in addition to ``width``, the queried region becomes
            rectangular, with the specified ``width`` and ``height``. If not
            `~astropy.Quantity`, assuming it is in degrees unit.
        columns : list of str, {'*', '**'}, optional
            The columns to be retrieved. The special column ``"*"`` requests
            just the default columns of a catalog; ``"**"`` (Default) would
            request all the columns. For sorting, use ``"+"`` in front of the
            column name. See the documentation:
            https://astroquery.readthedocs.io/en/latest/vizier/vizier.html#specifying-keywords-output-columns-and-constraints-on-columns
        column_filters : dict, optional
            The column filters for astroquery.vizier.
            Example can be ``{"gmag":"13.0..20.0", "e_gmag":"<0.10"}``.
        Return
        ------
        queried : astropy.table object
            The queried result.
        Notes
        -----
        All columns: http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=II/349
        """
        _params = dict(
            ra=ra,
            dec=dec,
            radius=radius,
            inner_radius=inner_radius,
            width=width,
            height=height,
        )

        for k, v in _params.items():
            if v is None:
                continue
            if not isinstance(v, u.Quantity):
                warn(f"{k} is not astropy Quantity: Assuming deg unit")
                _params[k] = v * u.deg

        self.ra = _params["ra"]
        self.dec = _params["dec"]
        self.radius = _params["radius"]
        self.inner_radius = _params["inner_radius"]
        self.width = _params["width"]
        self.height = _params["height"]
        self.frame = frame

        if isinstance(columns, str):
            if columns in ['*', '**']:
                self.columns = [columns]
            else:
                raise ValueError("If columns is str, it must be one of ['*', '**']")
        else:
            self.columns = columns

        self.column_filters = column_filters

    def query(self):
        vquery = Vizier(columns=self.columns, column_filters=self.column_filters, row_limit=-1)

        field = SkyCoord(ra=self.ra, dec=self.dec, frame=self.frame)

        self.queried = vquery.query_region(field,
                                           radius=self.radius,
                                           inner_radius=self.inner_radius,
                                           width=self.width,
                                           height=self.height,
                                           catalog="II/349/ps1")[0]

        return self.queried

    def select_xyinFOV(self, header=None, wcs=None, bezel=0, mode='all'):
        ''' Convert RA/DEC to xy (add columns) with rejection at bezels.
        Parameters
        ----------
        header : astropy.io.fits.Header, optional
            The header to extract WCS information. One and only one of `header`
            and `wcs` must be given.
        wcs : astropy.wcs.WCS, optional
            The WCS to convert the RA/DEC to XY. One and only one of `header`
            and `wcs` must be given.
        bezel: int or float, optional
            The bezel size to exclude stars at the image edges. If you want to
            keep some stars outside the edges, put negative values (e.g.,
            ``-5``).
        mode: 'all' or 'wcs', optional
            Whether to do the transformation including distortions (``'all'``)
            or only including only the core WCS transformation (``'wcs'``).
        '''
        self.queried = xyinFOV(table=self.queried, header=header, wcs=wcs,
                               ra_key="RAJ2000", dec_key="DEJ2000",
                               bezel=bezel, origin=0, mode=mode)

    def drop_for_diff_phot(self, del_flags=PS1_DR1_DEL_FLAG,
                           drop_by_Kron=True):
        ''' Drop objects which are not good for differential photometry.
        Parameters
        ----------
        del_flags : list of int, None, optional
            The flags to be used for dropping objects based on ``"f_objID"`` of
            Pan-STARRS1 query. These are the powers of 2 to identify the flag
            (e.g., 2 means ``2**2`` or flag ``4``). See Notes below for each
            flag. Set it to `None` to keep all the objects based on
            ``"f_objID"``.
        drop_by_Kron : bool, optional
            If `True` (default), drop the galaxies based on the Kron magnitude
            criterion suggested by PS1 (which works good only if i <~ 21):
            https://outerspace.stsci.edu/display/PANSTARRS/How+to+separate+stars+and+galaxies
        Notes
        -----
        H15 means Hernitschek+ 2015ApJ...801...45H.
          * FEW (1) : Used within relphot; skip star.
          * POOR (2) : Used within relphot; skip star.
          * ICRF_QSO (4) : object IDed with known ICRF quasar (may have ICRF
            position measurement)
          * HERN_QSO_P60 (8) : identified as likely QSO (H15), P_QSO >= 0.60
          * HERN_QSO_P05 (16) : identified as possible QSO (H15), P_QSO >= 0.05
          * HERN_RRL_P60 (32) : identified as likely RR Lyra (H15), P_RRLyra >=
            0.60
          * HERN_RRL_P05 (64) : identified as possible RR Lyra (H15), P_RRLyra
            >= 0.05
          * HERN_VARIABLE (128) : identified as a variable based on ChiSq (H15)
          * TRANSIENT (256) : identified as a non-periodic (stationary)
            transient
          * HAS_SOLSYS_DET (512) : at least one detection identified with a
            known solar-system object (asteroid or other).
          * MOST_SOLSYS_DET (1024 (2^10)) : most detections identified with a
            known solar-system object (asteroid or other).
          * LARGE_PM (2048) : star with large proper motion
          * RAW_AVE (4096) : simple weighted average position was used (no IRLS
            fitting)
          * FIT_AVE (8192) : average position was fitted
          * FIT_PM (16384) : proper motion model was fitted
          * FIT_PAR (32768) : parallax model was fitted
          * USE_AVE (65536) : average position used (not PM or PAR)
          * USE_PM (131072) : proper motion used (not AVE or PAR)
          * USE_PAR (262144) : parallax used (not AVE or PM)
          * NO_MEAN_ASTROM (524288) : mean astrometry could not be measured
          * STACK_FOR_MEAN (1048576 (2^20)) : stack position used for mean
            astrometry
          * MEAN_FOR_STACK (2097152) : mean astrometry used for stack position
          * BAD_PM (4194304) : failure to measure proper-motion model
          * EXT (8388608) : extended in our data (eg, PS)
          * EXT_ALT (16777216) : extended in external data (eg, 2MASS)
          * GOOD (33554432) : good-quality measurement in our data (eg,PS)
          * GOOD_ALT (67108864) : good-quality measurement in external data
            (eg, 2MASS)
          * GOOD_STACK (134217728) : good-quality object in the stack (>1 good
            stack measurement)
          * BEST_STACK (268435456) : the primary stack measurements are the
            best measurements
          * SUSPECT_STACK (536870912) : suspect object in the stack (no more
            than 1 good measurement, 2 or more suspect or good stack
            measurement)
          * BAD_STACK (1073741824 (2^30)) : poor-quality stack object (no more
            than 1 good or suspect measurement)
        Among the ``"f_objID"``, the following are better to be dropped because
        they are surely not usable for differential photometry:
            * 1, 2, 4, 8, 32, 128, 256, 512, 1024, 8388608, 16777216
        or in binary position (``del_flags``),
            * 0, 1, 2, 3, 5, 7, 8, 9, 10, 23, 24
        (plus maybe 2048(2^11) because centroiding may not work properly?)
        '''
        N_old = len(self.queried)
        if del_flags is not None:
            idx2remove = []
            for i, row in enumerate(self.queried):
                b_flag = list(f"{row['f_objID']:031b}")
                for bin_pos in del_flags:
                    if b_flag[-bin_pos] == '1':
                        idx2remove.append(i)
            self.queried.remove_rows(idx2remove)

            N_fobj = len(self.queried)
            mask_str(N_fobj, N_old, f"f_objID ({del_flags})")

            N_old = N_fobj

        if drop_by_Kron:
            dmag = (self.queried["imag"] - self.queried["iKmag"])
            mask = (dmag > 0.05)
            self.queried = self.queried[~mask]

            N_Kron = len(self.queried)
            mask_str(N_Kron, N_old, "the Kron magnitude criterion")

    def select_filters(
            self,
            filter_names=["g", "r", "i"],
            keep_columns=["_r", "objID", "f_objID", "RAJ2000", "DEJ2000", "x", "y"],
            n_mins=[0, 0, 0]
    ):
        ''' Abridges the columns depending on the specified filters.
        '''
        if not isinstance(filter_names, (list, tuple, np.ndarray)):
            filter_names = [filter_names]

        n_mins = np.atleast_1d(n_mins)
        if n_mins.shape[0] == 1:
            n_mins = np.repeat(n_mins, len(filter_names))
        elif n_mins.shape[0] != len(filter_names):
            raise ValueError("n_mins must be length 1 or same length as "
                             f"filter_names (now it's {len(filter_names)}).")

        selected_columns = keep_columns
        toremove_columns = []

        for filt in filter_names:
            selected_columns.append(f"N{filt}")
            selected_columns.append(f"{filt}mag")
            selected_columns.append(f"{filt}Kmag")
            selected_columns.append(f"{filt}Flags")
            selected_columns.append(f"{filt}PSFf")
            selected_columns.append(f"{filt}magStd")
            selected_columns.append(f"e_{filt}mag")
            selected_columns.append(f"e_{filt}Kmag")
            selected_columns.append(f"o_{filt}mag")
            selected_columns.append(f"b_{filt}mag")
            selected_columns.append(f"B_{filt}mag")

        for c in self.queried.columns:
            if c not in selected_columns:
                toremove_columns.append(c)

        self.queried.remove_columns(toremove_columns)

        N_old = len(self.queried)

        for i, filt in enumerate(filter_names):
            mask = np.array(self.queried[f"o_{filt}mag"]) < n_mins[i]
            self.queried = self.queried[~mask]

        N_new = len(self.queried)
        mask_str(N_new, N_old, f"o_{filter_names}mag >= {n_mins}")

    def check_nearby(self, minsep, maxmag=None, filter_names=["r"]):
        ''' Checkes whether there is any nearby object.
        Notes
        -----
        It checks the ``"_r"`` column of the ``PanSTARRS1`` queried result.
        Therefore, the query center should be the position where you want to
        check for any nearby object.
        Parameters
        ----------
        minsep : float or `~astropy.Quantity`
            The minimum separation to detect nearby object
        maxmag : int or float, optional
            The maximum magnitude value to mask objects. Objects fainter than
            this magnitude (Mean PSF magnitude) will be accepted even though it
            is nearby the search center.
        '''
        if isinstance(minsep, (float, int)):
            warn("minsep is not Quantity. Assuming degree unit.")
            minsep = minsep * u.deg
        elif not isinstance(minsep, u.Quantity):
            raise TypeError("minsep not understood.")

        if not isinstance(filter_names, (list, tuple, np.ndarray)):
            filter_names = [filter_names]

        chktab = self.queried.copy()

        if maxmag is not None:
            for filt in filter_names:
                chktab = chktab[chktab[filt] <= maxmag]
        minsep = minsep.to(chktab["_r"].unit).value
        isnear = (np.array(chktab["_r"]).min() <= minsep)
        return isnear

    def drop_star_groups(self, crit_separation):
        N_old = len(self.queried)
        grouped_rows = group_stars(table=self.queried, crit_separation=crit_separation,
                                   xcol="x", ycol="y", index_only=True)
        self.queried.remove_rows(grouped_rows)
        N_new = len(self.queried)
        mask_str(N_new, N_old,
                 (f"DAOGROUP with {crit_separation:.3f}-pixel critical separation."))


#%%
# Initialize PanSTARRS1 class
q = PanSTARRS1(
    ra=cent_coord.ra.value, 
    dec=cent_coord.dec.value, 
    radius=rad,
    column_filters={"gmag":"13.0..20.0", "e_gmag":"<0.10"}
)

# Query to the website (VizieR)
# This is where the most of the time is spent.
q.query()

#%%
# Only select the stars within 50-pixel bezel in the FOV.
q.select_xyinFOV(header=hdr, bezel=50)

# Remove objects not suitable for differential photometry (see description below)
q.drop_for_diff_phot(drop_by_Kron=True)

# Remove redundant columns, remove objects with too few observations:
q.select_filters(filter_names=['g', 'r', 'i'], n_mins=5)
# You can try a list of ``n_mins``:
# q.select_filters(filter_names=['g', 'r', 'i'], n_mins=[10, 3, 5])

# %%
q_stars_orig = q.queried.copy()
pos_stars_orig = np.array([q_stars_orig["x"], q_stars_orig["y"]]).T
q_stars_orig

#%%
def group_stars(table, crit_separation, xcol="x", ycol="y", index_only=True):
    ''' Group stars using DAOGROUP algorithm and return row indices.
    Parameters
    ----------
    table : astropy.table.Table
        The queried result table.
    crit_separation : float or int
        Distance, in units of pixels, such that any two stars separated by less
        than this distance will be placed in the same group.
    xcol, ycol : str, optional
        The column names for x and y positions. This is necessary since
        `~photutils.DAOGroup accepts a table which has x y positions designated
        as ``"x_0"`` and ``"y_0"``.
    index : bool, optional
        Whether to return only the index of the grouped rows (group information
        is lost) or the full grouped table (after group_by).
    Notes
    -----
    Assuming the psf fwhm to be known, ``crit_separation`` may be set to
    ``k * fwhm``, for some positive real k.
    See Also
    --------
    photutils.DAOStarFinder
    References
    ----------
    [1] Stetson, Astronomical Society of the Pacific, Publications,
        (ISSN 0004-6280), vol. 99, March 1987, p. 191-222.
        Available at: http://adsabs.harvard.edu/abs/1987PASP...99..191S
    Return
    ------
    gtab: Table
        Returned when ``index_only=False``. The table underwent
        ``.group_by("group_id")``.
    grouped_rows: list
        Returned when ``index_only=True``.
        The indices of the rows which are "grouped" stars. You may remove such
        rows using ``table.remove_rows(grouped_rows)``.
    '''
    from photutils.psf.groupstars import DAOGroup
    tab = table.copy()

    tab[xcol].name = "x_0"
    tab[ycol].name = "y_0"
    gtab = DAOGroup(crit_separation=crit_separation)(tab).group_by("group_id")

    if not index_only:
        return gtab
    else:
        gid, gnum = np.unique(gtab["group_id"], return_counts=True)
        gmask = gid[gnum != 1]  # group id with > 1 stars
        grouped_rows = []
        for i, gid in enumerate(gtab["group_id"]):
            if gid in gmask:
                grouped_rows.append(i)
        return grouped_rows

# %%
fwhm = 4
q.drop_star_groups(crit_separation=6*fwhm)
q_stars_dropped = q.queried.copy()  # This will be overridden: see below

#%%
def organize_ps1_and_isnear(
        ps1,
        header=None,
        bezel=0,
        nearby_obj_minsep=0*u.deg,
        group_crit_separation=0,
        select_filter_kw={},
        del_flags=PS1_DR1_DEL_FLAG,
        drop_by_Kron=True,
        calc_JC=True
):
    ''' Organizes the PanSTARRS1 object and check nearby objects.
    Parameters
    ----------
    ps1 : `~PanSTARRS1`
        The `~PanSTARRS1` object.
    header : `astropy.header.Header`, None, optional
        The header to extract WCS related information. If `None` (default), it
        will not drop any stars based on the field of view criterion.
    bezel : int, float, optional
        The bezel used to select stars inside the field of view.
    nearby_obj_minsep : float, `~astropy.Quantity`, optional.
        If there is any object closer than this value, a warning message will
        be printed.
    group_crit_separation : float, optional
        The critical separation parameter used in DAOGROUP algorithm
        (`~photutils.DAOGroup`) to select grouped stars.
    select_filter_kw : dict, optional
        The kwargs for `~PanSTARRS1.select_filter()` method.
    del_flags : list of int, optional
        The flags to be used for dropping objects based on ``"f_objID"`` of
        Pan-STARRS1 query.
    drop_by_Kron : bool, optional
        If `True` (default), drop the galaxies based on the Kron magnitude
        criterion suggested by PS1 (which works good only if i <~ 21):
        https://outerspace.stsci.edu/display/PANSTARRS/How+to+separate+stars+and+galaxies
    calc_JC : bool, optional
        Whether to calculate the Johnson-Cousins B V R_C filter magnitudes by
        the linear relationship given by Table 6 of Tonry J. et al. 2012, ApJ,
        750, 99., using g-r color. The following columns will be added to the
        table ``ps1.queried``:
          * ``"C_gr"``: The ``g-r`` color.
          * ``"dC_gr"``: The total error of ``g-r`` color. Not only the
            error-bar of the mean PSF magnitude (``"e_Xmag`` for filter
            ``X="g"`` and ``X="r"``), but also the intrinsic error-bar of each
            measurements (``"XmagStd" for filter ``X="g"`` and ``X="r"``) are
            considered, i.e., four error-bars are propagated by first order
            approximation (square sum and square rooted).
          * ``"Bmag"``, ``"Vmag"``, ``Rmag``:
            ``B = g + 0.213 + 0.587(g-r) (+- 0.034)``,
            ``V = r + 0.006 + 0.474(g-r) (+- 0.012)``,
            ``R = r - 0.138 -0.131(g-r) (+-0.015)``
          * ``"dBmag"``, ``"dVmag"``, ``"dRmag"``: The total error of above
            magnitudes. The scatter reported by Tonry et al. (e.g., 0.012 mag
            for V) is propagated with the first order error estimated from the
            magnitude calculation formula.
    Returns
    -------
    isnear : bool
        True if there is any nearby object from ``ps1.queried``.
    '''

    _ = ps1.query()

    # Select only those within FOV & bezel.
    # If you wanna keep those outside the edges, just set negative
    # ``bezel``.
    if header is not None:
        ps1.select_xyinFOV(header=header, bezel=bezel)

    # Check whether any object is near our target
    isnear = ps1.check_nearby(minsep=nearby_obj_minsep)
    if isnear:
        warn("There are objects near the target!")

    # Drop objects near to each other
    ps1.drop_star_groups(crit_separation=group_crit_separation)

    # Drop for preparing differential photometry
    ps1.drop_for_diff_phot(del_flags=del_flags, drop_by_Kron=drop_by_Kron)

    # remove columns that are of no interest
    ps1.select_filters(**select_filter_kw)

    ps1.queried["_r"] = ps1.queried["_r"].to(u.arcsec)
    ps1.queried["_r"].format = "%.3f"

    if calc_JC:
        ps1.queried["grcolor"] = ps1.queried["gmag"] - ps1.queried["rmag"]

        # Unfortunately the Std of color is unavailable from the columns
        # provided by PS1 DR1. But error of the mean can be estimated by
        # Gaussian error propagation
        var_g = ps1.queried["e_gmag"]**2
        var_r = ps1.queried["e_rmag"]**2
        dc_gr = np.sqrt(var_g + var_r)
        ps1.queried["e_grcolor"] = dc_gr

        pars = dict(Bmag=[0.213, 0.587, 0.034, "gmag"],
                    Vmag=[0.006, 0.474, 0.012, "rmag"],
                    Rmag=[-0.138, -0.131, 0.015, "rmag"])
        # filter_name = [B_0, B_1, B_sc of Tonry, mag used for conversion]
        for k, p in pars.items():
            ps1.queried[k] = (p[0] + p[1]*ps1.queried["grcolor"] + ps1.queried[p[3]])
            ps1.queried[f"e_{k}"] = np.sqrt((p[1]*ps1.queried["e_grcolor"])**2 + p[2]**2)

    return isnear


# %%
# Query sidereal objects (PS1)
cent_coord = center_radec(ccd_or_header=hdr, center_of_image=True)

rad = fov_radius(header=hdr, unit=u.deg)

# Initialize PanSTARRS1 class
ps1 = PanSTARRS1(
    #ra=cent_coord.ra.value, 
    #dec=cent_coord.dec.value,
    ra=cent_coord.ra, 
    dec=cent_coord.dec, 
    radius=rad,
    column_filters={"gmag":"13.0..20.0", "e_gmag":"<0.10"}
)

_ = organize_ps1_and_isnear(
    ps1=ps1,
    header=hdr,
    bezel=50,
    group_crit_separation=6*fwhm,
    select_filter_kw=dict(filter_names=['g', 'r', 'i'], n_mins=5),
    drop_by_Kron=True,
    calc_JC=True  # This uses g-r color, so both g and r must not be removed.
)
# %%
ps1.queried

#%%
q_stars = ps1.queried.copy()
pos_stars = np.array([q_stars["x"], q_stars["y"]]).T

#%%
ap_stars_orig = CAp(positions=pos_stars_orig, r=15)
ap_stars = CAp(positions=pos_stars, r=20)

fig, axs = plt.subplots(1, 1, figsize=(5, 6), sharex=False, sharey=False, gridspec_kw=None)
yvu.norm_imshow(axs, hdul[0].data, zscale=True)
ap_stars_orig.plot(color='w', lw=2)
ap_stars.plot(color='r', lw=2)

axs.set_title("PS1 Query: Dropping Nearby Stars")
plt.tight_layout()
