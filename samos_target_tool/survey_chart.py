"""Astronomical-image finding charts for SAMOS, using the original FITS celestial WCS.

An actual FITS file is mandatory. No simulated image is substituted.
"""
from pathlib import Path
import numpy as np


def render_chart(image_path, output_path, slits, center_ra, center_dec,
                 *, width_arcsec=0.54, length_arcsec=1.62,
                 field_arcsec=180.0, show_traces=False):
    from astropy.io import fits
    from astropy.wcs import WCS
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    from astropy.visualization import ImageNormalize, AsinhStretch, PercentileInterval

    with fits.open(image_path, memmap=False) as hdul:
        hdu = next((h for h in hdul if h.data is not None and h.data.ndim == 2), None)
        if hdu is None:
            raise ValueError('A real 2D FITS image is required')
        wcs = WCS(hdu.header).celestial
        if not wcs.has_celestial:
            raise ValueError('FITS image lacks a celestial WCS')
        data = np.array(hdu.data, dtype=float)
    if not np.isfinite(data).any():
        raise ValueError('Image has no finite pixels')
    x0, y0 = wcs.world_to_pixel_values(center_ra, center_dec)
    if not (0 <= x0 < data.shape[1] and 0 <= y0 < data.shape[0]):
        raise ValueError('Pointing star is outside image')
    fig = plt.figure(figsize=(9, 9))
    ax = fig.add_subplot(111, projection=wcs)
    norm = ImageNormalize(data, interval=PercentileInterval(99.5), stretch=AsinhStretch())
    ax.imshow(data, origin='lower', cmap='gray', norm=norm, interpolation='nearest')
    half = field_arcsec / 2
    def rect_vertices(ra, dec, half_w, half_l):
        # Tangent-plane offsets; longitude corrected by cos(dec).
        return [[ra + dx / (3600*np.cos(np.deg2rad(dec))), dec + dy/3600]
                for dx,dy in [(-half_w,-half_l),(half_w,-half_l),(half_w,half_l),(-half_w,half_l)]]
    def add_rectangle(ra, dec, hw, hl, **kwargs):
        corners = np.asarray(rect_vertices(ra,dec,hw,hl))
        xp, yp = wcs.world_to_pixel_values(corners[:,0],corners[:,1])
        ax.add_patch(Polygon(np.column_stack([xp,yp]),closed=True,fill=False,**kwargs))
    add_rectangle(center_ra,center_dec,half,half,color='yellow',linewidth=1.4)
    colors={'SCIENCE':'lime','FILLER':'orange','SKY':'cyan'}
    for row in slits.itertuples():
        c=colors.get(str(row.SLITTYPE),'lime')
        add_rectangle(float(row.RA),float(row.DEC),length_arcsec/2,width_arcsec/2,
                      edgecolor=c,linewidth=1.1)
        if not show_traces:
            xp,yp=wcs.world_to_pixel_values(float(row.RA),float(row.DEC))
            ax.annotate(str(row.SLITID).replace('SLIT',''),(xp,yp),xytext=(3,3),
                        textcoords='offset points',fontsize=5,color=c)
    # Do not draw invented detector trace stripes on sky-coordinate imagery.
    if show_traces:
        ax.set_title('SAMOS slit positions — trace geometry NOT calibrated')
    else:
        ax.set_title('SAMOS target/slit positions on survey FITS')
    ax.scatter([x0],[y0],marker='+',s=130,c='yellow',linewidths=1.6)
    # Constant center + 3.2 arcmin preview is appropriate for finder charts.
    ax.coords[0].set_axislabel('Right ascension (ICRS)')
    ax.coords[1].set_axislabel('Declination (ICRS)')
    ax.grid(color='white',alpha=0.25,linestyle=':')
    fig.text(.5,.015,'REAL FITS background; slits from astronomical RA/Dec; SKY apertures unverified for nebulosity',
             ha='center',fontsize=8)
    p=Path(output_path);p.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(p,dpi=180,bbox_inches='tight');plt.close(fig)
    return p
