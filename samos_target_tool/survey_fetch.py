"""Retrieve ORIGINAL sky-survey FITS cutouts for SAMOS finding charts.

SkyMapper DR4 SIAP: https://skymapper.anu.edu.au/how-to-access/
Pan-STARRS PS1: https://outerspace.stsci.edu/spaces/PANSTARRS/
No plotting or fabricated sky backgrounds in this module.
"""
import argparse
import csv
import io
from pathlib import Path
from urllib.parse import urljoin
import requests

SKYMAPPER = 'https://api.skymapper.nci.org.au/public/siap/dr4/query'
PS1 = 'https://ps1images.stsci.edu/cgi-bin/'

def retrieve(url, *, params=None, timeout=45):
    r=requests.get(url,params=params,timeout=timeout,headers={'User-Agent':'SAMOS-Target-Tool/1.0 (astronomy research)'})
    r.raise_for_status()
    return r

def _check_fits(content):
    # FITS primary HDU starts with SIMPLE, or extension-based file can begin XTENSION.
    if not (content.startswith(b'SIMPLE  =') or content.startswith(b'XTENSION=')):
        raise ValueError('The service did not return FITS bytes (HTML/error page?).')

def _save_fits(content,outfile,ra,dec):
    _check_fits(content)
    try:
        from astropy.io import fits
        from astropy.wcs import WCS
        import numpy as np
    except ImportError as e:
        raise RuntimeError('Install astropy to validate FITS WCS before saving') from e
    import io as _io
    with fits.open(_io.BytesIO(content),checksum=False) as hdul:
        h=next((h for h in hdul if getattr(h,'data',None) is not None and getattr(h.data,'ndim',0)==2),None)
        if h is None: raise ValueError('FITS has no 2D image HDU')
        w=WCS(h.header).celestial
        if not w.has_celestial: raise ValueError('FITS lacks valid celestial WCS')
        ny,nx=h.data.shape
        # Requested center must actually be included in cutout.
        x,y=w.world_to_pixel_values(ra,dec)
        if not (np.isfinite(x) and np.isfinite(y) and 0<=x<nx and 0<=y<ny):
            raise ValueError('Requested pointing star is outside returned FITS image')
    path=Path(outfile); path.parent.mkdir(parents=True,exist_ok=True)
    path.write_bytes(content)
    return str(path)

def skymapper(ra,dec,outfile,band='i',width_arcmin=3.5):
    """SkyMapper DR4 SIAP CSV metadata -> exact get_fits/get_image URL."""
    params={'POS':f'{ra:.10f},{dec:.10f}', 'SIZE':f'{width_arcmin/60:.9f}',
            'BAND':band, 'FORMAT':'image/fits','INTERSECT':'COVERS',
            'RESPONSEFORMAT':'CSV','VERB':3}
    result=retrieve(SKYMAPPER,params=params)
    rows=list(csv.DictReader(io.StringIO(result.text.lstrip('\ufeff'))))
    if not rows: raise RuntimeError('SkyMapper DR4 has no fully covering cutout at this center/band')
    def url_for(row):
        fields={k.strip().lower():v for k,v in row.items() if k}
        for key in ('get_fits','get_image','access_url','accessurl'):
            value=fields.get(key,'').strip()
            if value and ('fits' in value.lower() or key=='get_fits'):
                return urljoin(SKYMAPPER,value)
        return None
    # Preserve source ordering; first fully covering candidate selected. No undocumented image-quality ranking.
    urls=[url_for(row) for row in rows]
    urls=[url for url in urls if url]
    if not urls: raise RuntimeError('SkyMapper metadata did not provide a FITS cutout URL')
    errors=[]
    for url in urls:
        try: return _save_fits(retrieve(url).content,outfile,ra,dec)
        except (requests.RequestException,ValueError,RuntimeError) as e: errors.append(str(e))
    raise RuntimeError('No valid SkyMapper FITS candidate: '+ '; '.join(errors[:3]))

def panstarrs(ra,dec,outfile,band='i',width_arcmin=3.5):
    """MAST PS1 filenames list followed by FITS cutout (0.25 arcsec PS1 pixels)."""
    result=retrieve(PS1+'ps1filenames.py',params={'ra':ra,'dec':dec,'filters':band})
    # ASCII whitespace-delimited service response; use header field names.
    lines=[line.strip() for line in result.text.splitlines() if line.strip() and not line.startswith('#')]
    if len(lines)<2: raise RuntimeError('Pan-STARRS: no matching image found')
    header=lines[0].split(); vals=lines[1].split()
    if 'filename' not in header: raise RuntimeError('Pan-STARRS service response missing filename column')
    filename=vals[header.index('filename')]
    size=round(width_arcmin*60/0.25)
    r=retrieve(PS1+'fitscut.cgi',params={'ra':ra,'dec':dec,'size':size,'format':'fits','red':filename})
    return _save_fits(r.content,outfile,ra,dec)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ra',type=float,required=True);p.add_argument('--dec',type=float,required=True)
    p.add_argument('--survey',choices=['auto','skymapper','panstarrs'],default='auto')
    p.add_argument('--band',default='i');p.add_argument('--width-arcmin',type=float,default=3.5)
    p.add_argument('--out',required=True)
    a=p.parse_args()
    if not (0<a.width_arcmin<10):p.error('width must be >0 and <10 arcmin')
    if a.band not in ('u','v','g','r','i','z','y'):p.error('unsupported filter')
    choices=['skymapper','panstarrs'] if a.survey=='auto' else [a.survey]
    errors=[]
    for choice in choices:
        try:
            path=(skymapper if choice=='skymapper' else panstarrs)(a.ra,a.dec,a.out,a.band,a.width_arcmin)
            print(f'SURVEY={choice} FITS={path}');return
        except Exception as e:
            errors.append(f'{choice}: {type(e).__name__}: {e}')
    p.exit(2,'Survey FITS retrieval FAILED (no synthetic substitution):\n'+'\n'.join(errors)+'\n')
if __name__=='__main__':main()
