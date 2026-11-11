#!/usr/bin/env python3
"""Export reference catalogs for SAMOS preparation; see README and notebook V5.4 cells 85, 106, 114, 116.
Gaia --gaia-input accepts an archived source catalog; --gaia-online uses Gaia DR3 and propagates PM.
SkyMapper and Pan-STARRS catalogs are queried online. No synthetic photometry is generated.
"""
import argparse
from pathlib import Path
from io import StringIO
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import numpy as np
import pandas as pd

SM_COLS = ('r_cntr object_id ngood flags smss_j raj2000 dej2000 e_raj2000 e_dej2000 mean_epoch nimaflags class_star '
 'u_psf e_u_psf u_petro e_u_petro u_apc05 e_u_apc05 v_psf e_v_psf v_petro e_v_petro v_apc05 e_v_apc05 '
 'g_psf e_g_psf g_petro e_g_petro g_apc05 e_g_apc05 r_psf e_r_psf r_petro e_r_petro r_apc05 e_r_apc05 '
 'i_psf e_i_psf i_petro e_i_petro i_apc05 e_i_apc05 z_psf e_z_psf z_petro e_z_petro z_apc05 e_z_apc05 self_dist1 self_id1').split()

def in_square(df, ra_col, dec_col, ra, dec, size_arcmin=3):
    x = (pd.to_numeric(df[ra_col],errors='coerce') - ra) * np.cos(np.deg2rad(dec))*60
    y = (pd.to_numeric(df[dec_col],errors='coerce') - dec)*60
    return df.loc[(x.abs() <= size_arcmin/2)&(y.abs() <= size_arcmin/2)].copy()

def fetch_text(url):
    req=Request(url,headers={'User-Agent':'SAMOS-target-tool/1.1 reference-catalog-export'})
    with urlopen(req,timeout=60) as f:return f.read().decode('utf-8')

def skymapper(ra,dec):
    url='https://skymapper.anu.edu.au/sm-cone/public/query?'+urlencode({'RA':ra,'DEC':dec,'SR':0.05,'VERB':1})
    tables=pd.read_html(StringIO(fetch_text(url)))
    if not tables:raise RuntimeError('SkyMapper returned no HTML table')
    tab=tables[0]
    if isinstance(tab.columns,pd.MultiIndex):tab.columns=tab.columns.get_level_values(-1)
    if len(tab.columns)==len(SM_COLS) and 'raj2000' not in tab.columns:
        tab.columns=SM_COLS # historical notebook response structure
    if not {'raj2000','dej2000'}.issubset(tab.columns):
        raise RuntimeError('SkyMapper column format changed: '+str(list(tab.columns)))
    return in_square(tab,'raj2000','dej2000',ra,dec)

def panstarrs(ra,dec):
    url='https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean.csv?'+urlencode({
        'ra':ra,'dec':dec,'radius':0.036,'columns':'[objID,raMean,decMean,gMeanPSFMag,rMeanPSFMag,iMeanPSFMag,zMeanPSFMag,yMeanPSFMag]'})
    tab=pd.read_csv(StringIO(fetch_text(url)))
    if not {'raMean','decMean'}.issubset(tab.columns):
        raise RuntimeError('Pan-STARRS column format changed: '+str(list(tab.columns)))
    return in_square(tab,'raMean','decMean',ra,dec)

def gaia_online(ra,dec,epoch):
    from astroquery.gaia import Gaia
    from astropy.coordinates import SkyCoord
    from astropy.time import Time
    from astropy import units as u
    Gaia.ROW_LIMIT=-1
    q=Gaia.query_object_async(coordinate=SkyCoord(ra*u.deg,dec*u.deg),width=0.10*u.deg,height=0.10*u.deg)
    tab=q.to_pandas()
    # Retain Gaia DR3 original columns; propagate when proper motion is present.
    tab['ra_now']=tab['ra'].astype(float)
    tab['dec_now']=tab['dec'].astype(float)
    if epoch is not None:
        good=tab[['ra','dec','pmra','pmdec','ref_epoch']].notna().all(axis=1)
        if good.any():
            c=SkyCoord(ra=tab.loc[good,'ra'].to_numpy()*u.deg,
                dec=tab.loc[good,'dec'].to_numpy()*u.deg,
                pm_ra_cosdec=tab.loc[good,'pmra'].to_numpy()*u.mas/u.yr,
                pm_dec=tab.loc[good,'pmdec'].to_numpy()*u.mas/u.yr,
                obstime=Time(tab.loc[good,'ref_epoch'].to_numpy(),format='jyear'))
            p=c.apply_space_motion(new_obstime=Time(epoch,format='jyear'))
            tab.loc[good,'ra_now']=p.ra.deg
            tab.loc[good,'dec_now']=p.dec.deg
    return in_square(tab,'ra_now','dec_now',ra,dec)

def gaia_local(path,ra,dec):
    tab=pd.read_csv(path,dtype={'target_id':'string'})
    required={'target_id','RA','DEC','Gmag'}
    if not required.issubset(tab.columns):raise ValueError('Gaia input missing columns: '+str(required-set(tab.columns)))
    tab=tab.rename(columns={'target_id':'SOURCE_ID','RA':'ra_now','DEC':'dec_now','Gmag':'phot_g_mean_mag'})
    return in_square(tab,'ra_now','dec_now',ra,dec)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--name',required=True)
    ap.add_argument('--center',nargs=2,type=float,required=True,metavar=('RA_DEG','DEC_DEG'))
    ap.add_argument('--gaia-input',help='Local Gaia catalog; coordinates are NOT epoch-propagated')
    ap.add_argument('--gaia-online',action='store_true')
    ap.add_argument('--epoch',type=float,help='Julian observation year, for online Gaia PM propagation')
    ap.add_argument('--survey',choices=['skymapper','panstarrs'],default='skymapper')
    ap.add_argument('--outdir',default='reference_catalogs')
    a=ap.parse_args();ra,dec=a.center
    if a.gaia_online==bool(a.gaia_input):ap.error('Supply exactly one of --gaia-input or --gaia-online')
    out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True)
    stem=f'{a.name}_RADEC={ra:.7f}{dec:+.7f}'
    gaia=gaia_online(ra,dec,a.epoch) if a.gaia_online else gaia_local(a.gaia_input,ra,dec)
    gcols=['SOURCE_ID','ra_now','dec_now','phot_g_mean_mag']
    if 'source_id' in gaia.columns and 'SOURCE_ID' not in gaia.columns:gaia=gaia.rename(columns={'source_id':'SOURCE_ID'})
    gaia=gaia[gcols].dropna(subset=['ra_now','dec_now']).sort_values('phot_g_mean_mag')
    gp=out/f'{stem}_GAIA_in_field.csv';gaia.to_csv(gp,index=False)
    print(f'GAIA {len(gaia)} sources -> {gp}')
    if not a.gaia_online:print('NOTE: local Gaia coordinates are input-catalog epoch; no PM correction performed')
    try:
        phot=skymapper(ra,dec) if a.survey=='skymapper' else panstarrs(ra,dec)
        pp=out/f'{stem}_{a.survey.upper()}_phot_in_field.csv'
        phot.to_csv(pp,index=False)
        print(f'{a.survey.upper()} {len(phot)} sources -> {pp}')
    except Exception as e:
        print('PHOTOMETRIC CATALOG NOT WRITTEN: '+str(e))
        raise
if __name__=='__main__':main()
