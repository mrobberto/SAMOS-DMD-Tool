"""Schematic (NOT detector-calibrated) projected SAMI spectral stripes.

Operates on finalized ALL-mask CSVs and its derived EVEN/ODD subsets.
No assignment, aperture coordinates or slit identifiers are modified.
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd


def make_stripes(mask_csv, pointing_ra, pointing_dec, outdir=None,
                 slit_length_arcsec=1.62, spectral_length_pix=2883,
                 nominal_width_pix=4000, nominal_height_pix=1300):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    mask_csv = Path(mask_csv)
    out = Path(outdir) if outdir else mask_csv.parent
    out.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(mask_csv, dtype={'SLITID':'string', 'target_id':'string'})
    req = {'RA','DEC','SLITID','TRACESET','SLITTYPE'}
    missing = req - set(d.columns)
    if missing: raise ValueError(f'{mask_csv}: missing columns {sorted(missing)}')
    if len(d) == 0: raise ValueError('Empty mask')
    if not set(d.TRACESET).issubset({'EVEN','ODD'}): raise ValueError('Invalid TRACESET')
    if spectral_length_pix<=0 or nominal_width_pix<=0 or nominal_height_pix<=0: raise ValueError('Invalid dimensions')
    # RA is the nominal cross-dispersion coordinate in the SAMOS mask convention.
    cross_arcsec=(d.RA.to_numpy(dtype=float)-pointing_ra)*3600*np.cos(np.deg2rad(pointing_dec))
    disp_arcsec=(d.DEC.to_numpy(dtype=float)-pointing_dec)*3600
    y=nominal_height_pix/2 + cross_arcsec*(nominal_height_pix/180)
    x=160 + (disp_arcsec+90)*(750/180)
    h=slit_length_arcsec*(nominal_height_pix/180)
    if np.any((y-h/2<0)|(y+h/2>nominal_height_pix)):
        raise ValueError('Some slit footprints extend outside the nominal cross-dispersion field')
    collisions=[]
    for i in range(len(d)):
        for j in range(i+1,len(d)):
            if abs(y[i]-y[j])<h and min(x[i]+spectral_length_pix,x[j]+spectral_length_pix)>max(x[i],x[j]):
                collisions.append((str(d.iloc[i].SLITID),str(d.iloc[j].SLITID)))
    saved=[]
    for suffix,indices in [('ALL',np.arange(len(d))),('EVEN',np.flatnonzero(d.TRACESET.to_numpy()=='EVEN')),('ODD',np.flatnonzero(d.TRACESET.to_numpy()=='ODD'))]:
        fig,ax=plt.subplots(figsize=(12,4.4),dpi=170)
        ax.set_facecolor('#f8f8f8');ax.set_xlim(0,nominal_width_pix);ax.set_ylim(0,nominal_height_pix)
        for i in indices:
            r=d.iloc[i];col='#1f77b4' if r.SLITTYPE=='SKY' else ('#a4428f' if r.SLITTYPE=='REPEAT' else ('#da8c18' if r.SLITTYPE=='FILLER' else '#555555'))
            ax.add_patch(Rectangle((float(x[i]),float(y[i]-h/2)),spectral_length_pix,h,facecolor=col,edgecolor='black',lw=.25,alpha=.68))
            ax.text(max(0,float(x[i]-18)),float(y[i]),str(r.SLITID)[-3:],fontsize=5,ha='right',va='center',color='black')
        ax.set_xlabel('Schematic dispersion coordinate (not calibrated SAMI pixels)')
        ax.set_ylabel('Cross-dispersion (scaled from projected RA)')
        ax.set_title(f'{mask_csv.stem.split("_RADEC=")[0]} {suffix} · {len(indices)} apertures — nominal {nominal_width_pix} × {nominal_height_pix}')
        ax.text(.01,-.25,'Schematic trace footprint; detector origin, tilt and distortion not calibrated',transform=ax.transAxes,fontsize=8)
        fig.subplots_adjust(left=.08,right=.99,top=.89,bottom=.29)
        # Filename uses actual mask name, retaining pointing-coordinate tag.
        path=out/(mask_csv.stem+f'_{suffix}_spectral_stripes_schematic.png')
        fig.savefig(path,dpi=170);plt.close(fig);saved.append(path)
    return saved, collisions


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mask-csv',required=True)
    p.add_argument('--pointing-star',type=float,nargs=2,required=True)
    p.add_argument('--outdir')
    p.add_argument('--slit-length-arcsec',type=float,default=1.62)
    p.add_argument('--spectral-length-pix',type=int,default=2883)
    args=p.parse_args()
    paths,collisions=make_stripes(args.mask_csv,*args.pointing_star,outdir=args.outdir,slit_length_arcsec=args.slit_length_arcsec,spectral_length_pix=args.spectral_length_pix)
    print('SCHEMATIC only; not a calibrated collision guarantee')
    print('Rectangular footprint overlaps:',len(collisions))
    for p in paths:print(p)
if __name__=='__main__':main()
