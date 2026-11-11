"""Create per-T EVEN/ODD sky-region products and catalog-based inspection figures.

Not a detector/DMD command generator. Requires numpy,pandas,matplotlib.
Optional FITS backdrop requires astropy and a user-supplied WCS FITS file.
"""
import argparse
import shutil
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


def region_file(path, df, w_arcsec, l_arcsec):
    with open(path,'w') as f:
        f.write('# Region file format: DS9 version 4.1\n')
        f.write('global color=green width=1 font="helvetica 10 normal"\n')
        f.write('icrs\n')
        for r in df.itertuples():
            color={'SKY':'cyan','FILLER':'yellow','REPEAT':'magenta','SCIENCE':'green'}.get(r.SLITTYPE,'green')
            f.write(f'box({r.RA:.9f},{r.DEC:.9f},{l_arcsec:.3f}\",{w_arcsec:.3f}\",0) # color={color} text={{{r.SLITID}}}\n')


def offsets(df,ra,dec):
    return ((df.RA.to_numpy(dtype=float)-ra)*3600*np.cos(np.deg2rad(dec)),(df.DEC.to_numpy(dtype=float)-dec)*3600)


def sky_points(n, selected, reference, ra, dec, exclusion, minsep, half, seed):
    if n == 0:return []
    rng=np.random.default_rng(seed)
    sx,sy=offsets(selected,ra,dec)
    rx,ry=offsets(reference,ra,dec)
    result=[]
    for x,y in rng.uniform(-half+6,half-6,size=(15000,2)):
        if len(rx) and np.min((rx-x)**2+(ry-y)**2) < exclusion**2:continue
        if len(sx) and np.min((sx-x)**2+(sy-y)**2) < minsep**2:continue
        if result and min((x-p[0])**2+(y-p[1])**2 for p in result)<(minsep*2)**2:continue
        # current occupancy model rejects spectra in nearly same cross-dispersion stripe
        if len(sx) and np.min(np.abs(sx-x)) < minsep:continue
        if result and min(abs(x-p[0]) for p in result)<minsep:continue
        result.append((float(x),float(y)))
        if len(result)==n:break
    return result


def choose_gaia_fillers(science, reference, already_used, ra, dec, n_max, gmin, gmax, sep_ra, spread, seed):
    """Add conservative calibration-star candidates, without displacing science.

    Projected RA separation is checked against every selected slit, and sky
    (cross-dispersion) separation is checked against all slit types.  This is a
    geometric screening step, not a calibrated spectral collision validation.
    """
    if n_max <= 0:
        return []
    required = {'target_id','RA','DEC','Gmag'}
    if not required.issubset(reference.columns):
        raise ValueError('Filler catalog must include '+', '.join(sorted(required)))
    pool = reference.copy()
    pool['target_id'] = pool['target_id'].astype(str)
    pool['Gmag'] = pd.to_numeric(pool.Gmag, errors='coerce')
    pool = pool[pool.Gmag.between(gmin,gmax) & ~pool.target_id.isin(already_used)].copy()
    if pool.empty: return []
    fx,fy=offsets(pool,ra,dec)
    pool=pool[(np.abs(fx)<86)&(np.abs(fy)<86)].copy()
    if pool.empty: return []
    # Spread evenly in projected cross-dispersion (RA) while preferring
    # moderately bright stellar calibrators; deterministic seed breaks ties.
    rng=np.random.default_rng(seed)
    pool['x_offset'],pool['y_offset']=offsets(pool,ra,dec)
    pool['random_tie']=rng.random(len(pool))
    pool=pool.sort_values(['Gmag','random_tie']).reset_index(drop=True)
    selected=[]
    sx,sy=offsets(science,ra,dec)
    xlist=list(sx); ylist=list(sy)
    while len(selected)<n_max and not pool.empty:
        best=None; best_score=-1e50
        for i,r in pool.iterrows():
            x=float(r.x_offset); y=float(r.y_offset)
            # Minimum projected RA distance from ALL selected apertures.
            if xlist and min(abs(x-xx) for xx in xlist)<sep_ra:continue
            # A separate, conservative same-trace exclusion check.
            if xlist and min(abs(x-xx) for xx in xlist)<spread:continue
            # Prefer good spread in the cross-dispersion direction.
            clearance=min((abs(x-xx) for xx in xlist),default=90.0)
            score=min(clearance,30.0) - 0.35*(float(r.Gmag)-gmin) + 1e-4*float(r.random_tie)
            if score>best_score: best=i; best_score=score
        if best is None:break
        r=pool.loc[best]
        selected.append({'target_id':str(r.target_id),'RA':float(r.RA),'DEC':float(r.DEC),
                         'Gmag':float(r.Gmag),'priority':98,'SLITTYPE':'FILLER'})
        xlist.append(float(r.x_offset)); ylist.append(float(r.y_offset))
        pool=pool.drop(index=best)
    return selected



def choose_repeat_targets(current, assignments, mask_index, ra, dec, n_max, sep_ra, trace_sep):
    """Use earlier science targets only after the new science assignment is frozen.

    Candidates are sorted faintest first; accepted targets cannot move or displace
    existing science or reserved sky apertures. The simplified RA trace-clearance
    test is the same conservative condition used for Gaia fillers.
    """
    if mask_index == 0 or n_max == 0:
        return []
    pool = assignments[(assignments.assigned_mask >= 0) &
                       (assignments.assigned_mask < mask_index)].copy()
    if pool.empty:
        return []
    pool['Gmag'] = pd.to_numeric(pool.Gmag, errors='coerce')
    pool = pool.dropna(subset=['Gmag']).sort_values(['Gmag','RA','DEC'],
                                                     ascending=[False,True,True],kind='stable')
    used_ids = set(current.target_id.astype(str))
    xlist = list(offsets(current,ra,dec)[0])
    chosen = []
    for r in pool.itertuples():
        if str(r.target_id) in used_ids:
            continue
        x,y=offsets(pd.DataFrame([{'RA':r.RA,'DEC':r.DEC}]),ra,dec)
        x=float(x[0]);y=float(y[0])
        if abs(x) > 90-0.81 or abs(y) > 90-0.27:
            continue
        if xlist and min(abs(x-other) for other in xlist) < max(sep_ra,trace_sep):
            continue
        original=f'T{int(r.assigned_mask):02d}'
        chosen.append({'target_id':str(r.target_id), 'RA':float(r.RA), 'DEC':float(r.DEC),
                       'Gmag':float(r.Gmag), 'priority':int(r.priority),
                       'assigned_mask':mask_index, 'SLITTYPE':'REPEAT',
                       'FIRST_TMASK':original})
        xlist.append(x)
        used_ids.add(str(r.target_id))
        if len(chosen) >= n_max:
            break
    return chosen

def plot_one(out, selected, reference, ra,dec,point_ra,point_dec, width_arcsec, length_arcsec, traces=False, fits_image=None):
    fig,ax=plt.subplots(figsize=(8,8))
    if fits_image:
        try:
            from astropy.io import fits
            from astropy.wcs import WCS
            from astropy.coordinates import SkyCoord
            import astropy.units as u
            with fits.open(fits_image) as hdul:
                h=next(h for h in hdul if h.data is not None and h.data.ndim==2)
                arr=h.data.astype(float); w=WCS(h.header).celestial
                yy,xx=np.indices(arr.shape)
                # full image can be large; inspect with input coordinate grid limited for simplicity
                xs,ys=w.pixel_to_world_values(xx[::max(1,arr.shape[0]//700),::max(1,arr.shape[1]//700)],yy[::max(1,arr.shape[0]//700),::max(1,arr.shape[1]//700)])
                sl=arr[::max(1,arr.shape[0]//700),::max(1,arr.shape[1]//700)]
                ox=-(xs-ra)*3600*np.cos(np.deg2rad(dec)); oy=(ys-dec)*3600
                lo,hi=np.nanpercentile(sl,[10,99.5]);ax.pcolormesh(ox,oy,sl,cmap='gray',shading='auto',vmin=lo,vmax=hi,rasterized=True)
        except Exception as e:
            raise RuntimeError(f'FITS background cannot be loaded: {e}')
    else:
        ax.set_facecolor('#f3f3f3')
        ax.text(0,-83,'CATALOG POSITIONS ONLY — NO SURVEY IMAGE',ha='center',fontsize=8)
    rx,ry=offsets(reference,ra,dec);ax.scatter(-rx,ry,s=9,c='0.6',alpha=.45,label='Reference catalog')
    sx,sy=offsets(selected,ra,dec)
    for r,x,y in zip(selected.itertuples(),sx,sy):
        x=-x
        col={'SKY':'tab:cyan','FILLER':'tab:orange','REPEAT':'tab:purple','SCIENCE':'tab:red'}.get(r.SLITTYPE,'tab:red')
        ax.add_patch(Rectangle((x-length_arcsec/2,y-width_arcsec/2),length_arcsec,width_arcsec,fill=False,edgecolor=col,lw=1.1))
        if traces:ax.add_patch(Rectangle((x-length_arcsec/2-0.4,y-90),length_arcsec+0.8,180,fill=False,edgecolor=col,alpha=.25,lw=.7,linestyle='--'))
        if not traces:ax.text(x+2,y+2,r.SLITID[-3:],fontsize=5)
    px=-(point_ra-ra)*3600*np.cos(np.deg2rad(dec));py=(point_dec-dec)*3600
    ax.scatter([px],[py],marker='*',s=90,c='gold',edgecolors='black',zorder=10,label='Reference star')
    ax.add_patch(Rectangle((-90,-90),180,180,fill=False,lw=1.4,edgecolor='black'))
    ax.set_xlim(-94,94);ax.set_ylim(-94,94);ax.set_aspect('equal');ax.set_xlabel('Offset (arcsec), east to left');ax.set_ylabel('Offset north (arcsec)')
    ax.set_title(out.stem+' — '+('approximate trace clearance' if traces else 'slit positions'))
    if traces:ax.text(0,-83,'SCHEMATIC DISPERSION DIRECTION — NOT CALIBRATED SAMI FOOTPRINT',ha='center',fontsize=7)
    fig.tight_layout();fig.savefig(out,dpi=160);plt.close(fig)


def run(args):
    assignments=pd.read_csv(args.assignments,dtype={'target_id':str})
    reference=pd.read_csv(args.reference,dtype={'target_id':str})
    if not {'RA','DEC','assigned_mask'}.issubset(assignments):raise ValueError('Assignment table lacks required columns')
    out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True)
    point_ra,point_dec=args.pointing_star
    # SAMOS reference pointing star DEFINES the field center.
    ra,dec=point_ra,point_dec
    if args.field_center is not None and (abs(args.field_center[0]-ra)>1e-9 or abs(args.field_center[1]-dec)>1e-9):
        raise ValueError("Field center must coincide with pointing-star coordinates; omit --field-center or supply the same coordinates")
    if args.max_fillers < 0 or args.max_repeats < 0 or args.n_sky < 0: raise ValueError('Counts must be nonnegative')
    if args.filler_min_ra_sep < 0: raise ValueError('RA separation must be nonnegative')
    if args.filler_g_min > args.filler_g_max: raise ValueError('Invalid filler magnitude range')
    if 'target_id' not in reference.columns or 'target_id' not in assignments.columns: raise ValueError('Missing target_id')
    science_ids=set(assignments.target_id.astype(str))
    all_rows=[];summ=[]
    for mask in sorted(assignments.loc[assignments.assigned_mask>=0,'assigned_mask'].astype(int).unique()):
        t=f'T{mask:02d}'
        sub=assignments[assignments.assigned_mask==mask].copy().reset_index(drop=True)
        sub['SLITTYPE']='SCIENCE'
        # Science assignments were computed upstream. Re-center onto reference star;
        # flag any science targets outside the real 3x3 arcmin field.
        x,y=offsets(sub,ra,dec)
        if np.any((np.abs(x)>90-args.slit_length_arcsec/2)|(np.abs(y)>90-args.slit_width_arcsec/2)):
            bad=sub.loc[(np.abs(x)>90-args.slit_length_arcsec/2)|(np.abs(y)>90-args.slit_width_arcsec/2),'target_id'].tolist()
            raise ValueError(f'{t}: {len(bad)} assigned science targets outside star-centered field. Re-run optimization at the new pointing before finalizing. Sample: {bad[:5]}')
        sub['FIRST_TMASK']=t
        # Reserve the required blank sky apertures BEFORE secondary targets,
        # so repeats and Gaia fillers cannot displace sky calibration slits.
        sk=sky_points(args.n_sky,sub,reference,ra,dec,args.sky_exclusion_arcsec,
                      args.sky_trace_separation_arcsec,90,args.seed+mask)
        if len(sk) != args.n_sky:
            raise ValueError(f'{t}: could reserve only {len(sk)}/{args.n_sky} sky apertures')
        if sk:
            sky_rows = pd.DataFrame([{'target_id':f'SKY_{t}_{j:02d}',
                'RA':ra+xx/3600/np.cos(np.deg2rad(dec)),
                'DEC':dec+yy/3600,'priority':99,'assigned_mask':mask,
                'SLITTYPE':'SKY','FIRST_TMASK':''} for j,(xx,yy) in enumerate(sk)])
            sub=pd.concat([sub,sky_rows],ignore_index=True)
        repeats=choose_repeat_targets(sub,assignments,mask,ra,dec,
                  args.max_repeats if args.repeat_targets else 0,
                  args.filler_min_ra_sep,args.filler_trace_sep)
        if repeats:
            sub=pd.concat([sub,pd.DataFrame(repeats)],ignore_index=True)
        fillers=choose_gaia_fillers(sub,reference,science_ids,ra,dec,
                args.max_fillers if args.gaia_fillers else 0,
                args.filler_g_min,args.filler_g_max,args.filler_min_ra_sep,
                args.filler_trace_sep,args.seed+mask)
        if fillers:
            sub=pd.concat([sub,pd.DataFrame([dict(fill,assigned_mask=mask,FIRST_TMASK='') for fill in fillers])],ignore_index=True)
        sub=sub.sort_values(['RA','DEC'],kind='stable').reset_index(drop=True)
        # The ALL file is the authoritative science mask. Zero-based interleaved
        # numbering is derived ONLY after this full-mask collision/filler/sky pass.
        sub['SET_RANK']=np.arange(len(sub))//2
        sub['SLITID']=[f'SLIT{k:03d}' for k in range(len(sub))]
        sub['TRACESET']=np.where(np.arange(len(sub))%2==0,'EVEN','ODD')
        sub['TMASK']=t;sub['MASK_VER']='v000';sub['APERTURE_KEY']=sub['TMASK']+'_'+sub['SLITID']
        stem=f'{args.name}-{t}-all_RADEC={point_ra:.7f}{point_dec:+.7f}'
        # Temporary CSV used only as an input to stripe diagnostics.
        sub.to_csv(out/(stem+'.csv'),index=False)
        region_file(out/(stem+'.reg'),sub,args.slit_width_arcsec,args.slit_length_arcsec)
        if args.fits_image:
            from .survey_chart import render_chart
            render_chart(args.fits_image,out/(stem+'_survey_slits.png'),sub,ra,dec,
                         width_arcsec=args.slit_width_arcsec,length_arcsec=args.slit_length_arcsec)
        else:
            plot_one(out/(stem+'_CATALOG_ONLY.png'),sub,reference,ra,dec,point_ra,point_dec,
                     args.slit_width_arcsec,args.slit_length_arcsec)
        if args.spectral_stripes:
            from .spectral_stripes import make_stripes
            stripe_paths, collisions = make_stripes(out/(stem+'.csv'),ra,dec,outdir=out,
                slit_length_arcsec=args.slit_length_arcsec,
                spectral_length_pix=args.spectral_length_pix)
            if collisions:
                print(f'WARNING: {t} nominal spectral stripes have {len(collisions)} rectangular overlaps; these are schematic only')
        (out/(stem+'.csv')).unlink()
        summ.append({'T':t,'set':'ALL','science':int((sub.SLITTYPE=='SCIENCE').sum()),
                     'repeats':len(repeats),'fillers':len(fillers),'sky':len(sk),'filename':stem+'.reg'})
        all_rows.append(sub)
        for parity in ('EVEN','ODD'):
            part=sub[sub.TRACESET==parity].copy()
            p_stem=f'{args.name}-{t}-{parity}_RADEC={point_ra:.7f}{point_dec:+.7f}'
            region_file(out/(p_stem+'.reg'),part,args.slit_width_arcsec,args.slit_length_arcsec)
            summ.append({'T':t,'set':parity,'science':int((part.SLITTYPE=='SCIENCE').sum()),
                         'repeats':int((part.SLITTYPE=='REPEAT').sum()),
                         'fillers':int((part.SLITTYPE=='FILLER').sum()),
                         'sky':int((part.SLITTYPE=='SKY').sum()),'filename':p_stem+'.reg'})
    if not all_rows:raise ValueError('No assigned masks in input catalog')
    full = pd.concat(all_rows,ignore_index=True)
    # Cross-configuration identity: the same Gaia source may have a different
    # SLITID and parity in each T pattern. Never use SLITID as a source key.
    full['target_id'] = full['target_id'].astype(str)
    full['OBS_SEQUENCE'] = full.groupby('target_id',sort=False).cumcount()+1
    full['N_CONFIGS'] = full.groupby('target_id')['TMASK'].transform('nunique')
    presence = (full.groupby('target_id',sort=False)
                .agg(OBSERVED_IN=('TMASK',lambda x: ';'.join(dict.fromkeys(x))),
                     APERTURE_HISTORY=('APERTURE_KEY',lambda x: ';'.join(x)))
                .reset_index())
    full=full.merge(presence,on='target_id',how='left',validate='many_to_one')
    full.to_csv(out/'slit_manifest.csv',index=False)
    source_rows=[]
    for target_id,group in full[full.SLITTYPE!='SKY'].groupby('target_id',sort=False):
        group=group.sort_values(['TMASK','SLITID'],kind='stable')
        first=group.iloc[0]
        source_rows.append({'target_id':target_id,'RA':first.RA,'DEC':first.DEC,
            'Gmag':first.get('Gmag',np.nan),'priority':first.get('priority',np.nan),
            'FIRST_TMASK':first.FIRST_TMASK,'N_CONFIGS':group.TMASK.nunique(),
            'OBSERVED_IN':';'.join(group.TMASK),
            'APERTURE_HISTORY':';'.join(group.APERTURE_KEY),
            'SLIT_TYPES':';'.join(group.SLITTYPE)})
    pd.DataFrame(source_rows).to_csv(out/'source_observation_index.csv',index=False)
    pd.DataFrame(summ).to_csv(out/'mask_summary.csv',index=False)
    print(pd.DataFrame(summ).to_string(index=False))
    return pd.DataFrame(summ)


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--assignments',required=True);p.add_argument('--reference',required=True);p.add_argument('--outdir',required=True)
    p.add_argument('--name',required=True,help='Observer-supplied field name; never hardcoded');p.add_argument('--field-center',type=float,nargs=2,default=None)
    p.add_argument('--pointing-star',type=float,nargs=2,required=True)
    p.add_argument('--gaia-fillers',action=argparse.BooleanOptionalAction,default=True,
        help='Add optional Gaia stellar calibration targets after primary science targets')
    p.add_argument('--max-fillers',type=int,default=10,help='Maximum per T ALL science mask; 0 disables fillers')
    p.add_argument('--filler-g-min',type=float,default=12.0)
    p.add_argument('--filler-g-max',type=float,default=18.0)
    p.add_argument('--filler-min-ra-sep',type=float,default=1.5,
        help='Minimum projected RA separation in arcseconds')
    p.add_argument('--filler-trace-sep',type=float,default=2.4,
        help='Conservative cross-dispersion stripe separation in arcseconds')
    p.add_argument('--repeat-targets',action=argparse.BooleanOptionalAction,default=True,
        help='Fill remaining space with earlier science targets, faintest first')
    p.add_argument('--max-repeats',type=int,default=100,
        help='Maximum repeat science apertures per later T mask; 0 disables')
    p.add_argument('--n-sky',type=int,default=5);p.add_argument('--sky-exclusion-arcsec',type=float,default=8)
    p.add_argument('--sky-trace-separation-arcsec',type=float,default=2.4)
    p.add_argument('--slit-width-arcsec',type=float,default=0.54);p.add_argument('--slit-length-arcsec',type=float,default=1.62)
    p.add_argument('--fits-image',default=None,help='Actual real survey FITS with celestial WCS (required for genuine survey finding charts)')
    p.add_argument('--spectral-stripes',action=argparse.BooleanOptionalAction,default=True,
        help='Generate schematic ALL/EVEN/ODD spectral stripes for each T mask')
    p.add_argument('--spectral-length-pix',type=int,default=2883)
    p.add_argument('--seed',type=int,default=42)
    return p
if __name__=='__main__':run(parser().parse_args())
