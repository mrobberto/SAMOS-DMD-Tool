"""Compare pre-finalization and finalized SAMOS science region files by Gaia target ID."""
import argparse
import re
from pathlib import Path
import numpy as np

BOX = re.compile(r'^\s*box\(\s*([+-]?[\d.]+)\s*,\s*([+-]?[\d.]+)\s*,\s*([\d.]+)"\s*,\s*([\d.]+)"\s*,\s*([+-]?[\d.]+)\s*\)(.*)$', re.I)
LABEL = re.compile(r'text=\{([^}]+)\}')

def read_regions(path):
    rows=[]
    for line in Path(path).read_text().splitlines():
        m=BOX.match(line)
        if not m: continue
        lab=LABEL.search(m.group(6))
        rows.append(dict(RA=float(m.group(1)),DEC=float(m.group(2)),ra_width=float(m.group(3)),dec_height=float(m.group(4)),angle=float(m.group(5)),label=lab.group(1) if lab else ''))
    return rows

def audit(old_region,final_region,manifest, tolerance_mas=1.0):
    import pandas as pd
    old=read_regions(old_region);new=read_regions(final_region)
    manifest=pd.read_csv(manifest,dtype={'target_id':str,'SLITID':str})
    # labels in old regions are Gaia target IDs, labels in final regions are SLITIDs.
    m=manifest[manifest.SLITTYPE=='SCIENCE'].copy()
    mapping=dict(zip(m.SLITID.astype(str),m.target_id.astype(str)))
    a={str(r['label']):r for r in old}
    b={mapping[r['label']]:r for r in new if r['label'] in mapping}
    mismatches=[]
    for sid in sorted(set(a)&set(b)):
        ar,br=a[sid],b[sid]
        dra=(ar['RA']-br['RA'])*3600*1000*np.cos(np.deg2rad(ar['DEC']))
        ddec=(ar['DEC']-br['DEC'])*3600*1000
        if np.hypot(dra,ddec)>tolerance_mas:
            mismatches.append((sid,round(dra,3),round(ddec,3)))
    return {'matched':len(set(a)&set(b)),'old_only':sorted(set(a)-set(b)), 'final_only':sorted(set(b)-set(a)), 'coord_mismatches_mas':mismatches,'old_size':sorted({(r['ra_width'],r['dec_height']) for r in old}),'final_size':sorted({(r['ra_width'],r['dec_height']) for r in new})}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original-reg',required=True);p.add_argument('--final-reg',required=True);p.add_argument('--manifest',required=True);a=p.parse_args()
    report=audit(a.original_reg,a.final_reg,a.manifest)
    for k,v in report.items():print(k,':',v)
    if report['coord_mismatches_mas']:raise SystemExit('ERROR: Science-target coordinates changed')
if __name__=='__main__':main()
