#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, os, numpy as np, matplotlib
# backend interactivo
for _b in ("Qt5Agg","QtAgg","TkAgg","GTK3Agg","MacOSX"):
    try: matplotlib.use(_b); break
    except Exception: pass

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import proj3d
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

# --- categorías
ACID  = {"Asp","Glu"}
BASIC = {"Lys","Arg","His"}
POLAR = {"Ser","Thr","Asn","Gln","Tyr","Cys","CysDB"}
APOLAR= {"Ala","Val","Leu","Ile","Pro","Phe","Trp","Met","Gly"}
amino_acid_code = {'A':'Ala','R':'Arg','N':'Asn','D':'Asp','C':'Cys','E':'Glu','Q':'Gln','G':'Gly','H':'His','I':'Ile',
                   'L':'Leu','K':'Lys','M':'Met','F':'Phe','P':'Pro','S':'Ser','T':'Thr','W':'Trp','Y':'Tyr','V':'Val'}

def residue_category(r):
    if r in ACID: return "acid"
    if r in BASIC: return "basic"
    if r in POLAR: return "polar"
    if r in APOLAR: return "apolar"
    return "apolar"

def read_xyz_file(fname):
    with open(fname,'r') as f:
        nu = int(f.readline().strip()); _ = f.readline()
        idx, residues, cats = [], [], []
        xs, ys, zs = [], [], []; xcm=ycm=zcm=0.0
        for i in range(nu):
            s = f.readline().split()
            res, x, y, z = s[0], float(s[1]), float(s[2]), float(s[3])
            if ('_Nt' in res) or ('_Ct' in res): res = amino_acid_code[res.split('_')[0]]
            idx.append(i); residues.append(res); cats.append(residue_category(res))
            xs.append(x); ys.append(y); zs.append(z); xcm+=x; ycm+=y; zcm+=z
        if len(idx)!=nu: raise SystemExit("missing residues")
        xcm/=nu; ycm/=nu; zcm/=nu
        xyz = np.column_stack((np.asarray(xs)-xcm, np.asarray(ys)-ycm, np.asarray(zs)-zcm))
    return idx, residues, cats, xyz

def load_icons():
    files = {"basic":"sphere-basic.png","acid":"sphere-acid.png","polar":"sphere-polar.png","apolar":"sphere-apolar.png"}
    icons={}
    for k,fn in files.items():
        p = os.path.join("../beads", fn)
        icons[k] = plt.imread(p)
    return icons

class BillboardOverlay:
    """
    Pega PNGs en un AXES 2D superpuesto con coords en PIXELES.
    Se actualiza en cada draw para que los sprites miren a cámara y sigan la proyección.
    """
    def __init__(self, fig, ax3d, icons, coords, categories, size_data=0.03, zorder=1000):
        self.fig   = fig
        self.ax3d  = ax3d
        self.icons = icons
        self.xyz   = coords
        self.cats  = categories
        self.size  = float(size_data)
        # overlay 2D en píxeles
        self.ax2d  = fig.add_axes([0,0,1,1], zorder=zorder)
        self.ax2d.set_axis_off()
        self.art   = []
        # conectar eventos
        self.cid = fig.canvas.mpl_connect('draw_event', self._on_draw)

    def clear(self):
        for a in self.art:
            try: a.remove()
            except Exception: pass
        self.art.clear()

    def _on_draw(self, event):
        # actualizar límites del eje 2D a tamaño actual del canvas (en px)
        bb = self.fig.bbox
        self.ax2d.set_xlim(0, bb.width); self.ax2d.set_ylim(bb.height, 0)  # origen arriba-izq
        # borrar sprites previos
        self.clear()
        M = self.ax3d.get_proj()
        # dibujar cada PNG
        for (x,y,z), cat in zip(self.xyz, self.cats):
            x2, y2, _ = proj3d.proj_transform(x, y, z, M)
            # data->px del eje 3D
            xp, yp = self.ax3d.transData.transform((x2, y2))
            # zoom en px: proyectar un paso 'size' en x
            x2dx, y2dx, _ = proj3d.proj_transform(x+self.size, y, z, M)
            xp2, yp2 = self.ax3d.transData.transform((x2dx, y2dx))
            px_target = max(8.0, abs(xp2 - xp))  # ancho efectivo en px
            im = self.icons.get(cat, self.icons['apolar'])
            wpx = float(im.shape[1]) if im.ndim>=2 else 32.0
            zoom = max(0.15, px_target / max(1.0, wpx))
            oi = OffsetImage(im, zoom=zoom)
            ab = AnnotationBbox(oi, (xp, yp),
                                xycoords=self.ax2d.transData,
                                frameon=False, box_alignment=(0.5,0.5),
                                zorder=2000)
            ab.set_clip_on(False)
            self.ax2d.add_artist(ab); self.art.append(ab)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-xyz", required=True, help="input .xyz file")
    ap.add_argument("--out", default="protein", help="output name (no ext)")
    ap.add_argument("-fmt", default="png", choices=["png","svg"])
    args = ap.parse_args()

    idx, residues, cats, coords = read_xyz_file(args.xyz)
    icons = load_icons()

    fig = plt.figure(figsize=(10, 10))
    ax  = fig.add_subplot(111, projection='3d')
    ax.set_axis_off(); ax.set_facecolor((1,1,1,0)); fig.patch.set_alpha(0)

    # límites cúbicos + margen
    mins = coords.min(axis=0); maxs = coords.max(axis=0)
    span = (maxs - mins).max() or 1.0; center = (mins + maxs)/2.0; m = 0.1*span
    ax.set_xlim(center[0]-0.5*span-m, center[0]+0.5*span+m)
    ax.set_ylim(center[1]-0.5*span-m, center[1]+0.5*span+m)
    ax.set_zlim(center[2]-0.5*span-m, center[2]+0.5*span+m)
    try: ax.set_box_aspect((1,1,1))
    except Exception: pass
    ax.set_autoscale_on(False)

    # puntos de referencia (quita si no los querés)
    ax.scatter(coords[:,0], coords[:,1], coords[:,2],
               s=5, c='#BBBBBB', depthshade=False, zorder=5)

    # overlay de PNGs (SIEMPRE FRENTE y SOBRE CADA PUNTO)
    overlay = BillboardOverlay(fig, ax, icons, coords, cats, size_data=0.03)

    # primer draw para ver sprites
    fig.canvas.draw_idle(); plt.pause(0.05)

    def on_key(e):
        if e.key and e.key.lower()=='s':
            plt.savefig(f"{args.out}.{args.fmt}", dpi=300, transparent=True)
            print(f"[INFO] Saved {args.out}.{args.fmt}")

    fig.canvas.mpl_connect('key_press_event', on_key)
    print(f"[INFO] Backend: {matplotlib.get_backend()}")
    print("[INFO] Rotá con el mouse. Tecla 'S' para guardar.")
    plt.show()

if __name__ == "__main__":
    main()

