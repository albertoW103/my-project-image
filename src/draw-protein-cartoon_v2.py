#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D projection of a .xyz file placing a PNG icon per residue category
with transparent background and ordered drawing in sequence order.
Ahora: rotación (Euler ZYZ) sin cambiar tamaño de imagen.
- Rotación con mouse (drag) y/o teclado (flechas/teclas).
- Sin opción 'plane' (usa XY por defecto).
"""

import argparse
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from scipy.spatial.transform import Rotation as R

# ==========================
# Residue classification
# ==========================

ACID  = {"Asp", "Glu"}
BASIC = {"Lys", "Arg", "His"}
POLAR = {"Ser", "Thr", "Asn", "Gln", "Tyr", "Cys", "CysDB"}
APOLAR= {"Ala", "Val", "Leu", "Ile", "Pro", "Phe", "Trp", "Met", "Gly"}

amino_acid_code = {
    'A':'Ala', 'R':'Arg', 'N':'Asn', 'D':'Asp', 'C':'Cys',
    'E':'Glu', 'Q':'Gln', 'G':'Gly', 'H':'His', 'I':'Ile',
    'L':'Leu', 'K':'Lys', 'M':'Met', 'F':'Phe', 'P':'Pro',
    'S':'Ser', 'T':'Thr', 'W':'Trp', 'Y':'Tyr', 'V':'Val'
}

def residue_category(residue):
    if residue in ACID:  return "acid"
    if residue in BASIC: return "basic"
    if residue in POLAR: return "polar"
    if residue in APOLAR:return "apolar"
    print('I need a category')

def read_xyz_file(input_filename):
    with open(input_filename, 'r') as file:
        nu = int(file.readline().strip())
        _ = file.readline()  # header

        indices, residues, categories = [], [], []
        xcm = ycm = zcm = 0.0
        x_coord, y_coord, z_coord = [], [], []
        index = 0

        for _ in range(nu):
            line    = file.readline().split()
            residue = line[0]
            x       = float(line[1]); y = float(line[2]); z = float(line[3])

            if ('_Nt' in residue) or ('_Ct' in residue):
                residue = amino_acid_code[residue.split('_')[0]]

            residue_cat = residue_category(residue)

            residues.append(residue)
            x_coord.append(x); y_coord.append(y); z_coord.append(z)
            indices.append(index); categories.append(residue_cat)

            xcm += x; ycm += y; zcm += z
            index += 1

        xcm /= nu; ycm /= nu; zcm /= nu

        xsv = np.asarray(x_coord, dtype=float) - xcm
        ysv = np.asarray(y_coord, dtype=float) - ycm
        zsv = np.asarray(z_coord, dtype=float) - zcm

        if index != nu:
            print('missing residues'); exit(1)

    coords = np.column_stack((xsv, ysv, zsv))
    return indices, residues, categories, coords

def project_XY(coords):
    """ Proyección XY por defecto (sin opción plane). """
    return coords[:, [0, 1]], coords[:, 2]

def load_icons():
    icon_files = {
        "basic":  "sphere-basic.png",
        "acid":   "sphere-acid.png",
        "polar":  "sphere-polar.png",
        "apolar": "sphere-apolar.png"
    }
    icons = {}
    for category, filename in icon_files.items():
        path = os.path.join("../beads", filename)
        icons[category] = plt.imread(path)
    return icons

def add_icon(ax, x, y, icon_rgba, size_data, zorder=0):
    # ancho objetivo en píxeles para un ancho 'size_data' en unidades de datos:
    p0 = ax.transData.transform((x, y))
    p1 = ax.transData.transform((x + float(size_data), y))
    px_target = abs(p1[0] - p0[0])

    wpx = float(icon_rgba.shape[1])
    zoom = (px_target / wpx) if wpx > 0 else 1.0

    oi = OffsetImage(icon_rgba, zoom=zoom)
    ab = AnnotationBbox(oi, (x, y),
                        xycoords='data',
                        frameon=False,
                        box_alignment=(0.5, 0.5))
    ab.set_zorder(zorder)
    ax.add_artist(ab)

# ==========================
# MAIN
# ==========================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-xyz", required=True, help="input .xyz file")
    ap.add_argument("-fmt", default="png", choices=["png","svg"], help="output format")
    ap.add_argument("-mouse", action="store_true", help="Habilitar rotación con mouse (drag).")
    ap.add_argument("-keyboard", action="store_true", help="Habilitar rotación con teclado.")
    args = ap.parse_args()

    # Si no eligen nada, por defecto mouse:
    if not args.mouse and not args.keyboard:
        args.mouse = True

    indices, residues, categories, coords0 = read_xyz_file(args.xyz)

    # ángulos internos (Euler ZYZ)
    theta = 0.0; phi = 0.0; psi = 0.0

    # Sensibilidades
    SENS_MOUSE = 0.01         # rad por pixel
    STEP_PHI   = 0.05         # rad por tecla
    STEP_THETA = 0.05         # rad por tecla
    STEP_PSI   = 0.05         # rad por tecla

    # rota simplificada (sin args): usa los ángulos actuales
    def rota(xsv, ysv, zsv):
        coords = np.column_stack((xsv, ysv, zsv))
        rotation = R.from_euler('ZYZ', [phi, theta, psi], degrees=False)
        cr = rotation.apply(coords)
        return cr[:,0], cr[:,1], cr[:,2]

    # parámetros de dibujo (idénticos a tu script)
    REF_BOX_NM     = 10.0
    REF_SPHERE_NM  = 0.03
    box       = 10.0
    size_data = REF_SPHERE_NM * (REF_BOX_NM / box)

    icons = load_icons()

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.axis("off")
    fig.patch.set_alpha(0.0)
    ax.set_facecolor((1, 1, 1, 0))

    def redraw():
        ax.clear()
        ax.axis("off")
        ax.set_facecolor((1,1,1,0))
        xr, yr, zr = rota(coords0[:,0], coords0[:,1], coords0[:,2])
        coords = np.column_stack((xr, yr, zr))
        xy, _ = project_XY(coords)
        for idx in indices:
            cat = categories[idx]
            icon = icons.get(cat, icons['apolar'])
            add_icon(ax, xy[idx, 0], xy[idx, 1], icon, size_data, zorder=10+idx)
        ax.set_xlim(-box/2, box/2)
        ax.set_ylim(-box/2, box/2)
        fig.canvas.draw_idle()

    redraw()

    # --------------------
    # MODO MOUSE (drag)
    # --------------------
    if args.mouse:
        drag = {"active": False, "x": 0, "y": 0}

        def on_press(event):
            if event.inaxes != ax: return
            drag["active"] = True
            drag["x"], drag["y"] = event.x, event.y

        def on_release(event):
            drag["active"] = False

        def on_move(event):
            nonlocal phi, theta
            if not drag["active"] or event.inaxes != ax: return
            dx = event.x - drag["x"]
            dy = event.y - drag["y"]
            drag["x"], drag["y"] = event.x, event.y
            # actualizar ZYZ: dx ↦ phi (alrededor de Z), dy ↦ theta
            phi   += dx * SENS_MOUSE
            theta += dy * SENS_MOUSE
            redraw()

        fig.canvas.mpl_connect('button_press_event', on_press)
        fig.canvas.mpl_connect('button_release_event', on_release)
        fig.canvas.mpl_connect('motion_notify_event', on_move)

    # --------------------
    # MODO TECLADO
    # --------------------
    if args.keyboard:
        def on_key(event):
            nonlocal phi, theta, psi
            need_redraw = False
            if event.key in ('left', 'a'):
                phi -= STEP_PHI;   need_redraw = True
            elif event.key in ('right', 'd'):
                phi += STEP_PHI;   need_redraw = True
            elif event.key in ('up', 'w'):
                theta -= STEP_THETA; need_redraw = True
            elif event.key in ('down', 's'):
                theta += STEP_THETA; need_redraw = True
            elif event.key in ('q',):
                psi -= STEP_PSI;   need_redraw = True
            elif event.key in ('e',):
                psi += STEP_PSI;   need_redraw = True
            elif event.key and event.key.lower() == 'r':
                phi = theta = psi = 0.0; need_redraw = True
            elif event.key and event.key.lower() == 's':
                plt.savefig("protein.%s" % args.fmt, dpi=300, transparent=True)
                print("[INFO] Saved current view.")
            if need_redraw:
                redraw()

        fig.canvas.mpl_connect('key_press_event', on_key)

    # Info de ayuda
    print("[INFO] Rotación:", end=" ")
    if args.mouse:    print("mouse (drag)", end="")
    if args.mouse and args.keyboard: print(" + ", end="")
    if args.keyboard: print("teclado (← → ↑ ↓, A/D/W/S, Q/E, R reset, S guardar)", end="")
    print()

    plt.show()

if __name__ == "__main__":
    main()

