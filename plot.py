#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Proyección 2D de un archivo .xyz colocando una imagen PNG por residuo
con color por categoría (basic, acid, polar, apolar), FONDO transparente,
y ORDEN DE DIBUJO encadenado por secuencia (cada residuo se dibuja encima del anterior).
"""

import argparse
import os
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

# ==========================
# Clasificación de residuos
# ==========================

# clasification of residues:
ACID  = {"Asp", "Glu"}
BASIC = {"Lys", "Arg", "His"}
POLAR = {"Ser", "Thr", "Asn", "Gln", "Tyr", "Cys", "CysDB"}
APOLAR= {"Ala", "Val", "Leu", "Ile", "Pro", "Phe", "Trp", "Met", "Gly"}

ICON_FILES = {
    "basic":  "sphere-basic.png",
    "acid":   "sphere-acid.png",
    "polar":  "sphere-polar.png",
    "apolar": "sphere-apolar.png"
}

amino_acid_code = {
    'A':'Ala', 'R':'Arg', 'N':'Asn', 'D':'Asp', 'C':'Cys',
    'E':'Glu', 'Q':'Gln', 'G':'Gly', 'H':'His', 'I':'Ile',
    'L':'Leu', 'K':'Lys', 'M':'Met', 'F':'Phe', 'P':'Pro',
    'S':'Ser', 'T':'Thr', 'W':'Trp', 'Y':'Tyr', 'V':'Val'
}

def residue_category(residue):
    """
    """
    if residue in ACID:
        return "acid"
    if residue in BASIC:
        return "basic"
    if residue in POLAR:
        return "polar"
    if residue in APOLAR:
        return "apolar"
    else:
        print('I need a cathegory')

def read_xyz_file(input_filename):
    """
    """
    with open(input_filename, 'r') as file:
        # number of atoms/points:
        nu = int(file.readline().strip())
        
        # header line:
        header = file.readline().rstrip("\n")
        
        # accumulators:
        indices, residues, categories = [], [], []
        
        # read N coordinate lines:
        xcm = ycm = zcm = 0.0           # for centre of mass
        x_coord, y_coord, z_coord = [], [], []
        index = 0                       # python index starts at 0
        
        for _ in range(nu):
            line    = file.readline().split()
            residue = line[0]
            x       = float(line[1])
            y       = float(line[2])
            z       = float(line[3])
            
            # convert from one-letter with suffix to three-letter if corresponde
            if ('_Nt' in residue) or ('_Ct' in residue):
                residue_one_letter_code = residue.split('_')[0]
                residue = amino_acid_code[residue_one_letter_code] 
            
            # get category:
            residue_cat = residue_category(residue)
            
            # append all information:
            residues.append(residue)
            x_coord.append(x)
            y_coord.append(y)
            z_coord.append(z)
            indices.append(index)
            categories.append(residue_cat)
            
            # acumular centro de masa
            xcm += x
            ycm += y
            zcm += z

            index += 1
        
        # get centre of mass:
        xcm /= nu
        ycm /= nu
        zcm /= nu
        
        # normalize:
        xsv = np.asarray(x_coord, dtype=float) - xcm
        ysv = np.asarray(y_coord, dtype=float) - ycm
        zsv = np.asarray(z_coord, dtype=float) - zcm
        
        # chequeo mínimo
        if index != nu:
            print('missing residues')
            exit(1)
    
    # coords como numpy array:
    coords = np.column_stack((xsv, ysv, zsv))
    return indices, residues, categories, coords

# ==========================
# Proyección / plot
# ==========================
def choose_axes(coords, plane):
    p = plane.upper()
    if p == "XY":
        return coords[:, [0, 1]], coords[:, 2]
    elif p == "XZ":
        return coords[:, [0, 2]], coords[:, 1]
    elif p == "YZ":
        return coords[:, [1, 2]], coords[:, 0]
    else:
        raise ValueError("Plano inválido. Use XY, XZ o YZ.")

def load_icons():
    '''
    '''
    icons = {}
    for cat, fname in ICON_FILES.items():
        path = os.path.join("beads", fname)  # <--- acá se fija en beads/
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Falta el icono: {path}")
        icons[cat] = plt.imread(path)
    return icons

def add_icon(ax, x, y, icon_rgba, size_px):
    '''
    size_px ahora se interpreta en UNIDADES DE DATOS (ancho del icono),
    para que escale con la caja fija y no con píxeles de la figura.
    '''
    hpx, wpx = icon_rgba.shape[:2]
    aspect = hpx / float(wpx)
    half_w = 0.5 * float(size_px)
    half_h = half_w * aspect
    extent = [x - half_w, x + half_w, y - half_h, y + half_h]
    ax.imshow(icon_rgba, extent=extent, origin='upper',
              interpolation='antialiased')

# ==========================
# MAIN
# ==========================
def main():
    '''
    '''
    # read arguments from terminal:
    ap = argparse.ArgumentParser()
    ap.add_argument("--xyz", required=True, help="file .xyz (into the current directory)")
    ap.add_argument("--plane", default="XY", choices=["XY", "XZ", "YZ"], help="projection plane")
    ap.add_argument("--outfile", default="protein.png")
    args = ap.parse_args()
    
    # read xyz file:
    indices, residues, categories, coords = read_xyz_file(args.xyz)
    xy, depth = choose_axes(coords, args.plane)
        
    # tamaño de esfera COMO:
    box     = 10.0
    size_px =  0.45
    #size_px = 0.10 * (40.0 / box)
    
    # get icons:
    icons = load_icons()
    
    # plot:
    fig, ax = plt.subplots(figsize=(6, 6), dpi=300)
    ax.axis("off")
    fig.patch.set_alpha(0.0)
    ax.set_facecolor((1, 1, 1, 0))
    
    # plot in sequence order:
    for idx in indices:
        residue_cat = categories[idx]   # select cathegory
        icon = icons[residue_cat]       # select icon
        add_icon(ax, xy[idx, 0], xy[idx, 1], icon, size_px=size_px)  # add icon
    
    ax.set_xlim(-box/2, box/2)
    ax.set_ylim(-box/2, box/2)
    plt.savefig(args.outfile, dpi=300, transparent=True)

if __name__ == "__main__":
    main()

