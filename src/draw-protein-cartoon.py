#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D projection of a .xyz file placing a PNG icon per residue category
with transparent background and ordered drawing in sequence order.
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
# Residue classification
# ==========================

# categories of residues:
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
    """
    Returns the residue category (acid, basic, polar or apolar)
    according to the dictionaries defined above.
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
        print('I need a category')

def read_xyz_file(input_filename):
    """
    Reads a .xyz file, extracts residue names and coordinates,
    and returns indices, residue categories and coordinates
    (recentred at the center of mass).
    """
    with open(input_filename, 'r') as file:
        # number of atoms/points:
        nu = int(file.readline().strip())
        
        # read header (ignored content, but keeps file format)
        header = file.readline().rstrip("\n")
        
        # accumulators:
        indices, residues, categories = [], [], []
        
        # read coordinate lines:
        xcm = ycm = zcm = 0.0               # centre of mass accumulators
        x_coord, y_coord, z_coord = [], [], []
        index = 0
        
        for _ in range(nu):
            line    = file.readline().split()
            residue = line[0]
            x       = float(line[1])
            y       = float(line[2])
            z       = float(line[3])
            
            # convert one-letter + suffix to 3-letter name when necessary
            if ('_Nt' in residue) or ('_Ct' in residue):
                residue_one_letter_code = residue.split('_')[0]
                residue = amino_acid_code[residue_one_letter_code] 
            
            # classify residue
            residue_cat = residue_category(residue)
            
            # append:
            residues.append(residue)
            x_coord.append(x)
            y_coord.append(y)
            z_coord.append(z)
            indices.append(index)
            categories.append(residue_cat)
            
            # accumulate for the center of mass
            xcm += x
            ycm += y
            zcm += z

            index += 1
        
        # compute center of mass:
        xcm /= nu
        ycm /= nu
        zcm /= nu
        
        # recenter coords:
        xsv = np.asarray(x_coord, dtype=float) - xcm
        ysv = np.asarray(y_coord, dtype=float) - ycm
        zsv = np.asarray(z_coord, dtype=float) - zcm
        
        # quick consistency check
        if index != nu:
            print('missing residues')
            exit(1)
    
    # pack in numpy array:
    coords = np.column_stack((xsv, ysv, zsv))
    return indices, residues, categories, coords

# ==========================
# Projection utilities
# ==========================
def choose_axes(coords, plane):
    """
    Reorders the coordinate axes according to the selected projection plane (XY, XZ or YZ).
    """
    p = plane.upper()
    if p == "XY":
        return coords[:, [0, 1]], coords[:, 2]
    elif p == "XZ":
        return coords[:, [0, 2]], coords[:, 1]
    elif p == "YZ":
        return coords[:, [1, 2]], coords[:, 0]
    else:
        raise ValueError("Invalid plane. Use XY, XZ or YZ.")

def load_icons():
    '''
    Loads the PNG icons for each residue category and returns a dictionary
    mapping category -> image array.
    '''
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
    """
    Draws the icon centered at (x, y) as a sprite using annotation,
    scaling the icon according to the desired size in data units.
    """
    # original PNG width in pixels
    hpx, wpx = icon_rgba.shape[:2]

    # compute how many display pixels correspond to `size_data` data units
    p0 = ax.transData.transform((x, y))
    p1 = ax.transData.transform((x + float(size_data), y))
    px_target = abs(p1[0] - p0[0])

    # compute zoom factor
    zoom = (px_target / float(wpx)) if wpx > 0 else 1.0

    oi = OffsetImage(icon_rgba, zoom=zoom)
    ab = AnnotationBbox(oi, (x, y),
                        xycoords='data',
                        frameon=False,
                        box_alignment=(0.5, 0.5))
    ab.set_zorder(zorder)
    ax.add_artist(ab)

# ==========================
# MAIN PROGRAM
# ==========================
def main():
    '''
    Main driver function: parses arguments, reads xyz, projects
    and draws residue category icons on a fixed-size box.
    '''
    # read arguments from terminal:
    ap = argparse.ArgumentParser()
    ap.add_argument("-xyz", required=True, help="input .xyz file")
    ap.add_argument("-plane", default="XY", choices=["XY", "XZ", "YZ"], help="projection plane")
    args = ap.parse_args()
    
    # read xyz file:
    indices, residues, categories, coords = read_xyz_file(args.xyz)
    xy, depth = choose_axes(coords, args.plane)
    
    # sphere size relative to the chosen box size:
    REF_BOX_NM     = 10.0
    REF_SPHERE_NM  = 0.03
    
    box       = 10.0                               # [default] in nm
    size_data = REF_SPHERE_NM * (REF_BOX_NM / box) # [default] in nm
    
    # load icons:
    icons = load_icons()
    
    # create figure:
    fig, ax = plt.subplots(figsize=(10, 10))  # lienzo rectangular
    ax.axis("off")
    fig.patch.set_alpha(0.0)          # make background fully transparent
    ax.set_facecolor((1, 1, 1, 0))    # make background fully transparent
    
    # plot each residue icon in sequence order:
    for idx in indices:
        residue_cat = categories[idx]
        icon = icons[residue_cat]
        add_icon(ax, xy[idx, 0], xy[idx, 1], icon, size_data, zorder=10+idx)
    
    ax.set_xlim(-box/2, box/2)
    ax.set_ylim(-box/2, box/2)
    plt.savefig('protein.png', dpi=50, transparent=True)

if __name__ == "__main__":
    main()

