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
ACID  = {"Asp", "Glu", "D", "E"}
BASIC = {"Lys", "Arg", "His", "K", "R", "H"}
POLAR = {"Ser", "Thr", "Asn", "Gln", "Tyr", "Cys", "S", "T", "N", "Q", "Y", "C"}
APOLAR= {"Ala", "Val", "Leu", "Ile", "Pro", "Phe", "Trp", "Met", "Gly",
         "A", "V", "L", "I", "P", "F", "W", "M", "G"}

ICON_FILES = {
    "basic":  "sphere-blue.png",
    "acid":   "sphere-red.png",
    "polar":  "sphere-green.png",
    "apolar": "sphere-gray.png",
}

def normalize_residue_name(token: str) -> str:
    """
    'ASP123', 'D_Nt', 'Arg', 'd' -> 'Asp', 'D', 'Arg', 'D'
    (quita sufijos, conserva letras iniciales)
    """
    t = token.strip()
    if "_" in t:
        t = t.split("_")[0]
    # conservar sólo letras iniciales
    i = 0
    while i < len(t) and t[i].isalpha():
        i += 1
    if i > 0:
        t_letters = t[:i]
    else:
        t_letters = t
    if len(t_letters) >= 3:
        t3 = t_letters[:3]
        return t3[0].upper() + t3[1:].lower()
    elif len(t_letters) == 1:
        return t_letters.upper()
    return t_letters

def residue_category(res3: str) -> str | None:
    if res3 in ACID:   return "acid"
    if res3 in BASIC:  return "basic"
    if res3 in POLAR:  return "polar"
    if res3 in APOLAR: return "apolar"
    return None

def residue_index(token: str) -> int | None:
    """
    Extrae un índice numérico de residuo si está presente (p.ej. ASP123 -> 123).
    Si no hay número, devuelve None.
    """
    m = re.search(r'(\d+)', token)
    return int(m.group(1)) if m else None

# ==========================
# Lectura XYZ
# ==========================
def read_xyz(path):
    coords, names = [], []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    start = 0
    try:
        _ = int(lines[0].strip())  # N
        start = 2                  # comentario
    except Exception:
        start = 0

    for line in lines[start:]:
        parts = line.strip().split()
        if len(parts) < 4:
            continue
        name = parts[0]
        x, y, z = map(float, parts[1:4])
        names.append(name)
        coords.append([x, y, z])

    if not coords:
        raise ValueError(f"No se encontraron coordenadas en {path}")
    return np.asarray(coords, dtype=float), names

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

def load_icons(icons_root: str):
    icons = {}
    for cat, fname in ICON_FILES.items():
        path = os.path.join(icons_root, fname)
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Falta el icono para '{cat}': {path}")
        icons[cat] = plt.imread(path)  # PNG original SIN tinte
    return icons

def add_icon(ax, x, y, icon_rgba, size_px):
    hpx, wpx = icon_rgba.shape[:2]
    zoom = float(size_px) / float(wpx)
    oi = OffsetImage(icon_rgba, zoom=zoom)
    ab = AnnotationBbox(oi, (x, y), frameon=False)
    ax.add_artist(ab)

def autoscale(ax, xy, pad_ratio=0.05):
    xmin, ymin = np.min(xy, axis=0)
    xmax, ymax = np.max(xy, axis=0)
    dx, dy = xmax - xmin, ymax - ymin
    if dx == 0: dx = 1.0
    if dy == 0: dy = 1.0
    ax.set_xlim(xmin - pad_ratio * dx, xmax + pad_ratio * dx)
    ax.set_ylim(ymin - pad_ratio * dy, ymax + pad_ratio * dy)
    ax.set_aspect("equal", adjustable="box")

# ==========================
# MAIN
# ==========================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xyz", required=True, help="archivo .xyz (en el directorio actual)")
    ap.add_argument("--icons-root", default="beads", help="carpeta con sphere-*.png")
    ap.add_argument("--plane", default="XY", choices=["XY", "XZ", "YZ"])
    ap.add_argument("--size", type=int, default=44, help="tamaño del icono (px)")
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--outfile", default="protein.png")
    ap.add_argument("--order", choices=["sequence", "depth"], default="sequence",
                    help="sequence: dibuja siguiendo el índice/residuo en orden; depth: de lejos→cerca")
    args = ap.parse_args()

    # XYZ en el CWD
    xyz_path = args.xyz
    if not os.path.isfile(xyz_path):
        raise FileNotFoundError(f"No se encontró el archivo XYZ en el directorio actual: {xyz_path}")

    # Datos
    coords, names = read_xyz(xyz_path)
    xy, depth = choose_axes(coords, args.plane)

    # Categoría y (posible) índice de residuo
    cats = []
    rids = []
    for nm in names:
        cats.append(residue_category(normalize_residue_name(nm)))
        rids.append(residue_index(nm))
    cats = np.array(cats, dtype=object)
    rids = np.array(rids, dtype=object)

    # Orden de dibujo
    if args.order == "sequence":
        # 1) si hay número de residuo, ordenar por ese número (estable)
        # 2) si no, usar el orden del archivo
        has_idx = np.array([ri is not None for ri in rids], dtype=bool)
        order = np.arange(len(names))
        if has_idx.any():
            # ordenar por rids, y mantener estabilidad para los None (quedan como en el archivo)
            idx_with = np.where(has_idx)[0]
            idx_without = np.where(~has_idx)[0]
            # ordenar sólo los que tienen número
            idx_sorted_with = idx_with[np.argsort(rids[idx_with].astype(int), kind="stable")]
            # concatenar respetando el orden archivo de los que no tienen
            order = np.concatenate([idx_sorted_with, idx_without], axis=0)
        # si nadie tiene número, order ya es el del archivo
    else:  # depth
        order = np.argsort(depth)

    # Cargar íconos originales
    icons = load_icons(args.icons_root)

    # Figura limpia
    fig, ax = plt.subplots(figsize=(6, 6), dpi=args.dpi)
    ax.axis("off")
    fig.patch.set_alpha(0.0)
    ax.set_facecolor((1, 1, 1, 0))

    # Dibujo encadenado: UNO POR UNO siguiendo 'order'
    for i in order:
        cat = cats[i]
        if cat is None:
            # si no se reconoce la categoría, usar apolar por defecto
            cat = "apolar"
        icon = icons[cat]
        add_icon(ax, xy[i, 0], xy[i, 1], icon, size_px=args.size)

    autoscale(ax, xy, pad_ratio=0.05)
    plt.savefig(args.outfile, dpi=args.dpi, transparent=True, bbox_inches="tight", pad_inches=0)
    print(f"[OK] Guardado: {args.outfile}")

if __name__ == "__main__":
    main()

