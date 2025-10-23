# Protein 2D Projection (Residue-Colored Beads)

This repository provides a small visualization utility to project a 3D protein structure stored in `.xyz` format onto a 2D plane (XY, XZ or YZ), placing a colored bead **per residue**. The beads are categorized by amino acid type (acidic, basic, polar, apolar) and drawn as transparent PNG icons.

Unlike pixel-based approaches, bead size is controlled in **data units (nm)** so that it remains physically meaningful relative to the size of the protein/domain. This ensures consistent scaling when switching between different protein sizes.

---

## ✨ Features

- ⚛️ One bead per residue
- 🎨 Automatic coloring by category (acid / basic / polar / apolar)
- 📦 Transparent background (ready for overlays or figure assembly)
- 📐 Real-size scaling in nanometers (not DPI-dependent)
- 🔄 Choice of projection plane: XY / XZ / YZ
- 🧭 Recentered at center-of-mass for consistent visualization
- ✅ Works with any protein size

---

## 📁 Repository structure


