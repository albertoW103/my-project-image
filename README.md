# Protein Draw Protein Cartoon

This project plots a 2D projection of a protein from an `.xyz` file, placing one colored bead per residue. The background is transparent and the bead size adapts to the size of the box.

---

## How it works

- Reads `.xyz` coordinates
- Recenters by center of mass
- Projects onto XY / XZ / YZ plane
- Colors each residue by category:
  - red: acidic (Asp, Glu) 
  - blue: basic (Lys, Arg, His)
  - green: polar (Ser, Thr, Asn, Gln, Tyr, Cys, CysDB)
  - gray: apolar (Ala, Val, Leu, Ile, Pro, Phe, Trp, Met, Gly)
- Draws a PNG (or SVG) bead centered on each residue
- Bead size scales according to the box size (in order to compared with other proteins)

---



The code run as following:

`python3 draw-protein-cartoon.py -xyz protein.xyz -plane XY -fmt svg` 


## Example Results

<p align="center"> <img src="figures/protein_4F5S.png" alt="BSA" width="45%"/> <img src="figures/protein_1J05.png" alt="2D5" width="45%"/> </p>

**Figure 1.** Coarce graned of proteins BSA and 2D5.







