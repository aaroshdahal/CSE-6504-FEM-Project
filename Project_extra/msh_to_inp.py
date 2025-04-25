#!/usr/bin/env python3
"""
convert_msh_to_inp.py

Reads a Gmsh .msh and writes your FEM input format with sections:
  * COORDINATES
  * ELEMENTS
  * SECTIONS
  * TRACTION FORCES
  * BOUNDARY

Usage:
    python convert_msh_to_inp.py  <in.msh>  <out.txt>
"""
import meshio
import sys, os

def main(mesh_file, out_file):
    mesh = meshio.read(mesh_file)
    pts2d = mesh.points[:, :2]

    # pick triangles if present, else quads
    if "triangle" in mesh.cells_dict:
        elems = mesh.cells_dict["triangle"]
        elem_keyword = "TRI"
    elif "quad" in mesh.cells_dict:
        elems = mesh.cells_dict["quad"]
        elem_keyword = "QUAD"
    else:
        raise RuntimeError("No TRI or QUAD cells found in mesh.")

    # start writing
    with open(out_file, "w") as f:
        # INFO
        f.write("* INFO\n")
        f.write(f"[{elem_keyword}] imported from {os.path.basename(mesh_file)}\n\n")

        # ANALYSIS TYPE (you can change this)
        f.write("* ANALYSIS TYPE\n")
        f.write("Plane Stress\n\n")

        # COORDINATES
        f.write("* COORDINATES\n")
        for i,(x,y) in enumerate(pts2d, start=1):
            f.write(f"{i:<4d} {x:>10.6f} {y:>10.6f}\n")
        f.write("\n")

        # ELEMENTS
        f.write("* ELEMENTS\n")
        for i, conn in enumerate(elems, start=1):
            # Gmsh→0-based, your format→1-based
            conn1 = conn + 1
            if elem_keyword=="TRI":
                f.write(f"{i:<4d}  1   {conn1[0]}   {conn1[1]}   {conn1[2]}\n")
            else:  # QUAD
                f.write(f"{i:<4d}  1   {conn1[0]}   {conn1[1]}   {conn1[2]}   {conn1[3]}\n")
        f.write("\n")

        # SECTIONS (stub — fill E, ν, t as needed)
        f.write("* SECTIONS\n")
        f.write("1   210000.0   0.3   1.0\n\n")

        # TRACTION FORCES (stub — can be extended to parse Physical groups)
        f.write("* TRACTION FORCES\n")
        # e.g. f.write("2   1   2.5   0.0\n")
        f.write("\n")

        # BOUNDARY (stub — similarly, parse tags if you like)
        f.write("* BOUNDARY\n")
        # e.g. f.write("1   1   0.0\n4   2   0.0\n")

    print(f"Wrote converted input → {out_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_msh_to_inp.py <in.msh> <out.txt>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
