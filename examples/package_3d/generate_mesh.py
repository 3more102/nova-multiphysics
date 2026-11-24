"""Generate a conforming, tagged flip-chip test coupon in ASCII Gmsh 2.2, SI metres."""
import argparse
from collections import Counter
from itertools import permutations
from pathlib import Path


def generate(path, level=1):
    if not 1 <= level <= 64:
        raise ValueError("subdivision multiplier must be between 1 and 64")
    nx, ny = 4 * level, 2 * level
    layers = (0, .4, .5, .8, 1., 1.1)  # mm, bottom to top
    z = [((1-i/level)*a+i/level*b)*1e-3
         for a, b in zip(layers, layers[1:]) for i in range(level)] + [layers[-1]*1e-3]
    nz = len(z)-1
    nodes = [(i*.004/nx, j*.002/ny, zk)
             for zk in z for j in range(ny+1) for i in range(nx+1)]

    def node(i, j, k):
        return 1+i+(nx+1)*(j+(ny+1)*k)

    elements, faces = [], Counter()
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                layer = k//level
                tag = (1, 2, 3, 4, 5 if i < level or i >= 3*level else 6)[layer]
                for order in permutations(range(3)):
                    xyz = [i, j, k]
                    tet = [node(*xyz)]
                    for axis in order:
                        xyz[axis] += 1
                        tet.append(node(*xyz))
                    # Positive orientation for each Freudenthal simplex.
                    if order in ((0, 2, 1), (1, 0, 2), (2, 1, 0)):
                        tet[1], tet[2] = tet[2], tet[1]
                    elements.append((4, tag, tet))
                    for skip in range(4):
                        faces[tuple(sorted(tet[:skip]+tet[skip+1:]))] += 1
    for face, count in faces.items():
        if count != 1:
            continue
        pts = [nodes[v-1] for v in face]
        tag = 15
        if all(p[2] == 0 for p in pts):
            tag = 11
        elif all(p[2] == z[-1] for p in pts):
            x = sum(p[0] for p in pts)/3
            tag = 12 if x < .001 else (13 if x > .003 else 14)
        elements.append((2, tag, face))

    names = [(3, i+1, s) for i, s in enumerate(
        ("spreader", "bond", "substrate", "die", "contact", "passivation"))]
    names += [(2, t, s) for t, s in ((11, "base"), (12, "terminal_left"),
                                       (13, "terminal_right"), (14, "lid"), (15, "side"))]
    with Path(path).open("w", encoding="ascii", newline="\n") as f:
        f.write(f"$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$PhysicalNames\n{len(names)}\n")
        for dim, tag, name in names:
            f.write(f'{dim} {tag} "{name}"\n')
        f.write(f"$EndPhysicalNames\n$Nodes\n{len(nodes)}\n")
        for i, (x, y, height) in enumerate(nodes, 1):
            f.write(f"{i} {x:.17g} {y:.17g} {height:.17g}\n")
        f.write(f"$EndNodes\n$Elements\n{len(elements)}\n")
        for i, (typ, tag, conn) in enumerate(elements, 1):
            f.write(f"{i} {typ} 2 {tag} {tag} {' '.join(map(str,conn))}\n")
        f.write("$EndElements\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh", type=Path)
    parser.add_argument("--level", type=int, default=1)
    args = parser.parse_args()
    generate(args.mesh, args.level)
