"""Gera cubos (.cube) de HOMO, LUMO e do par NTO dominante a partir de resultados/orbitais/<id>.npz (M3).

Não refaz o cálculo: reconstrói a molécula (símbolos, coordenadas, base) e usa os coeficientes salvos.
Estados exportados: S1 e o estado mais intenso acima de 250 nm.

  "<venv>/bin/python" orbitais_cubos.py [id ...] [--pontos 60]   → resultados/cubos/<id>__<orbital>.cube
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from pyscf import gto
from pyscf.tools import cubegen

from quimiolib.qm import nto_par

AQUI = Path(__file__).resolve().parent
ORB = AQUI / "resultados" / "orbitais"
CUB = AQUI / "resultados" / "cubos"


def cubos(npz: Path, pontos: int = 60) -> list[str]:
    z = np.load(npz)
    mol = gto.M(atom=[(str(s), tuple(c)) for s, c in zip(z["simbolos"], z["coords_ang"])],
                basis=str(z["base"]), verbose=0)
    occ = z["mo_occ"] > 0
    ihomo = int(np.where(occ)[0][-1])
    C = z["mo_coeff"]
    orbo, orbv = C[:, occ], C[:, ~occ]
    lam = 1239.84193 / z["e_ev"]
    estados = {"S1": 0}
    acima = [i for i in range(len(lam)) if lam[i] >= 250]
    if acima:
        estados["Sforte"] = max(acima, key=lambda i: z["f"][i])
    alvos = {"HOMO": C[:, ihomo], "LUMO": C[:, ihomo + 1]}
    meta = {"homo_idx": ihomo, "estados": {}}
    for nome, k in estados.items():
        w, o, v = nto_par(z["X"][k], orbo, orbv)
        alvos[f"{nome}_NTO_ocupado"], alvos[f"{nome}_NTO_virtual"] = o, v
        meta["estados"][nome] = {"estado": k + 1, "lambda_nm": round(float(lam[k]), 1),
                                 "f": round(float(z["f"][k]), 4), "nto_peso": round(w, 3)}
    CUB.mkdir(parents=True, exist_ok=True)
    feitos = []
    for nome, coef in alvos.items():
        arq = CUB / f"{npz.stem}__{nome}.cube"
        cubegen.orbital(mol, str(arq), coef, nx=pontos, ny=pontos, nz=pontos)
        feitos.append(arq.name)
    (CUB / f"{npz.stem}__meta.json").write_text(json.dumps(meta, indent=1))
    return feitos


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("ids", nargs="*")
    ap.add_argument("--pontos", type=int, default=60)
    a = ap.parse_args()
    arqs = [ORB / f"{i}.npz" for i in a.ids] if a.ids else sorted(ORB.glob("*.npz"))
    for f in arqs:
        print(f.stem, cubos(f, a.pontos))
