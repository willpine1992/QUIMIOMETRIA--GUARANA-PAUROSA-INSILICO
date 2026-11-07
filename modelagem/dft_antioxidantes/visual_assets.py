"""Recursos visuais do painel: estruturas 2D (RDKit) com o sítio ativo destacado, geometrias 3D de todas as espécies
e volumes (orbitais, densidade de spin, MEP) compactados para o visualizador 3D (3Dmol.js).

Volumes: grade 48³ (margem 4 bohr), valores recortados e quantizados em int8 (base64), gravados como
edson/painel_assets/vol/<chave>.js  →  window.VOL["<chave>"] = {...}. O painel carrega cada arquivo sob demanda
com <script src>, o que funciona abrindo o HTML direto do disco (file://).

Fontes
  M3  modelagem/tddft_uv/resultados/orbitais/<mol>.npz   (RKS CAM-B3LYP/6-31+G(d), etanol) → HOMO, LUMO, densidade+MEP
  M2b modelagem/dft_antioxidantes/especies/orbitais/<espécie>.npz (UKS M06-2X/6-311+G(d,p), metanol) → spin
"""
from __future__ import annotations

import base64
import json
import re
from pathlib import Path

import numpy as np

from dft_antioxidantes import AQUI, ESP, MOLS

EDSON = AQUI.parents[1]
ASSETS = EDSON / "painel_assets"
VOL = ASSETS / "vol"
M3 = EDSON / "modelagem" / "tddft_uv"
N = 48
CLIP = {"orb": 0.15, "spin": 0.02, "dens": 0.02, "mep": 0.10}   # au; recorte antes da quantização


# ------------------------------------------------------------------ 2D
def svg_2d(mol_id: str, bde_sitio: dict[str, float]) -> str:
    """SVG 2D da molécula; cada sítio X–H calculado ganha um halo cuja intensidade cresce quanto menor a BDE,
    e uma nota com o valor. Traços em currentColor (segue o tema do painel)."""
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D

    smi, rot, _ = MOLS[mol_id]
    m = Chem.MolFromSmiles(smi)
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(m)
    hl, cores, raios = [], {}, {}
    vals = [v for v in bde_sitio.values() if v is not None]
    lo, hi = (min(vals), max(vals)) if vals else (0, 1)
    for a in m.GetAtoms():
        mp = a.GetAtomMapNum()
        if mp:
            nome = rot.get(mp)
            v = bde_sitio.get(nome)
            a.SetAtomMapNum(0)
            if v is None:
                continue
            t = 1.0 if hi == lo else (hi - v) / (hi - lo)          # 1 = menor BDE (melhor doador)
            hl.append(a.GetIdx())
            cores[a.GetIdx()] = (0.165, 0.471, 0.839, 0.18 + 0.55 * t)   # azul s1 com alfa
            raios[a.GetIdx()] = 0.45 + 0.25 * t
            a.SetProp("atomNote", f"{v:.1f}".replace(".", ","))
    d = rdMolDraw2D.MolDraw2DSVG(420, 300)
    o = d.drawOptions()
    o.clearBackground = False
    o.annotationFontScale = 0.95
    o.minFontSize = 14
    o.bondLineWidth = 2
    o.padding = 0.12
    o.useBWAtomPalette()
    d.DrawMolecule(m, highlightAtoms=hl, highlightAtomColors=cores, highlightAtomRadii=raios, highlightBonds=[])
    d.FinishDrawing()
    s = d.GetDrawingText()
    s = re.sub(r"<\?xml[^>]*>\s*", "", s)
    s = s.replace("#000000", "currentColor").replace("#000", "currentColor")
    s = re.sub(r"<rect style='opacity:1.0;fill:#FFFFFF[^>]*>", "", s)
    return s


# ------------------------------------------------------------------ 3D
def xyz_str(simb, xyz, titulo=""):
    return f"{len(simb)}\n{titulo}\n" + "\n".join(f"{s} {x:.4f} {y:.4f} {z:.4f}" for s, (x, y, z) in zip(simb, xyz))


def _mol(z):
    from pyscf import gto
    ch = int(z["charge"]) if "charge" in z.files else 0
    sp = int(z["spin"]) if "spin" in z.files else 0
    return gto.M(atom=[(str(s), tuple(c)) for s, c in zip(z["simbolos"], z["coords_ang"])],
                 basis=str(z["base"]), charge=ch, spin=sp, verbose=0)


def _grade(mol):
    from pyscf.tools import cubegen
    cc = cubegen.Cube(mol, N, N, N, margin=4.0)
    return cc, cc.get_coords()


def _pacote(chave, tipo, cc, mol, valores, meta):
    lim = CLIP[tipo]
    q = np.clip(np.round(np.asarray(valores).reshape(N, N, N) / lim * 127), -127, 127).astype(np.int8)
    pac = {"n": N, "origem": cc.boxorig.tolist(), "passo": [float(cc.box[i, i] / N) for i in range(3)],
           "escala": lim / 127, "tipo": tipo, "z": [int(mol.atom_charge(i)) for i in range(mol.natm)],
           "xyz_bohr": mol.atom_coords().round(4).tolist(), "dados": base64.b64encode(q.tobytes()).decode(),
           **meta}
    VOL.mkdir(parents=True, exist_ok=True)
    (VOL / f"{chave}.js").write_text(f"(window.VOL=window.VOL||{{}})[{json.dumps(chave)}]={json.dumps(pac)};")
    return chave


def volumes_m3(npz: Path, com_mep=True) -> list[str]:
    """HOMO, LUMO e (opcional) densidade eletrônica + MEP de uma molécula do M3."""
    from pyscf import dft
    from pyscf.tools import cubegen
    z = np.load(npz)
    mol = _mol(z)
    cc, coords = _grade(mol)
    ao = mol.eval_gto("GTOval", coords)
    C, occ = z["mo_coeff"], z["mo_occ"] > 0
    ih = int(np.where(occ)[0][-1])
    e = z["mo_energy"] * 27.211386
    feitos = []
    for nome, k in (("HOMO", ih), ("LUMO", ih + 1)):
        feitos.append(_pacote(f"{npz.stem}__{nome}", "orb", cc, mol, ao @ C[:, k],
                              {"rotulo": nome, "energia_ev": round(float(e[k]), 3),
                               "nivel": "CAM-B3LYP/6-31+G(d), PCM etanol"}))
    if com_mep:
        dm = (C[:, occ] * 2) @ C[:, occ].T
        rho = dft.numint.NumInt().eval_rho(mol, ao, dm)
        feitos.append(_pacote(f"{npz.stem}__dens", "dens", cc, mol, rho, {"rotulo": "densidade"}))
        tmp = VOL / f"_{npz.stem}_mep.cube"
        cubegen.mep(mol, str(tmp), dm, nx=N, ny=N, nz=N, margin=4.0)
        mep = _ler_cube(tmp); tmp.unlink()
        feitos.append(_pacote(f"{npz.stem}__MEP", "mep", cc, mol, mep,
                              {"rotulo": "MEP", "nivel": "CAM-B3LYP/6-31+G(d), PCM etanol"}))
    return feitos


def volume_spin(npz: Path) -> str:
    """Densidade de spin (α − β) de um radical do M2b."""
    from pyscf import dft
    z = np.load(npz)
    mol = _mol(z)
    cc, coords = _grade(mol)
    ao = mol.eval_gto("GTOval", coords)
    C, occ = z["mo_coeff"], z["mo_occ"]
    da = (C[0][:, occ[0] > 0]) @ C[0][:, occ[0] > 0].T
    db = (C[1][:, occ[1] > 0]) @ C[1][:, occ[1] > 0].T
    ni = dft.numint.NumInt()
    spin = ni.eval_rho(mol, ao, da) - ni.eval_rho(mol, ao, db)
    # populações de spin de Löwdin por átomo (menos sensíveis às funções difusas que as de Mulliken)
    from scipy.linalg import sqrtm
    S = mol.intor("int1e_ovlp")
    Sh = np.real(sqrtm(S))
    ps = np.diag(Sh @ (da - db) @ Sh)
    pop = [float(ps[p0:p1].sum()) for _, _, p0, p1 in mol.aoslice_by_atom()]
    return _pacote(f"{npz.stem}__spin", "spin", cc, mol, spin,
                   {"rotulo": "densidade de spin", "nivel": "UKS M06-2X/6-311+G(d,p), PCM metanol",
                    "spin_lowdin": [round(x, 3) for x in pop]})


def _ler_cube(arq: Path) -> np.ndarray:
    linhas = arq.read_text().split("\n")
    nat = abs(int(linhas[2].split()[0]))
    return np.array(" ".join(linhas[6 + nat:]).split(), float)


def gerar_volumes(forcar=False) -> dict:
    """Gera os volumes que ainda não existem. Retorna o manifesto {chave: metadados leves}."""
    VOL.mkdir(parents=True, exist_ok=True)
    for npz in sorted((M3 / "resultados" / "orbitais").glob("*.npz")):
        if forcar or not (VOL / f"{npz.stem}__HOMO.js").exists():
            volumes_m3(npz)
    for npz in sorted((ESP / "orbitais").glob("*.npz")):
        if forcar or not (VOL / f"{npz.stem}__spin.js").exists():
            volume_spin(npz)
    return manifesto()


def manifesto() -> dict:
    man = {}
    for f in sorted(VOL.glob("*.js")):
        if f.name.startswith("_"):
            continue
        txt = f.read_text()
        meta = json.loads(txt[txt.index("=", txt.index("]")) + 1:].rstrip(";"))
        meta.pop("dados"); meta.pop("xyz_bohr"); meta.pop("z")
        man[f.stem] = {k: meta[k] for k in ("tipo", "rotulo") if k in meta} | (
            {"energia_ev": meta["energia_ev"]} if "energia_ev" in meta else {}) | (
            {"nivel": meta["nivel"]} if "nivel" in meta else {}) | (
            {"spin_lowdin": meta["spin_lowdin"]} if "spin_lowdin" in meta else {})
    return man


if __name__ == "__main__":
    print(json.dumps(gerar_volumes(), indent=1, ensure_ascii=False))
