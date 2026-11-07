"""DFT dos antioxidantes — mecanismos HAT, SET-PT e SPLET (Edson: guaraná × pau-rosa × padrões).

Pergunta: por que, na tese, só o extrato ETANÓLICO de guaraná sequestra DPPH, enquanto DCM, hexânico e o
óleo de pau-rosa não? (uso apenas QUALITATIVO do resultado experimental; nenhum dado da tese é entrada).

Protocolo (triagem, "M06-2X/6-311+G(d,p)//GFN2-xTB", metanol = solvente do ensaio DPPH;
  SCF direto, sem density fitting; PCM com Lebedev 17 = 110 pontos/átomo):
  1. conformero: ETKDGv3 (20) + MMFF → menor energia GFN2-xTB/ALPB(MeOH)
  2. geometria + frequências (RRHO 298.15 K): GFN2-xTB/ALPB(MeOH) para neutro, radicais, ânions, cátion
  3. energia: single-point M06-2X/6-311+G(d,p), DF, IEF-PCM (ε = 32.613, metanol)
  4. entalpias de reação (kcal/mol):
       HAT    BDE = H(R•)  + H(H•)  − H(RH)
       SET-PT IP  = H(RH•+) + H(e−) − H(RH)        PDE = H(R•) + H(H+) − H(RH•+)
       SPLET  PA  = H(R−)  + H(H+)  − H(RH)        ETE = H(R•) + H(e−) − H(R−)
     H(H•): calculado no mesmo nível (+ 5/2 RT).
     H(H+) = 5/2RT(gás, 6.197 kJ/mol) + ΔsolvH(H+, MeOH) = −1038 kJ/mol
             [IEF-PCM B3LYP/aug-cc-pVTZ; Rottmannová et al. 2013, Acta Chim Slov 6:60, doi:10.2478/acs-2013-0011]
     H(e−) = 3.145 kJ/mol (gás, Bartmess 1994) + ΔsolvH(e−, MeOH) = −86 kJ/mol
             [B3LYP/6-311++G**, Rimarčík et al. 2010, conforme citado na literatura posterior]
     IP/PDE/PA/ETE absolutos dependem dessas constantes → também são reportados Δ em relação ao fenol.

Modelos: procianidinas (B1–B4, C1, A2) representadas pela unidade (−)-epicatequina (mesmo anel B catecol);
(+)-catequina é diastereoisômero da epicatequina (energias praticamente iguais na literatura) → não repetida.

Rodar (demora horas; retomável — cada espécie vira um JSON em ./especies/):
  "QUIMIOINFORMÁTICA/.venv/bin/python" dft_antioxidantes.py [--workers 6]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

from quimiolib import qm

AQUI = Path(__file__).resolve().parent
ESP = AQUI / "especies"
ESP.mkdir(exist_ok=True)

RT = qm.KB_KCAL * 298.15
KJ = 1 / 4.184
H_PROTON = (6.197 - 1038.0) * KJ       # kcal/mol
H_ELETRON = (3.145 - 86.0) * KJ        # kcal/mol

# nome: (SMILES com mapas nos sítios, {mapa: rótulo}, tipo do sítio por mapa, origem na tese)
MOLS = {
    "fenol": ("c1ccc([OH:1])cc1", {1: "O–H"}, "calibração"),
    "epicatequina": ("C1[C@H]([C@H](OC2=CC(=CC(=C21)[OH:5])[OH:7])C3=CC(=C(C=C3)[OH:14])[OH:13])[OH:3]",
                     {3: "3-OH", 5: "5-OH", 7: "7-OH", 13: "3'-OH", 14: "4'-OH"}, "guaraná · EtOH"),
    "cafeina": ("Cn1c(=O)c2c(n[cH:8]n2C)n(C)c1=O", {8: "C8–H"}, "guaraná · DCM/EtOH"),
    "teobromina": ("Cn1[cH:8]nc2c1c(=O)[nH:1]c(=O)n2C", {1: "N1–H", 8: "C8–H"}, "guaraná · DCM/EtOH"),
    "teofilina": ("Cn1c(=O)c2[nH:7][cH:8]nc2n(C)c1=O", {7: "N7–H", 8: "C8–H"}, "guaraná · DCM/EtOH"),
    "linalol": ("CC(C)=C[CH2:4]CC(C)([OH:1])C=C", {1: "O–H", 4: "C4–H (alílico)"}, "pau-rosa · OE"),
    "quercetina": ("O=C1C([OH:3])=C(Oc2cc([OH:7])cc([OH:5])c12)c1ccc([OH:14])c([OH:13])c1",
                   {3: "3-OH", 5: "5-OH", 7: "7-OH", 13: "3'-OH", 14: "4'-OH"}, "padrão (F6)"),
    "acido_galico": ("OC(=O)c1cc([OH:3])c([OH:4])c(O)c1", {3: "3-OH (=5-OH)", 4: "4-OH"}, "padrão Folin"),
}


# ------------------------------------------------------------------ jobs (top-level p/ multiprocessing)
def job_neutro(nome, smiles):
    m, simb, xyz = qm.conformero_xtb(smiles)
    xyz, e, ok = qm.xtb_otimizar(simb, xyz)
    termo = qm.xtb_termo(simb, xyz)
    sitios = {}
    for a in m.GetAtoms():
        if a.GetAtomMapNum():
            hs = [n.GetIdx() for n in a.GetNeighbors() if n.GetAtomicNum() == 1]
            sitios[a.GetAtomMapNum()] = {"atomo": a.GetIdx(), "elemento": a.GetSymbol(), "hs": hs}
    return {"simbolos": simb, "xyz": xyz.tolist(), "charge": 0, "spin": 0, "e_xtb": e,
            "opt_ok": ok, **termo, "sitios": {str(k): v for k, v in sitios.items()}}


def job_derivado(nome, simbolos, xyz, remover, charge, spin):
    simb = [s for i, s in enumerate(simbolos) if i != remover]
    xyz = np.array([c for i, c in enumerate(xyz) if i != remover])
    xyz, e, ok = qm.xtb_otimizar(simb, xyz, charge=charge, uhf=spin)
    termo = qm.xtb_termo(simb, xyz, charge=charge, uhf=spin)
    return {"simbolos": simb, "xyz": xyz.tolist(), "charge": charge, "spin": spin,
            "e_xtb": e, "opt_ok": ok, **termo}


def job_dft(nome, simbolos, xyz, charge, spin):
    return qm.dft_sp(simbolos, xyz, charge=charge, spin=spin, df=False, max_memory=500)  # 500 MB força SCF direto (sem ERI em RAM)


# ------------------------------------------------------------------ pipeline
def main(workers: int):
    xtb_dir = ESP / "xtb"; xtb_dir.mkdir(exist_ok=True)
    dft_dir = ESP / "dft"; dft_dir.mkdir(exist_ok=True)

    print("1/3 neutros (conformeros + xTB)…", flush=True)
    neutros = qm.rodar_paralelo(job_neutro, [{"nome": n, "smiles": s} for n, (s, _, _) in MOLS.items()],
                                workers, xtb_dir)

    print("2/3 radicais, ânions e cátions (xTB)…", flush=True)
    tarefas = []
    for n, (smi, rot, _) in MOLS.items():
        nt = neutros[n]
        tarefas.append({"nome": f"{n}__cation", "simbolos": nt["simbolos"], "xyz": nt["xyz"],
                        "remover": -1, "charge": 1, "spin": 1})
        for mp, s in nt["sitios"].items():
            h = s["hs"][0]
            tarefas.append({"nome": f"{n}__rad_{mp}", "simbolos": nt["simbolos"], "xyz": nt["xyz"],
                            "remover": h, "charge": 0, "spin": 1})
            if s["elemento"] in ("O", "N"):   # C–H não é ácido: sem ânion
                tarefas.append({"nome": f"{n}__an_{mp}", "simbolos": nt["simbolos"], "xyz": nt["xyz"],
                                "remover": h, "charge": -1, "spin": 0})
    derivados = qm.rodar_paralelo(job_derivado, tarefas, workers, xtb_dir)

    print("3/3 DFT single points…", flush=True)
    todas = {**{n: neutros[n] for n in MOLS}, **derivados}
    dft_t = [{"nome": k, "simbolos": v["simbolos"], "xyz": v["xyz"], "charge": v["charge"], "spin": v["spin"]}
             for k, v in todas.items()]
    dft_t.append({"nome": "H_atomo", "simbolos": ["H"], "xyz": [[0.0, 0.0, 0.0]], "charge": 0, "spin": 1})
    # fila por molécula (menores primeiro): cada molécula fecha completa e já dá resultado parcial
    ordem = ["H_atomo", "fenol", "cafeina", "teobromina", "teofilina", "acido_galico", "linalol",
             "epicatequina", "quercetina"]
    dft_t.sort(key=lambda t: (ordem.index(t["nome"].split("__")[0]), t["nome"]))
    dfts = qm.rodar_paralelo(job_dft, dft_t, workers, dft_dir)

    analisar(todas, dfts)


def carregar():
    """Lê os JSONs já calculados: (xtb por espécie, dft por espécie)."""
    todas = {f.stem: json.loads(f.read_text()) for f in (ESP / "xtb").glob("*.json")}
    dfts = {f.stem: json.loads(f.read_text()) for f in (ESP / "dft").glob("*.json")}
    return todas, dfts


def calcular(todas, dfts):
    """Descritores por sítio (DataFrame, kcal/mol) + lista de avisos. Só moléculas com neutro+cátion+radical."""
    if "H_atomo" not in dfts:  # calcula na hora (< 1 s)
        dfts["H_atomo"] = job_dft("H_atomo", ["H"], [[0.0, 0.0, 0.0]], 0, 1)
        (ESP / "dft" / "H_atomo.json").write_text(json.dumps(dfts["H_atomo"], indent=1))
    H = lambda k: dfts[k]["e_dft"] * qm.HARTREE_KCAL + todas[k]["h_corr_kcal"]
    H_H = dfts["H_atomo"]["e_dft"] * qm.HARTREE_KCAL + 2.5 * RT

    linhas, avisos = [], []
    for n, (_, rot, origem) in MOLS.items():
        if n not in dfts or f"{n}__cation" not in dfts:
            continue
        ip = H(f"{n}__cation") + H_ELETRON - H(n)
        for mp, lab in rot.items():
            r, a = f"{n}__rad_{mp}", f"{n}__an_{mp}"
            if r not in dfts:
                continue
            d = {"molecula": n, "origem_tese": origem, "sitio": lab,
                 "BDE": H(r) + H_H - H(n), "IP": ip, "PDE": H(r) + H_PROTON - H(f"{n}__cation")}
            if a in dfts:
                d["PA"] = H(a) + H_PROTON - H(n)
                d["ETE"] = H(r) + H_ELETRON - H(a)
            d["S2_radical"] = dfts[r]["s2"]
            linhas.append(d)
        for k in [n, f"{n}__cation"] + [f"{n}__rad_{m}" for m in rot] + [f"{n}__an_{m}" for m in rot]:
            if k in todas and todas[k].get("n_imag", 0):
                avisos.append(f"{k}: {todas[k]['n_imag']} freq. imaginária(s) (min {todas[k]['freq_min']:.0f} cm-1)")
            if k in dfts and not dfts[k].get("convergiu", True):
                avisos.append(f"{k}: SCF não convergiu")

    df = pd.DataFrame(linhas)
    if df.empty:
        return df, avisos
    for c in ("PA", "ETE"):
        if c not in df:
            df[c] = np.nan
    if "fenol" in set(df.molecula):
        ref = df[df.molecula == "fenol"].iloc[0]
        for c in ("BDE", "IP", "PA"):
            df[f"Δ{c}_vs_fenol"] = df[c] - ref[c]
    return df.round(2), avisos


def analisar(todas=None, dfts=None):
    if todas is None:
        todas, dfts = carregar()
    print(f"DFT pronto para {len(dfts)} espécies (de {len(todas) + 1}).")
    df, avisos = calcular(todas, dfts)
    if df.empty:
        print("Nenhuma molécula completa ainda (precisa de neutro + cátion + ao menos um radical).")
        return
    faltando = [n for n in MOLS if n not in set(df.molecula)]
    if faltando:
        print("Incompletas (ainda sem DFT suficiente):", ", ".join(faltando))
    df.to_csv(AQUI / "descritores_por_sitio.csv", index=False)

    resumo = df.groupby(["molecula", "origem_tese"]).agg(
        BDE_min=("BDE", "min"), IP=("IP", "first"), PA_min=("PA", "min")).reset_index()
    resumo["sitio_BDE_min"] = resumo["molecula"].map(df.loc[df.groupby("molecula")["BDE"].idxmin()]
                                                     .set_index("molecula")["sitio"])
    resumo = resumo.sort_values("BDE_min").round(1)
    resumo.to_csv(AQUI / "resumo_por_molecula.csv", index=False)
    (AQUI / "avisos.txt").write_text("\n".join(avisos) or "nenhum")
    print(resumo.to_string(index=False))
    print("\navisos:", avisos or "nenhum")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--so-analise", action="store_true")
    a = ap.parse_args()
    analisar() if a.so_analise else main(a.workers)
