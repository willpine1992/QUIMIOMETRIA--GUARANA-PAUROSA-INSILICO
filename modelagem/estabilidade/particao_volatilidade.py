"""M7a — Onde cada constituinte fica na emulsão O/A e quanto ele tende a escapar para o ar (estudo Edson).

Pergunta: as alterações vistas no estudo de estabilidade preliminar (FONSECA JÚNIOR, 2021) — perda do odor do OE de
pau-rosa (F2, F8, F9, F10) e mudança de cor com o extrato etanólico (F4, F9, F10) — são compatíveis com a partição
das moléculas entre gota de óleo, água e ar?

Protocolo (solvatação implícita, qualitativo)
  Para cada fase φ ∈ {gás, água, 1-octanol, n-hexadecano}: 20 confôrmeros ETKDGv3+MMFF → menor energia GFN2-xTB/ALPB(φ)
  → otimização GFN2-xTB/ALPB(φ). ΔG_solv(φ) = E(φ) − E(gás), cada um na sua geometria relaxada.
  log P(oct/água)  = [ΔG(água) − ΔG(oct)] / (2,303 RT)      → validação contra log Kow experimental (PubChem)
  log K(óleo/água) = [ΔG(água) − ΔG(hexadecano)] / (2,303 RT) → n-hexadecano como modelo do óleo mineral/IPM
  log K(ar/água)   = ΔG(água) / (2,303 RT)                    → constante de Henry adimensional (Cg/Cw)
  log K(ar/óleo)   = ΔG(hexadecano) / (2,303 RT)              → quanto escapa da gota de óleo para o espaço livre
  Estado-padrão 1 mol/L nas duas fases (convenção "gsolv" do ALPB), 298,15 K.
  pH das emulsões 3,2–4,2: os fenóis (pKa > 7) e as xantinas (forma protonada só abaixo de pH ~1) estão neutros,
  então log D ≈ log P nessa faixa.

Uso
  "<venv>/bin/python" particao_volatilidade.py [--workers 3]   → particao_volatilidade.csv
"""
from __future__ import annotations

import argparse
import json
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

for v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(v, "1")

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen

from quimiolib import qm

AQUI = Path(__file__).resolve().parent
EDSON = AQUI.parents[1]
BIB = EDSON / "biblioteca" / "saida" / "biblioteca_edson.csv"
CACHE = AQUI / "xtb_fases"
RT = qm.KB_KCAL * 298.15
FASES = {"gas": None, "agua": "water", "octanol": "octanol", "hexadecano": "hexadecane"}

# id → (nome, origem na formulação). "lib:<id>" vem da biblioteca; "smi:<SMILES>" é estrutura nova.
MOLS = {
    "linalol": ("Linalol", "OE pau-rosa (majoritário)", "lib:AR016"),
    "oxido_linalol_f": ("cis-Óxido de linalol (furanoide)", "OE pau-rosa", "lib:AR013"),
    "alfa_terpineol": ("α-Terpineol", "OE pau-rosa", "lib:AR021"),
    "geraniol": ("Geraniol", "OE pau-rosa", "lib:AR015"),
    "limoneno": ("Limoneno", "OE pau-rosa", "lib:AR008"),
    "benzoato_benzila": ("Benzoato de benzila", "OE pau-rosa", "lib:AR002"),
    "linalol_6_ooh": ("Hidroperóxido de linalol (6-OOH)", "produto de oxidação do linalol", "lib:AR004"),
    "linalol_7_ooh": ("Hidroperóxido de linalol (7-OOH)", "produto de oxidação do linalol", "lib:AR005"),
    "cafeina": ("Cafeína", "extratos DCM/EtOH de guaraná", "lib:PC004"),
    "teobromina": ("Teobromina", "extratos DCM/EtOH de guaraná", "lib:PC005"),
    "teofilina": ("Teofilina", "extratos DCM/EtOH de guaraná", "lib:PC006"),
    "epicatequina": ("(−)-Epicatequina", "extrato EtOH de guaraná", "lib:PC003"),
    "catequina": ("(+)-Catequina", "extrato EtOH de guaraná", "lib:PC001"),
    "epicatequina_oquinona": ("o-Quinona da epicatequina", "produto de oxidação da epicatequina",
                              "smi:C1[C@H]([C@H](OC2=CC(=CC(=C21)O)O)C3=CC(=O)C(=O)C=C3)O"),
    "quercetina": ("Quercetina", "padrão (F6)", "lib:REF004"),
}

# log Kow e Henry experimentais (PubChem, seções "LogP" e "Henry's Law Constant", consultadas em 04/11/2026).
# Henry convertido de atm·m³/mol para adimensional: K_aw = H / (R T), R T = 0,024465 atm·m³/mol.
EXP_LOGP = {"cafeina": -0.07, "teobromina": -0.78, "teofilina": -0.02, "linalol": 2.97, "quercetina": 2.13,
            "alfa_terpineol": 2.98, "geraniol": 3.56, "benzoato_benzila": 3.97}
EXP_LOGKAW = {"linalol": round(math.log10(2.15e-5 / 0.024465), 2),
              "alfa_terpineol": round(math.log10(2.30e-6 / 0.024465), 2)}


def smiles(mid):
    o = MOLS[mid][2]
    if o.startswith("smi:"):
        return o[4:]
    return pd.read_csv(BIB).set_index("id").loc[o[4:], "smiles"]


def fase(mid: str, nome_fase: str) -> dict:
    """Energia GFN2-xTB relaxada numa fase (cache em JSON)."""
    CACHE.mkdir(exist_ok=True)
    arq = CACHE / f"{mid}__{nome_fase}.json"
    if arq.exists():
        return json.loads(arq.read_text())
    solv = FASES[nome_fase]
    _, simb, xyz = qm.conformero_xtb(smiles(mid), solvente=solv)
    xyz, e, ok = qm.xtb_otimizar(simb, xyz, solvente=solv)
    r = {"mid": mid, "fase": nome_fase, "e_h": e, "convergiu": bool(ok), "simbolos": simb, "xyz": xyz.tolist()}
    arq.write_text(json.dumps(r))
    return r


def tarefa(mid, nome_fase):
    r = fase(mid, nome_fase)
    return mid, nome_fase, r["e_h"], r["convergiu"]


def main(workers: int):
    tarefas = [(m, f) for m in MOLS for f in FASES]
    E, conv = {}, {}
    with ProcessPoolExecutor(workers) as ex:
        futs = [ex.submit(tarefa, m, f) for m, f in tarefas]
        for fu in as_completed(futs):
            m, f, e, ok = fu.result()
            E[(m, f)], conv[(m, f)] = e, ok
            print(f"  {m:24s} {f:10s} ok={ok}", flush=True)
    L = 2.302585 * RT
    linhas = []
    for m, (nome, origem, _) in MOLS.items():
        dG = {f: (E[(m, f)] - E[(m, "gas")]) * qm.HARTREE_KCAL for f in FASES if f != "gas"}
        mol = Chem.MolFromSmiles(smiles(m))
        d = {"id": m, "nome": nome, "origem": origem,
             "dGsolv_agua": round(dG["agua"], 2), "dGsolv_octanol": round(dG["octanol"], 2),
             "dGsolv_hexadecano": round(dG["hexadecano"], 2),
             "logP_xtb": round((dG["agua"] - dG["octanol"]) / L, 2),
             "logP_crippen": round(Crippen.MolLogP(mol), 2),
             "logP_exp": EXP_LOGP.get(m),
             "logK_oleo_agua": round((dG["agua"] - dG["hexadecano"]) / L, 2),
             "logK_ar_agua": round(dG["agua"] / L, 2),
             "logK_ar_agua_exp": EXP_LOGKAW.get(m),
             "logK_ar_oleo": round(dG["hexadecano"] / L, 2),
             "todas_convergiram": all(conv[(m, f)] for f in FASES)}
        # leitura qualitativa da localização na emulsão O/A (heurística; ver proposta, M7)
        k = d["logK_oleo_agua"]
        d["local_provavel"] = "gota de óleo" if k > 1 else ("fase aquosa" if k < -1 else "interface / as duas fases")
        linhas.append(d)
    df = pd.DataFrame(linhas)
    # correção sistemática do xTB: offset médio frente ao experimental (só para leitura; colunas brutas mantidas)
    val = df.dropna(subset=["logP_exp"])
    if len(val) >= 3:
        err = val["logP_xtb"] - val["logP_exp"]
        print(f"\nlogP xTB vs exp (n={len(val)}): MAE = {err.abs().mean():.2f} · viés = {err.mean():+.2f}")
        cr = val["logP_crippen"] - val["logP_exp"]
        print(f"logP Crippen vs exp (n={len(val)}): MAE = {cr.abs().mean():.2f} · viés = {cr.mean():+.2f}")
    df.to_csv(AQUI / "particao_volatilidade.csv", index=False)
    print(df[["nome", "logP_xtb", "logP_exp", "logK_oleo_agua", "logK_ar_agua", "logK_ar_agua_exp",
              "logK_ar_oleo", "local_provavel"]].to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=3)
    main(ap.parse_args().workers)
