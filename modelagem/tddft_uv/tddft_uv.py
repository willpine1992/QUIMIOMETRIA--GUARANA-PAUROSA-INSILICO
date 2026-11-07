"""M3 — TD-DFT: espectros UV simulados dos constituintes (estudo Edson).

Pergunta: quais constituintes dos extratos absorvem em UVB (290–320 nm) e UVA (320–400 nm), e como se comparam
à benzofenona-3? Comparação qualitativa com a absorção UV observada por Fonseca Júnior (2021): extratos EtOH e DCM
absorvem fortemente, o óleo essencial pouco.

Protocolo
  Geometria: GFN2-xTB/ALPB. Moléculas do M2 reaproveitam a geometria neutra já otimizada (ALPB metanol);
             as novas (benzofenona-3, óxido de linalol) são otimizadas aqui (ALPB etanol). O efeito do solvente
             na geometria entre metanol e etanol é desprezível para excitação vertical.
  Excitação: TD-CAM-B3LYP/6-31+G(d) com ajuste de densidade (RI-JK; desde 05/11/2026), IEF-PCM etanol (ε = 24,852; solvente do ensaio de Mansur), 15 singletos.
  Espectro:  alargamento gaussiano em energia (FWHM 0,40 eV); λmax; áreas UVB/UVA2/UVA1; λ crítico (Diffey).
  M7:        o-quinona da epicatequina e desidrodicatequina A (fim da fila) testam se a oxidação explica a mudança
             de cor das emulsões com extrato EtOH (absorção deslocada para o visível, 400–500 nm).
  Orbitais:  HOMO, LUMO, gap e descritores de Koopmans (η, μ, ω) no resumo; por estado, a transição MO principal e o
             peso do NTO dominante. resultados/orbitais/<id>.npz guarda MOs e vetores X para gerar cubos
             (HOMO, LUMO, NTOs) depois com orbitais_cubos.py, sem refazer o cálculo.
  Calibração: benzofenona-3 roda primeiro; o desvio de λmax frente ao experimental mede o erro sistemático.

Uso
  "<venv>/bin/python" tddft_uv.py --preparar     só geometrias (rápido)
  "<venv>/bin/python" tddft_uv.py                fila TD-DFT (1 molécula por vez; retomável)
  "<venv>/bin/python" tddft_uv.py --analise      espectros + resumo a partir dos JSONs prontos
"""
from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

for v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(v, "1")

import numpy as np
import pandas as pd

from quimiolib import qm

AQUI = Path(__file__).resolve().parent
EDSON = AQUI.parents[1]
XTB_M2 = AQUI.parent / "dft_antioxidantes" / "especies" / "xtb"
GEO = AQUI / "geometrias"
RES = AQUI / "resultados"
BIB = EDSON / "biblioteca" / "saida" / "biblioteca_edson.csv"

# id → (nome, papel, origem da geometria: "m2:<arquivo>" ou "lib:<id da biblioteca>")
MOLS = {
    "benzofenona_3": ("Benzofenona-3", "calibração (filtro UV comercial)", "lib:REF006"),
    "linalol": ("Linalol", "óleo essencial de pau-rosa (majoritário)", "m2:linalol"),
    "oxido_linalol": ("cis-Óxido de linalol (furanoide)", "óleo essencial de pau-rosa", "lib:AR013"),
    "cafeina": ("Cafeína", "extratos DCM/EtOH de guaraná", "m2:cafeina"),
    "teobromina": ("Teobromina", "extratos DCM/EtOH de guaraná", "m2:teobromina"),
    "teofilina": ("Teofilina", "extratos DCM/EtOH de guaraná", "m2:teofilina"),
    "epicatequina": ("(−)-Epicatequina", "extrato EtOH de guaraná (unidade das procianidinas)", "m2:epicatequina"),
    "quercetina": ("Quercetina", "padrão (formulação F6)", "m2:quercetina"),
    # M7 · mecanismo da mudança de cor (F4, F9, F10): produtos de oxidação da epicatequina/catequina
    "epicatequina_oquinona": ("o-Quinona da epicatequina", "M7 · produto de oxidação (cor)",
                              "smi:C1[C@H]([C@H](OC2=CC(=CC(=C21)O)O)C3=CC(=O)C(=O)C=C3)O"),
    # opcional (66 átomos, o mais caro da fila): PubChem CID 182270, dímero amarelo de oxidação da catequina
    "desidrodicatequina_a": ("Desidrodicatequina A", "M7 · produto de oxidação (cor) · opcional",
                             "smi:C1[C@@H]([C@H](OC2=C1C(=CC3=C2C4=CC(=O)[C@@]5(CC4(O3)[C@@H]6[C@@H](O5)CC7=C(C=C(C=C7O6)O)O)O)O)C8=CC(=C(C=C8)O)O)O"),
}
OPCIONAIS = ["desidrodicatequina_a"]   # caro (66 átomos): roda por último, depois do M2b (rodar_m2b.sh)
ORDEM = [m for m in MOLS if m not in OPCIONAIS]  # calibração primeiro, depois do menor para o maior

FAIXAS = {"UVC": (200, 290), "UVB": (290, 320), "UVA2": (320, 340), "UVA1": (340, 400), "VIS": (400, 500)}
FWHM_EV = 0.40


def geometria(mid: str) -> dict:
    GEO.mkdir(exist_ok=True)
    arq = GEO / f"{mid}.json"
    if arq.exists():
        return json.loads(arq.read_text())
    origem = MOLS[mid][2]
    if origem.startswith("m2:"):
        x = json.loads((XTB_M2 / f"{origem[3:]}.json").read_text())
        g = {"simbolos": x["simbolos"], "xyz": x["xyz"], "fonte": f"M2 xTB/ALPB(metanol) · n_imag={x.get('n_imag')}"}
    else:
        smi = origem[4:] if origem.startswith("smi:") else pd.read_csv(BIB).set_index("id").loc[origem[4:], "smiles"]
        _, simb, xyz = qm.conformero_xtb(smi, solvente="ethanol")
        xyz, _, ok = qm.xtb_otimizar(simb, xyz, solvente="ethanol")
        g = {"simbolos": simb, "xyz": np.asarray(xyz).tolist(), "smiles": smi, "opt_ok": bool(ok),
             "fonte": "conformero ETKDGv3+MMFF → GFN2-xTB/ALPB(etanol)"}
    arq.write_text(json.dumps(g))
    return g


def rodar(fila=None):
    RES.mkdir(exist_ok=True)
    fila = fila or ORDEM
    for i, mid in enumerate(fila, 1):
        arq = RES / f"{mid}.json"
        if arq.exists():
            continue
        g = geometria(mid)
        print(f"[{i}/{len(fila)}] {mid} ({len(g['simbolos'])} átomos) …", flush=True)
        t = time.time()
        (RES / "logs").mkdir(exist_ok=True)
        r = qm.tddft_uv(g["simbolos"], g["xyz"], nstates=15, salvar=RES / "orbitais" / f"{mid}.npz",
                        log=RES / "logs" / f"{mid}.log")   # ajuste de densidade automático (qm.tddft_uv)
        r.update({"id": mid, "nome": MOLS[mid][0], "papel": MOLS[mid][1], "geometria": g["fonte"],
                  "n_atomos": len(g["simbolos"])})
        arq.write_text(json.dumps(r, indent=1, ensure_ascii=False))
        print(f"  ok em {(time.time() - t) / 60:.1f} min · S1 = {r['lambda_nm'][0]:.0f} nm", flush=True)


def espectro(e_ev, f, grade_nm):
    """ε relativo (soma de gaussianas em energia, amplitude ∝ f)."""
    s = FWHM_EV / (2 * np.sqrt(2 * np.log(2)))
    E = 1239.84193 / grade_nm
    return sum(fi * np.exp(-((E - ei) ** 2) / (2 * s ** 2)) for ei, fi in zip(e_ev, f))


def analisar():
    grade = np.arange(200, 601, 1.0)
    linhas, curvas = [], {"lambda_nm": grade}
    for mid in list(MOLS):
        arq = RES / f"{mid}.json"
        if not arq.exists():
            continue
        r = json.loads(arq.read_text())
        y = espectro(r["e_ev"], r["f"], grade)
        curvas[mid] = np.round(y, 5)
        d = {"id": mid, "nome": r["nome"], "papel": r["papel"], "n_atomos": r["n_atomos"],
             "convergiu": r["scf_convergiu"] and r["td_convergiu"], "t_min": round(r["t_total_s"] / 60, 1)}
        m = grade >= 250
        d["lambda_max_nm"] = float(grade[m][np.argmax(y[m])]) if y[m].max() > 1e-4 else None
        # estado mais intenso acima de 250 nm
        acima = [(l, fi) for l, fi in zip(r["lambda_nm"], r["f"]) if l >= 250]
        d["estado_forte_nm"], d["estado_forte_f"] = max(acima, key=lambda t: t[1]) if acima else (None, None)
        if "homo_ev" in r:  # orbitais de fronteira (CAM-B3LYP) e descritores globais de Koopmans
            h, l = r["homo_ev"], r["lumo_ev"]
            eta, mu = (l - h) / 2, (h + l) / 2
            d.update({"homo_ev": round(h, 3), "lumo_ev": round(l, 3), "gap_ev": round(l - h, 3),
                      "dureza_eta_ev": round(eta, 3), "potencial_quimico_mu_ev": round(mu, 3),
                      "eletrofilicidade_omega_ev": round(mu ** 2 / (2 * eta), 3)})
            if acima:
                k = max((i for i, lam in enumerate(r["lambda_nm"]) if lam >= 250), key=lambda i: r["f"][i])
                d["estado_forte_carater"] = r["transicao_principal"][k]
                d["estado_forte_nto_peso"] = r["nto_peso"][k]
        for k, (a, b) in FAIXAS.items():
            mk = (grade >= a) & (grade <= b)
            d[f"area_{k}"] = round(float(np.trapezoid(y[mk], grade[mk])), 3)
        mk = (grade >= 290) & (grade <= 400)
        acum = np.cumsum(y[mk]); tot = acum[-1]
        d["lambda_critico_nm"] = float(grade[mk][np.searchsorted(acum, 0.9 * tot)]) if tot > 1e-3 else None
        linhas.append(d)
    if not linhas:
        print("nenhum resultado ainda"); return
    res = pd.DataFrame(linhas)
    res.to_csv(AQUI / "resumo_uv.csv", index=False)
    pd.DataFrame(curvas).to_csv(AQUI / "espectros_uv.csv", index=False)
    print(res.to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--preparar", action="store_true")
    ap.add_argument("--analise", action="store_true")
    ap.add_argument("--opcionais", action="store_true", help="roda só as moléculas opcionais (dímero)")
    a = ap.parse_args()
    if a.preparar:
        for mid in MOLS:
            g = geometria(mid)
            print(f"{mid:15s} {len(g['simbolos']):3d} átomos · {g['fonte']}")
    elif a.analise:
        analisar()
    elif a.opcionais:
        rodar(OPCIONAIS)
        analisar()
    else:
        rodar()
        analisar()
