"""M2b — refino das geometrias com frequência imaginária + funções de onda dos radicais-chave (estudo Edson).

Etapa A · reotimização (pendência de qualidade do M2)
  Para cada espécie com n_imag > 0 no xTB: desloca a geometria ±0,15 Å (amplitude máxima) ao longo do modo
  imaginário, reotimiza (GFN2-xTB/ALPB metanol) as duas, fica com a de menor energia e recalcula as frequências;
  repete até 3 vezes. Depois refaz o single point DFT (mesmo nível do M2). Os JSONs antigos vão para
  especies/_antes_reotimizacao/. As BDE/PA/IP da cafeína e da teobromina são então recalculadas pelo próprio M2.

Etapa B · funções de onda para as figuras
  Single point UKS M06-2X/6-311+G(d,p), IEF-PCM metanol (idêntico ao M2), gravando os orbitais em
  especies/orbitais/<espécie>.npz: radical de menor BDE de fenol, epicatequina (3'-O• e 4'-O•), quercetina,
  linalol e cafeína. Daí saem a densidade de spin (densidades.py), o MEP e o NCI.

  "<venv>/bin/python" m2b_refino_densidades.py [--so A|B]
"""
from __future__ import annotations

import argparse
import json
import shutil
import time

import numpy as np

from dft_antioxidantes import ESP, calcular, carregar, job_dft
from quimiolib import qm

XTB, DFT = ESP / "xtb", ESP / "dft"
ANTES = ESP / "_antes_reotimizacao"
ORB = ESP / "orbitais"
LOG = ESP.parent / "m2b.log"


def log(msg):
    linha = f"[{time.strftime('%d/%m %H:%M')}] {msg}"
    print(linha, flush=True)
    with open(LOG, "a") as f:
        f.write(linha + "\n")


def reotimizar(nome: str, x: dict, amp=0.15, max_iter=3) -> dict:
    simb, xyz, ch, uhf = x["simbolos"], np.asarray(x["xyz"]), x["charge"], x["spin"]
    termo = qm.xtb_termo(simb, xyz, charge=ch, uhf=uhf, modo=True)
    hist = [round(termo["freq_min"], 1)]
    for _ in range(max_iter):
        if termo["n_imag"] == 0:
            break
        modo = np.asarray(termo["modo_imag"])
        cand = []
        for sinal in (1, -1):
            g, e, ok = qm.xtb_otimizar(simb, xyz + sinal * amp * modo, charge=ch, uhf=uhf, maxiter=500)
            cand.append((e, g, ok))
        e, xyz, ok = min(cand, key=lambda c: c[0])
        termo = qm.xtb_termo(simb, xyz, charge=ch, uhf=uhf, modo=True)
        hist.append(round(termo["freq_min"], 1))
    termo.pop("modo_imag", None)
    novo = {**x, "xyz": np.asarray(xyz).tolist(), "e_xtb": e if len(hist) > 1 else x["e_xtb"], "opt_ok": True,
            **termo, "reotimizado": True, "freq_min_historico": hist}
    return novo


def etapa_a():
    todas, _ = carregar()
    alvo = [k for k, v in sorted(todas.items()) if v.get("n_imag", 0) > 0]
    log(f"Etapa A: {len(alvo)} espécie(s) com freq. imaginária: {', '.join(alvo) or '—'}")
    ANTES.mkdir(exist_ok=True)
    for k in alvo:
        for pasta in (XTB, DFT):
            f = pasta / f"{k}.json"
            if f.exists() and not (ANTES / f"{pasta.name}__{k}.json").exists():
                shutil.copy2(f, ANTES / f"{pasta.name}__{k}.json")
        t = time.time()
        novo = reotimizar(k, todas[k])
        (XTB / f"{k}.json").write_text(json.dumps(novo))
        log(f"  {k}: freq. mín. {novo['freq_min_historico']} cm-1 → n_imag={novo['n_imag']} "
            f"({(time.time() - t) / 60:.1f} min xTB)")
        t = time.time()
        d = job_dft(k, novo["simbolos"], novo["xyz"], novo["charge"], novo["spin"])
        (DFT / f"{k}.json").write_text(json.dumps(d))
        log(f"  {k}: DFT refeito em {(time.time() - t) / 60:.1f} min (convergiu={d['convergiu']})")


def radicais_chave() -> list[str]:
    todas, dfts = carregar()
    df, _ = calcular(todas, dfts)
    escolha = []
    for mol in ("fenol", "epicatequina", "quercetina", "linalol", "cafeina"):
        sub = df[df.molecula == mol]
        if sub.empty:
            continue
        sitio = sub.loc[sub.BDE.idxmin(), "sitio"]
        mp = next(k for k, v in __import__("dft_antioxidantes").MOLS[mol][1].items() if v == sitio)
        escolha.append(f"{mol}__rad_{mp}")
    for extra in ("epicatequina__rad_13", "epicatequina__rad_14"):   # os dois O–H do catecol
        if extra not in escolha:
            escolha.append(extra)
    return escolha


def etapa_b():
    todas, _ = carregar()
    alvo = radicais_chave()
    log(f"Etapa B: funções de onda de {', '.join(alvo)}")
    for k in alvo:
        npz = ORB / f"{k}.npz"
        if npz.exists():
            continue
        x = todas[k]
        t = time.time()
        d = qm.dft_sp(x["simbolos"], x["xyz"], charge=x["charge"], spin=x["spin"], df=False, max_memory=500,
                      salvar=npz)
        log(f"  {k}: {(time.time() - t) / 60:.1f} min · ⟨S²⟩={d['s2']:.3f} · convergiu={d['convergiu']}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--so", choices=["A", "B"])
    a = ap.parse_args()
    if a.so in (None, "A"):
        etapa_a()
    if a.so in (None, "B"):
        etapa_b()
    log("M2b concluído.")
