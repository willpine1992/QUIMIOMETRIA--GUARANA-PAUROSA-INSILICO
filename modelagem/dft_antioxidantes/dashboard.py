"""Gera o dashboard HTML autocontido do DFT de antioxidantes (Edson).

Lê os JSONs já calculados (parciais ou completos), calcula os descritores com a MESMA função do pipeline
(`dft_antioxidantes.calcular`) e grava um único .html com dados, CSS e JS embutidos (sem CDN, abre offline).

  "<venv>/bin/python" dashboard.py            → edson/dashboard_dft_antioxidantes.html
O layout (storytelling: enigma → capítulos com questão central) fica em template_dashboard.html.
"""
from __future__ import annotations

import datetime as dt
import json
import math
from pathlib import Path

from dft_antioxidantes import AQUI, ESP, H_ELETRON, H_PROTON, MOLS, calcular, carregar

SAIDA = AQUI.parents[1] / "dashboard_dft_antioxidantes.html"

ORDEM = ["H_atomo", "fenol", "cafeina", "teobromina", "teofilina", "acido_galico", "linalol",
         "epicatequina", "quercetina"]
META = {
    "fenol": ("Fenol", "calib", "calibração", "—"),
    "epicatequina": ("(−)-Epicatequina", "flavanol", "extrato EtOH de guaraná",
                     "extrato EtOH: ATIVO no DPPH"),
    "cafeina": ("Cafeína", "xantina", "extratos DCM/EtOH de guaraná", "extrato DCM: inativo no DPPH"),
    "teobromina": ("Teobromina", "xantina", "extratos DCM/EtOH de guaraná", "extrato DCM: inativo no DPPH"),
    "teofilina": ("Teofilina", "xantina", "extratos DCM/EtOH de guaraná", "extrato DCM: inativo no DPPH"),
    "linalol": ("Linalol", "oe", "óleo essencial de pau-rosa", "óleo essencial: inativo no DPPH"),
    "quercetina": ("Quercetina", "padrao", "padrão antioxidante (formulação F6)", "padrão: ATIVO no DPPH"),
    "acido_galico": ("Ácido gálico", "padrao", "padrão de Folin-Ciocalteu", "não testado no DPPH"),
}
GRUPOS = {
    "flavanol": "Guaraná · flavan-3-ol",
    "xantina": "Guaraná · metilxantina",
    "oe": "Pau-rosa · óleo essencial",
    "padrao": "Padrão de referência",
    "calib": "Calibração",
}
TIPO = {"neutro": "neutro", "cation": "cátion-radical", "rad": "radical", "an": "ânion"}


def fila(todas, dfts):
    nomes = sorted(list(todas) + ["H_atomo"], key=lambda n: (ORDEM.index(n.split("__")[0]), n))
    linhas = []
    for i, n in enumerate(nomes, 1):
        mol = n.split("__")[0]
        suf = n.split("__")[1] if "__" in n else "neutro"
        t = suf.split("_")[0]
        mp = suf.split("_")[1] if "_" in suf else None
        sitio = MOLS[mol][1].get(int(mp), "") if (mp and mol in MOLS) else ""
        if n == "H_atomo":
            nat, ch, mult, imag = 1, 0, 2, None
        else:
            x = todas[n]
            nat, ch, mult, imag = len(x["simbolos"]), x["charge"], x["spin"] + 1, x["n_imag"]
        d = dfts.get(n)
        linhas.append({
            "ordem": i, "id": n, "molecula": META.get(mol, ("Átomo de H",))[0] if mol != "H_atomo" else "Átomo de H",
            "mol": mol, "especie": "átomo H" if n == "H_atomo" else TIPO[t], "sitio": sitio or "—",
            "carga": ch, "mult": mult, "atomos": nat, "imag": imag,
            "pronto": d is not None, "t_min": round(d["t_dft_s"] / 60, 1) if d else None,
            "s2": round(d["s2"], 3) if d and mult == 2 else None, "conv": d["convergiu"] if d else None,
        })
    return linhas


def estimativa(q):
    """Minutos estimados por espécie pendente {id: min}. Ajusta t = k·N^a nas espécies prontas (a livre se houver
    ≥2 tamanhos; senão a=3), com fator 1,4 para camada aberta (UKS). None se ainda não há base."""
    pts = [(r["atomos"], r["t_min"] / (1.4 if r["mult"] == 2 else 1.0))
           for r in q if r["pronto"] and r["atomos"] > 5 and r["t_min"]]
    if not pts:
        return None
    xs = [math.log(n) for n, _ in pts]; ys = [math.log(t) for _, t in pts]
    if max(xs) - min(xs) > 0.2:
        mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
        a = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)
        a = min(max(a, 2.5), 4.5)
    else:
        a = 3.0
    k = sorted(t / n ** a for n, t in pts)[len(pts) // 2]
    modelo = lambda r: k * r["atomos"] ** a * (1.4 if r["mult"] == 2 else 1.0)
    med = lambda v: sorted(v)[len(v) // 2]
    feitas = [r for r in q if r["pronto"] and r["atomos"] > 5 and r["t_min"]]
    # correção empírica por multiplicidade: razão observado/modelo nas espécies prontas da maior molécula
    corr = {}
    for mult in (1, 2):
        g = [r for r in feitas if r["mult"] == mult]
        if g:
            nmax = max(r["atomos"] for r in g)
            corr[mult] = med([r["t_min"] / modelo(r) for r in g if r["atomos"] >= nmax - 3])
    out = {}
    for r in q:
        if r["pronto"]:
            continue
        mesmas = [x["t_min"] for x in feitas if x["mol"] == r["mol"] and x["mult"] == r["mult"]]
        out[r["id"]] = med(mesmas) if mesmas else modelo(r) * corr.get(r["mult"], 1.0)
    return out


def relogio(q, est):
    """Âncoras do relógio de término (epoch ms). Com 1 worker em série, a espécie atual começou quando o último
    JSON do DFT foi gravado; o término = esse instante + soma das estimativas pendentes."""
    if not est:
        return None
    jsons = [f for f in (ESP / "dft").glob("*.json") if f.stem != "H_atomo"]
    if not jsons:
        return None
    ultimo = max(f.stat().st_mtime for f in jsons)
    pend = [r for r in q if not r["pronto"]]
    atual = pend[0]
    return {
        "inicio_atual_ms": int(ultimo * 1000),
        "fim_ms": int((ultimo + sum(est.values()) * 60) * 1000),
        "atual": atual["id"], "atual_nome": f'{atual["molecula"]} · {atual["especie"]}'
                                            + (f' ({atual["sitio"]})' if atual["sitio"] != "—" else ""),
        "atual_est_min": round(est[atual["id"]], 1),
        "pendentes": len(pend),
    }


M3 = AQUI.parent / "tddft_uv"
M3_EXTRAS = {"epicatequina_oquinona": "o-Quinona da epicatequina", "benzofenona_3": "Benzofenona-3",
             "oxido_linalol": "cis-Óxido de linalol", "desidrodicatequina_a": "Desidrodicatequina A"}


def visual(todas, df):
    """Estruturas 2D, geometrias 3D, dados do M3 (orbitais, espectros) e manifesto de volumes."""
    import visual_assets as va
    bde = {m: {} for m in MOLS}
    if not df.empty:
        for r in df.itertuples():
            bde[r.molecula][r.sitio] = float(r.BDE)
    svg = {}
    for m in MOLS:
        try:
            svg[m] = va.svg_2d(m, bde.get(m, {}))
        except Exception as e:   # nunca derruba o painel por causa de uma figura
            svg[m] = f"<!-- erro: {e} -->"
    geo = {k: va.xyz_str(v["simbolos"], v["xyz"], k) for k, v in todas.items()}
    for mid in M3_EXTRAS:
        g = M3 / "geometrias" / f"{mid}.json"
        if g.exists():
            x = json.loads(g.read_text())
            geo[mid] = va.xyz_str(x["simbolos"], x["xyz"], mid)
    sitio_atomo = {m: {MOLS[m][1][int(mp)]: s["atomo"] for mp, s in todas[m]["sitios"].items()}
                   for m in MOLS if m in todas}
    # M3: orbitais de fronteira e espectros
    uv = []
    grade = list(range(200, 601, 2))
    try:
        import sys
        sys.path.insert(0, str(M3))
        from tddft_uv import espectro
        for f in sorted((M3 / "resultados").glob("*.json")):
            r = json.loads(f.read_text())
            y = espectro(r["e_ev"], r["f"], __import__("numpy").array(grade, float))
            uv.append({"id": r["id"], "nome": r["nome"], "papel": r["papel"], "homo": r.get("homo_ev"),
                       "lumo": r.get("lumo_ev"), "gap": r.get("gap_ev"),
                       "estados": [{"nm": round(l, 1), "f": round(fi, 4), "carater": c}
                                   for l, fi, c in zip(r["lambda_nm"], r["f"],
                                                       r.get("transicao_principal", [""] * len(r["f"])))],
                       "curva": [round(float(v), 4) for v in y]})
    except Exception as e:
        print("aviso M3:", e)
    try:
        vols = va.manifesto() if va.VOL.exists() else {}
    except Exception as e:
        print("aviso volumes:", e); vols = {}
    return {"svg2d": svg, "geo": geo, "sitio_atomo": sitio_atomo, "uv": uv, "uv_grade": grade,
            "vols": vols, "extras": M3_EXTRAS}


def _seg(s: str) -> float:
    """'1-02:03:04', '15:21:56', '303:22.40' (ps etime/time do macOS) → segundos."""
    d, _, s = s.strip().rpartition("-")
    p = [float(x) for x in s.split(":")]
    while len(p) < 3:
        p.insert(0, 0.0)
    return (int(d) if d else 0) * 86400 + p[0] * 3600 + p[1] * 60 + p[2]


def m3_status():
    """Fila do TD-DFT (M3) para o relógio do painel: prontas, molécula atual, tempo de CPU × relógio do processo
    e estimativa de término (t ∝ N³, calibrada nas moléculas prontas; None até a 1ª terminar)."""
    import subprocess
    import sys
    try:
        sys.path.insert(0, str(M3))
        import tddft_uv as tu
    except Exception as e:
        print("aviso M3 status:", e)
        return None
    itens = []
    for mid in tu.ORDEM:
        g, r = M3 / "geometrias" / f"{mid}.json", M3 / "resultados" / f"{mid}.json"
        d = json.loads(r.read_text()) if r.exists() else None
        simb = json.loads(g.read_text())["simbolos"] if g.exists() else None
        itens.append({"id": mid, "nome": tu.MOLS[mid][0], "atomos": len(simb) if simb else None,
                      # nº de funções de base 6-31+G(d) (esférica): 18 por átomo pesado (C, N, O) + 2 por H
                      "nao": sum(2 if s == "H" else 18 for s in simb) if simb else None,
                      "df": d.get("df", False) if d else None,
                      "pronto": d is not None, "t_min": round(d["t_total_s"] / 60, 1) if d else None,
                      "fim_s": r.stat().st_mtime if d else None})
    # processo Python do tddft_uv.py (o .pid guarda o caffeinate, que não acumula CPU)
    proc = None
    try:
        out = subprocess.run(["ps", "-Ao", "pid=,time=,etime=,command="], capture_output=True, text=True).stdout
        for ln in out.splitlines():
            pid, cpu, et, cmd = ln.split(None, 3)
            if "tddft_uv.py" in cmd and not cmd.startswith(("caffeinate", "nice", "zsh", "/bin/zsh")):
                proc = {"pid": int(pid), "cpu_s": _seg(cpu), "wall_s": _seg(et)}
                break
    except Exception as e:
        print("aviso M3 ps:", e)
    agora = dt.datetime.now().timestamp()
    pend = [x for x in itens if not x["pronto"]]
    atual = pend[0] if (pend and proc) else None
    fins = [x["fim_s"] for x in itens if x["fim_s"]]
    inicio = max(fins) if fins else (agora - proc["wall_s"] if proc else None)
    # custo ∝ nao³; calibra só no modo em uso (ajuste de densidade, se já houver molécula pronta com ele) e pelo
    # MENOR k observado: o tempo é de relógio, e o Mac dormindo ou um solver lento só aumentam k.
    feitos = [x for x in itens if x["pronto"] and x["nao"] and x["t_min"]]
    if any(x["df"] for x in feitos):
        feitos = [x for x in feitos if x["df"]]
    est = None
    if feitos and pend:
        k = min(x["t_min"] / x["nao"] ** 3 for x in feitos)
        est = {x["id"]: k * x["nao"] ** 3 for x in pend if x["nao"]}
    return {
        "total": len(itens), "prontas": len(itens) - len(pend), "itens": [{k: v for k, v in x.items() if k != "fim_s"}
                                                                          for x in itens],
        "rodando": proc is not None, "cpu_s": proc and round(proc["cpu_s"]), "wall_s": proc and round(proc["wall_s"]),
        "atual": atual and atual["id"], "atual_nome": atual and atual["nome"], "atual_atomos": atual and atual["atomos"],
        "atual_ordem": atual and itens.index(atual) + 1,
        "inicio_atual_ms": int(inicio * 1000) if (atual and inicio) else None,
        "atual_est_min": round(est[atual["id"]], 1) if (est and atual and atual["id"] in est) else None,
        "fim_ms": int((inicio + sum(est.values()) * 60) * 1000) if (est and inicio) else None,
    }


M7 = AQUI.parent / "estabilidade"


def estabilidade():
    """M7 (estabilidade das emulsões): partição/volatilidade (M7a) e compatibilidade com o emulsionante (M7c)."""
    import pandas as pd
    out = {}
    for k, f in [("part", "particao_volatilidade.csv"), ("val", "compat_validacao.csv"),
                 ("ost", "compat_ostwald.csv"), ("afin", "compat_afinidade.csv")]:
        try:
            out[k] = json.loads(pd.read_csv(M7 / f).to_json(orient="records", force_ascii=False))
        except Exception as e:
            print(f"aviso M7 ({f}):", e)
            out[k] = []
    return out


M56 = AQUI.parent / "seguranca_espaco"


def seguranca_espaco():
    """M5 (segurança tópica) + M6 (espaço químico): tabelas prontas de analise_m5_m6.py."""
    import pandas as pd
    out = {}
    for k, f in [("m5", "m5_seguranca.csv"), ("sim", "m6_similaridade.csv"), ("pca", "m6_pca.csv")]:
        try:
            df = pd.read_csv(M56 / f)
            out[k] = json.loads(df.drop(columns=[c for c in ("smiles",) if c in df]).to_json(
                orient="records", force_ascii=False))
        except Exception as e:
            print(f"aviso M5/M6 ({f}):", e)
            out[k] = []
    try:  # var, cargas, conj_ref, alertas, n_ref, ref_lista (gerados por analise_m5_m6.py)
        out.update(json.loads((M56 / "m56_extra.json").read_text(encoding="utf-8")))
    except Exception as e:
        print("aviso M5/M6 (m56_extra.json):", e)
    return out


M8 = AQUI.parent / "interacoes"


def interacoes():
    """M8 (cafeína–epicatequina): complexos únicos, energias, NCI. Os volumes NCI vão para painel_assets/vol/nci__<id>.js
    (carregados sob demanda pelo painel); o resto entra embutido. Retorna {"rodando": bool} enquanto não há resumo."""
    import subprocess
    rodando = False
    try:
        out = subprocess.run(["ps", "-Ao", "command="], capture_output=True, text=True).stdout
        rodando = any("interacoes.py" in ln and not ln.startswith(("caffeinate", "nice", "zsh", "/bin/zsh"))
                      for ln in out.splitlines())
    except Exception:
        pass
    n_part = len(list((M8 / "resultados" / "partidas").glob("*.json"))) if (M8 / "resultados" / "partidas").exists() else 0
    f = M8 / "resultados" / "resumo.json"
    if not f.exists():
        return {"rodando": rodando, "n_otimizadas": n_part}
    try:
        r = json.loads(f.read_text(encoding="utf-8"))
        import numpy as np
        from rdkit import Chem
        pt = Chem.GetPeriodicTable()
        simb = r["simbolos"]
        z = [pt.GetAtomicNumber(s) for s in simb]
        vol_dir = AQUI.parents[1] / "painel_assets" / "vol"
        vol_dir.mkdir(parents=True, exist_ok=True)
        cx = {c["id"]: c for c in r["complexos"]}
        nci = {}
        for cid, d in r["nci"].items():
            alvo = vol_dir / f"nci__{cid}.js"
            if not alvo.exists() or alvo.stat().st_mtime < f.stat().st_mtime:
                pac = d["vol"] | {"z": z, "xyz_bohr": (np.asarray(cx[cid]["xyz"]) / 0.52917721092).round(4).tolist()}
                alvo.write_text(f"(window.VOL=window.VOL||{{}})[{json.dumps('nci__' + cid)}]={json.dumps(pac)};")
            nci[cid] = {"scatter": d["scatter"], "integrais": d["integrais"]}
        mostrar = {c["id"] for c in r["complexos"][:12]} | set(nci)
        comp = []
        for c in r["complexos"]:
            item = {k: c.get(k) for k in ("id", "partida", "tipo", "dE", "dH", "dG", "TdS", "n_imag", "pi", "hb",
                                          "catecol_livre", "convergiu")}
            if c["id"] in mostrar:
                item["xyz"] = f"{len(simb)}\n{c['id']}\n" + "\n".join(
                    f"{s} {x:.4f} {y:.4f} {w:.4f}" for s, (x, y, w) in zip(simb, c["xyz"]))
            comp.append(item)
        fd = M8 / "resultados" / "refino_dft.json"
        dft_ref = json.loads(fd.read_text(encoding="utf-8")) if fd.exists() else None
        n_dft = len(list((M8 / "resultados" / "dft").glob("*.json"))) if (M8 / "resultados" / "dft").exists() else 0
        return {"dft": dft_ref, "n_dft": n_dft, "rodando": rodando, "n_otimizadas": n_part, "n_partidas": r["n_partidas"], "n_validos": r["n_validos"],
                "n_caf": r["n_caf"], "nivel": r["nivel"], "complexos": comp, "nci": nci, "boltz": r["boltzmann_dE"],
                "gerado": dt.datetime.fromtimestamp(f.stat().st_mtime).strftime("%d/%m/%Y %H:%M")}
    except Exception as e:
        print("aviso M8:", e)
        return {"rodando": rodando, "n_otimizadas": n_part, "erro": str(e)}


def main():
    todas, dfts = carregar()
    df, avisos = calcular(todas, dfts)
    q = fila(todas, dfts)
    n_esp = {m: sum(1 for r in q if r["mol"] == m) for m in MOLS}
    n_ok = {m: sum(1 for r in q if r["mol"] == m and r["pronto"]) for m in MOLS}
    mols = []
    for m in [x for x in ORDEM if x in MOLS]:
        nome, grupo, amostra, dpph = META[m]
        mols.append({"id": m, "nome": nome, "grupo": grupo, "amostra": amostra, "dpph": dpph,
                     "feitas": n_ok[m], "total": n_esp[m], "completa": n_ok[m] == n_esp[m]})
    sitios = [] if df.empty else json.loads(df.to_json(orient="records", force_ascii=False))
    est = estimativa(q)
    resto_h = None if est is None else sum(est.values()) / 60
    dados = {
        "gerado": dt.datetime.now().strftime("%d/%m/%Y %H:%M"),
        "pronto": sum(r["pronto"] for r in q), "total": len(q),
        "restante_h": None if resto_h is None else round(resto_h, 1),
        "relogio": relogio(q, est),
        "m3": m3_status(),
        "mols": mols, "sitios": sitios, "fila": q, "grupos": GRUPOS, "avisos": avisos,
        "const": {"H_proton": round(H_PROTON, 2), "H_eletron": round(H_ELETRON, 2)},
        "vis": visual(todas, df),
        "m7": estabilidade(),
        "m56": seguranca_espaco(),
        "m8": interacoes(),
        "metodo": __import__("metodologia").M,
    }
    html = TEMPLATE.replace("__3DMOL__", (AQUI / "vendor" / "3Dmol-min.js").read_text(encoding="utf-8")
                            .replace("</script", "<\\/script")).replace("__DADOS__", json.dumps(dados, ensure_ascii=False, allow_nan=False,
                                                      default=lambda o: None).replace("</", "<\\/"))
    SAIDA.write_text(html, encoding="utf-8")
    print(f"dashboard: {SAIDA}  ({dados['pronto']}/{dados['total']} espécies)")


TEMPLATE = (Path(__file__).resolve().parent / "template_dashboard.html").read_text(encoding="utf-8")


if __name__ == "__main__":
    main()
