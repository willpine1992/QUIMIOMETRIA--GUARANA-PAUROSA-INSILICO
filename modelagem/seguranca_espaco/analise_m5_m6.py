"""M5 (segurança tópica in silico) + M6 (espaço químico e similaridade) — estudo Edson.

Entradas
  biblioteca/saida/biblioteca_edson.csv            constituintes de guaraná (PC) e pau-rosa (AR)
  _compartilhado/bibliotecas/cosmeticos_referencia filtros UV orgânicos + antioxidantes cosméticos (COS)
Saídas (nesta pasta)
  m5_seguranca.csv          1 linha por constituinte: permeação, alertas, alérgeno UE
  m6_similaridade.csv       vizinho mais próximo (Tanimoto ECFP4) entre filtros e antioxidantes + cromóforo
  m6_pca.csv / m6_pca_cargas.csv
  para_modelos_externos.smi entrada pronta para VEGA / Pred-Skin / ProTox-3 (rodar manualmente)
  ../../dashboard_seguranca_espaco.html (via template_m5_m6.html)
Uso:  "<venv>/bin/python" analise_m5_m6.py
"""
from __future__ import annotations

import datetime as dt
import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from quimiolib import biblioteca, seguranca

AQUI = Path(__file__).resolve().parent
EDSON = AQUI.parents[1]
SAIDA_HTML = EDSON / "dashboard_seguranca_espaco.html"

DESC_PCA = [("massa_molar", "Massa molar"), ("logp_crippen", "logP"), ("tpsa", "TPSA"), ("hbd", "Doadores de LH"),
            ("hba", "Aceptores de LH"), ("ligacoes_rotacionaveis", "Lig. rotacionáveis"),
            ("aneis_aromaticos", "Anéis aromáticos"), ("fsp3", "Fsp3")]


def grupo(r) -> str:
    if str(r["classe"]).startswith("hidroperóxido"):
        return "oxid"
    return "guarana" if str(r["especie"]).startswith("Paullinia") else "paurosa"


def maior_conjugado(smiles: str) -> int:
    """Nº de átomos do maior sistema π conjugado (proxy de cromóforo; o M3 quantifica a absorção)."""
    m = Chem.MolFromSmiles(smiles)
    adj: dict[int, set[int]] = {}
    for b in m.GetBonds():
        if b.GetIsConjugated() or b.GetIsAromatic():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            adj.setdefault(i, set()).add(j); adj.setdefault(j, set()).add(i)
    vistos, melhor = set(), 0
    for a in adj:
        if a in vistos:
            continue
        pilha, comp = [a], 0
        vistos.add(a)
        while pilha:
            x = pilha.pop(); comp += 1
            for y in adj[x] - vistos:
                vistos.add(y); pilha.append(y)
        melhor = max(melhor, comp)
    if melhor == 0 and any(b.GetBondType() == Chem.BondType.DOUBLE for b in m.GetBonds()):
        melhor = 2  # ligação dupla isolada
    return melhor


def main():
    ed = pd.read_csv(EDSON / "biblioteca/saida/biblioteca_edson.csv")
    ed = ed[ed["especie"] != "padrão"].copy()
    ed["grupo"] = ed.apply(grupo, axis=1)
    cos = biblioteca.carregar("cosmeticos_referencia")
    cos["grupo"] = np.where(cos["classe"] == "filtro UV", "filtro", "antiox")

    # ---------------- M5 ----------------
    tab = seguranca.tabela_alertas().set_index("codigo")
    ed["alertas"] = ed["smiles"].map(seguranca.alertas)
    ed = seguranca.marcar_alergenos(ed)
    m5 = ed[["id", "nome", "grupo", "especie", "classe", "fracao_tese", "teor_max_pct", "massa_molar", "logp_crippen",
             "tpsa", "logKp_PottsGuy_cm_h", "alergeno_ue", "alergeno_tipo", "smiles"]].copy()
    m5["abaixo_500Da"] = m5["massa_molar"] < 500
    m5["alertas"] = ed["alertas"].map(lambda xs: "; ".join(xs))
    m5["sens_dominios"] = ed["alertas"].map(
        lambda xs: "; ".join(sorted({tab.loc[c, "dominio"] for c in xs if tab.loc[c, "desfecho"] == "sensibilização"})))
    m5["genotox"] = ed["alertas"].map(
        lambda xs: "; ".join(tab.loc[c, "rotulo"] for c in xs if tab.loc[c, "desfecho"] == "genotoxicidade"))
    m5["atencao"] = ed["alertas"].map(lambda xs: "; ".join(tab.loc[c, "rotulo"] for c in xs if tab.loc[c, "desfecho"] == "atenção"))
    m5.to_csv(AQUI / "m5_seguranca.csv", index=False)
    with open(AQUI / "para_modelos_externos.smi", "w") as f:
        for r in m5.itertuples():
            f.write(f"{r.smiles}\t{r.id}\n")

    # ---------------- M6 ----------------
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp = lambda s: gen.GetFingerprint(Chem.MolFromSmiles(s))
    fe, fc = [fp(s) for s in ed["smiles"]], [fp(s) for s in cos["smiles"]]
    sim = np.array([DataStructs.BulkTanimotoSimilarity(a, fc) for a in fe])
    linhas = []
    for i, r in enumerate(ed.itertuples()):
        d = {"id": r.id, "nome": r.nome, "grupo": r.grupo}
        for g, rot in (("filtro", "filtro"), ("antiox", "antiox")):
            idx = np.where(cos["grupo"].values == g)[0]
            j = idx[np.argmax(sim[i, idx])]
            d[f"sim_{rot}"] = round(float(sim[i, j]), 3)
            d[f"viz_{rot}"] = cos["nome"].iloc[j]
        d["conjugado"] = maior_conjugado(r.smiles)
        linhas.append(d)
    m6 = pd.DataFrame(linhas)
    m6.to_csv(AQUI / "m6_similaridade.csv", index=False)

    todos = pd.concat([ed.assign(origem="constituinte"), cos.assign(origem="referência")], ignore_index=True)
    X = StandardScaler().fit_transform(todos[[c for c, _ in DESC_PCA]].astype(float))
    pca = PCA(n_components=2, random_state=0).fit(X)
    Z = pca.transform(X)
    pts = pd.DataFrame({"id": todos["id"], "nome": todos["nome"], "grupo": todos["grupo"],
                        "pc1": Z[:, 0].round(3), "pc2": Z[:, 1].round(3)})
    pts.to_csv(AQUI / "m6_pca.csv", index=False)
    cargas = pd.DataFrame(pca.components_.T, columns=["pc1", "pc2"], index=[n for _, n in DESC_PCA]).round(3)
    cargas.to_csv(AQUI / "m6_pca_cargas.csv")
    conj_ref = {g: sorted(int(maior_conjugado(s)) for s in cos.loc[cos["grupo"] == g, "smiles"]) for g in ("filtro", "antiox")}

    # ---------------- dashboard ----------------
    def rec(df):
        return json.loads(df.to_json(orient="records", force_ascii=False))
    dados = {
        "gerado": dt.datetime.now().strftime("%d/%m/%Y %H:%M"),
        "m5": rec(m5.drop(columns=["smiles"])),
        "m6": rec(m6), "pca": rec(pts),
        "var": [round(float(v) * 100, 1) for v in pca.explained_variance_ratio_],
        "cargas": rec(cargas.reset_index().rename(columns={"index": "desc"})),
        "conj_ref": conj_ref,
        "alertas": rec(seguranca.tabela_alertas().drop(columns=["smarts"])),
        "n_ref": {"filtro": int((cos["grupo"] == "filtro").sum()), "antiox": int((cos["grupo"] == "antiox").sum())},
        "ref_lista": rec(cos[["id", "nome", "grupo", "papel"]]),
    }
    tpl = (AQUI / "template_m5_m6.html").read_text(encoding="utf-8")
    SAIDA_HTML.write_text(tpl.replace("__DADOS__", json.dumps(dados, ensure_ascii=False).replace("</", "<\\/")),
                          encoding="utf-8")
    print(f"M5: {len(m5)} constituintes · M6: {len(cos)} referências · PCA {dados['var']} % · → {SAIDA_HTML.name}")


if __name__ == "__main__":
    main()
