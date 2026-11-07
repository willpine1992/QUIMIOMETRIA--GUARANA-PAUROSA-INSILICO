"""Biblioteca do projeto Edson (guaraná + pau-rosa em emulsões) a partir das bibliotecas compartilhadas.

As bibliotecas ficam em QUIMIOINFORMÁTICA/_compartilhado/bibliotecas/ (pacote `quimiolib`).
Aqui só se acrescenta o que é específico desta dissertação:
  - fracao_tese: em qual amostra da tese (Hex → DCM → EtOH, OE, formulação) a molécula provavelmente está.
    INFERÊNCIA por classe/polaridade + CCD da tese (Dragendorff: cafeína no DCM e EtOH;
    FeCl3/AlCl3: fenólicos e flavonoides no EtOH). Não é dado medido.
Para atualizar as bibliotecas, rode o build.py de cada uma; depois rode este script.

Rodar:  "QUIMIOINFORMÁTICA/.venv/bin/python" montar_biblioteca_edson.py
"""
from pathlib import Path

import pandas as pd

from quimiolib import biblioteca, descritores

OUT = Path(__file__).resolve().parent / "saida"
OUT.mkdir(exist_ok=True)

FRACAO = {
    "proantocianidina (tanino condensado)": "EtOH",
    "flavan-3-ol": "EtOH",
    "metilxantina (alcaloide)": "DCM/EtOH",
    "ácido graxo (óleo da semente)": "Hex",
    "volátil da semente (traço)": "Hex",
    "hidroperóxido (alérgeno de contato)": "formulação (estresse térmico/luz)",
}

lib = biblioteca.carregar("paullinia_cupana_semente", "aniba_rosaeodora_OE", "padroes_referencia")
lib["especie"] = lib["biblioteca"].map({"paullinia_cupana_semente": "Paullinia cupana",
                                        "aniba_rosaeodora_OE": "Aniba rosaeodora",
                                        "padroes_referencia": "padrão"})
lib["fracao_tese"] = [
    "padrão" if b == "padroes_referencia" else FRACAO.get(c, "OE" if b == "aniba_rosaeodora_OE" else "?")
    for b, c in zip(lib["biblioteca"], lib["classe"])]

# mesma molécula em mais de uma espécie (ex.: cariofileno) → uma linha, espécies/frações unidas
juntar = lambda s: "; ".join(sorted({str(x) for x in s.dropna()}))
extra = lib.groupby("inchikey").agg(especies=("especie", juntar), fracoes=("fracao_tese", juntar),
                                    ids_compartilhados=("id", juntar))
edson = lib.drop_duplicates("inchikey").merge(extra, on="inchikey")
edson = edson.drop(columns=["especie", "fracao_tese"]).rename(
    columns={"especies": "especie", "fracoes": "fracao_tese"})
front = ["id", "ids_compartilhados", "biblioteca", "nome", "especie", "classe", "fracao_tese", "teor_max_pct"]
edson = edson[front + [c for c in edson.columns if c not in front]]

edson.to_csv(OUT / "biblioteca_edson.csv", index=False)
descritores.sdf_3d(edson, OUT / "biblioteca_edson_3d.sdf")
descritores.grade(edson[edson["biblioteca"] != "padroes_referencia"], OUT / "grade_estruturas.svg")
pd.concat([pd.read_csv(biblioteca.BIBLIOTECAS / b / "excluidos.csv").assign(biblioteca=b)
           for b in ("paullinia_cupana_semente", "aniba_rosaeodora_OE")]).to_csv(OUT / "excluidos.csv", index=False)

print(f"{len(edson)} estruturas únicas → {OUT}")
print(edson.groupby(["especie", "fracao_tese"]).size().to_string())
