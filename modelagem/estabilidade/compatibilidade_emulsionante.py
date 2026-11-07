"""M7c — Compatibilidade dos ativos com o sistema emulsionante (estudo Edson).

Pergunta: a separação de fase após o estresse térmico e o aumento de viscosidade com os ativos (FONSECA JÚNIOR, 2021)
são compatíveis com (i) a solubilidade em água dos constituintes do OE, que alimenta o amadurecimento de Ostwald, e
(ii) a afinidade dos polifenóis pelas cadeias de PEG do tensoativo?

Protocolo (contribuição de grupos; qualitativo)
  1. EHL de Griffin da Ceteareth-20 a partir da estrutura (C16/C18 + 20 EO): EHL = 20 · M(PEG) / M(total).
  2. UNIFAC original (VLE; grupos por ugropy, modelo thermo) a 298,15 K, soluto a diluição infinita.
     Validação: log Kow = log[(γ∞_água / γ∞_octanol úmido) · (V_água / V_octanol úmido)], octanol com x_água = 0,275,
     contra log Kow experimental (PubChem). O Dortmund foi testado e descartado (erro maior nos terpenoides).
  3. Amadurecimento de Ostwald (teoria LSW): taxa ∝ D · c∞ · Vm. Com D semelhante entre moléculas pequenas,
     a taxa relativa ∝ solubilidade em fração de volume φ = x_sol · Vm / V_água. x_sol experimental (PubChem) quando
     existe; senão 1/γ∞_água (UNIFAC). Vm do modelo de Abdulelah-Gani (ugropy).
  4. Afinidade: ln γ∞ do soluto em cada componente da formulação (quanto menor, mais favorável):
     água, óleo mineral (n-hexadecano), miristato de isopropila, dimeticona (trissiloxano), álcool cetílico,
     segmento de PEG da Ceteareth-20 (tetraglima), propilenoglicol e glicerina.
  Limitação: o UNIFAC é pouco confiável para polifenóis polifuncionais (efeitos de proximidade entre OH); xantinas
  e quercetina não têm grupos UNIFAC e ficam de fora (ver M7a).

Uso
  "<venv>/bin/python" compatibilidade_emulsionante.py   → compat_validacao.csv, compat_ostwald.csv, compat_afinidade.csv
"""
from __future__ import annotations

import math
from pathlib import Path

import pandas as pd
import ugropy as ug
from thermo.unifac import UFIP, UFSG, UNIFAC

AQUI = Path(__file__).resolve().parent
T = 298.15
V_AGUA = 0.018069          # L/mol
V_OCT_UMIDO = 0.725 * 0.1576 + 0.275 * V_AGUA

SOLV = {  # componentes da formulação (Tabela 10 de FONSECA JÚNIOR, 2021) e o modelo molecular usado
    "água": "O",
    "óleo mineral (hexadecano)": "CCCCCCCCCCCCCCCC",
    "miristato de isopropila": "CCCCCCCCCCCCCC(=O)OC(C)C",
    "dimeticona (trissiloxano)": "C[Si](C)(C)O[Si](C)(C)O[Si](C)(C)C",
    "álcool cetílico": "CCCCCCCCCCCCCCCCO",
    "PEG da Ceteareth-20 (tetraglima)": "COCCOCCOCCOCCOC",
    "propilenoglicol": "CC(O)CO",
    "glicerina": "OCC(O)CO",
}
SOLUTOS = {
    "Linalol": "CC(=CCCC(C)(C=C)O)C",
    "cis-Óxido de linalol (furanoide)": "C[C@]1(CC[C@@H](O1)C(C)(C)O)C=C",
    "α-Terpineol": "CC1=CCC(CC1)C(C)(C)O",
    "Geraniol": "CC(=CCC/C(=C/CO)/C)C",
    "Limoneno": "CC1=CCC(CC1)C(=C)C",
    "Benzoato de benzila": "C1=CC=C(C=C1)COC(=O)C2=CC=CC=C2",
    "Hidroperóxido de linalol (7-OOH)": "CC(C)(/C=C/CC(C)(C=C)O)OO",
    "(−)-Epicatequina": "C1[C@H]([C@H](OC2=CC(=CC(=C21)O)O)C3=CC(=C(C=C3)O)O)O",
    "(+)-Catequina": "C1[C@@H]([C@H](OC2=CC(=CC(=C21)O)O)C3=CC(=C(C=C3)O)O)O",
    "o-Quinona da epicatequina": "C1[C@H]([C@H](OC2=CC(=CC(=C21)O)O)C3=CC(=O)C(=O)C=C3)O",
    "Cafeína": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Quercetina": "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O",
}
OLEOS = {"Óleo mineral (hexadecano)": "CCCCCCCCCCCCCCCC", "Miristato de isopropila": "CCCCCCCCCCCCCC(=O)OC(C)C"}

# PubChem (seções LogP e Solubility), consultado em 04/11/2026
EXP_LOGP = {"Linalol": 2.97, "Geraniol": 3.56, "α-Terpineol": 2.98, "Benzoato de benzila": 3.97, "Limoneno": 4.57,
            "Fenol": 1.46, "Catecol": 0.88, "1-Hexanol": 2.03}
VALID_EXTRA = {"Fenol": "Oc1ccccc1", "Catecol": "Oc1ccccc1O", "1-Hexanol": "CCCCCCO"}
EXP_SOL_MG_L = {"Linalol": (1590, 154.25), "Geraniol": (100, 154.25), "Limoneno": (7.57, 136.24),
                "Óleo mineral (hexadecano)": (2.1e-5, 226.45)}


def grupos(smi):
    r = ug.unifac.get_groups(smi, "smiles")
    return r.subgroups_num if r.subgroups else None


def vm_l(smi):
    return ug.abdulelah_gani.get_groups(smi, "smiles").liquid_molar_volume.to("L/mol").magnitude


def gama_inf(g_sol, comps, xs):
    m = UNIFAC.from_subgroups(T=T, xs=[1e-9] + [x * (1 - 1e-9) for x in xs], chemgroups=[g_sol] + comps,
                              version=0, interaction_data=UFIP, subgroups=UFSG)
    return m.gammas()[0]


def main():
    G = {k: grupos(s) for k, s in SOLV.items()}
    sem = [k for k, g in G.items() if not g]
    assert not sem, f"solventes sem grupos UNIFAC: {sem}"
    g_w, g_oct = G["água"], grupos("CCCCCCCCO")

    # 1. EHL de Griffin
    C, H, O = 12.011, 1.008, 15.999
    peg = 20 * (2 * C + 4 * H + O)
    ehl = {n: 20 * peg / (n * C + (2 * n + 2) * H + O + peg) for n in (16, 18)}
    print(f"EHL Griffin Ceteareth-20: C16 {ehl[16]:.2f} · C18 {ehl[18]:.2f}")

    # 2. validação log Kow
    val = []
    for nome, smi in {**SOLUTOS, **VALID_EXTRA}.items():
        if nome not in EXP_LOGP or not (g := grupos(smi)):
            continue
        kow = math.log10(gama_inf(g, [g_w], [1]) / gama_inf(g, [g_oct, g_w], [0.725, 0.275]) * V_AGUA / V_OCT_UMIDO)
        val.append({"molecula": nome, "logKow_unifac": round(kow, 2), "logKow_exp": EXP_LOGP[nome]})
    val = pd.DataFrame(val)
    val["erro"] = (val.logKow_unifac - val.logKow_exp).round(2)
    val.to_csv(AQUI / "compat_validacao.csv", index=False)
    print(val.to_string(index=False))
    print(f"MAE = {val.erro.abs().mean():.2f} · viés = {val.erro.mean():+.2f}\n")

    # 3. Ostwald (taxa relativa ∝ solubilidade em fração de volume)
    ost = []
    for nome, smi in {**{k: v for k, v in SOLUTOS.items() if k in ("Linalol", "cis-Óxido de linalol (furanoide)",
                                                                   "α-Terpineol", "Geraniol", "Limoneno",
                                                                   "Benzoato de benzila")}, **OLEOS}.items():
        g = grupos(smi)
        x_unifac = 1 / gama_inf(g, [g_w], [1]) if g else None
        x_exp = None
        if nome in EXP_SOL_MG_L:
            mg, mm = EXP_SOL_MG_L[nome]
            x_exp = (mg / 1000 / mm) / (1 / V_AGUA)
        x = x_exp if x_exp is not None else x_unifac
        vm = vm_l(smi)
        ost.append({"molecula": nome, "x_sol_unifac": x_unifac, "x_sol_exp": x_exp, "Vm_L_mol": round(vm, 4),
                    "phi_sol": x * vm / V_AGUA, "fonte_x": "PubChem" if x_exp is not None else "UNIFAC"})
    ost = pd.DataFrame(ost)
    ref = ost.loc[ost.molecula == "Óleo mineral (hexadecano)", "phi_sol"].item()
    ost["log_taxa_rel_oleo_mineral"] = [round(math.log10(p / ref), 1) for p in ost.phi_sol]
    ost.to_csv(AQUI / "compat_ostwald.csv", index=False)
    print(ost.to_string(index=False), "\n")

    # 4. afinidade (ln γ∞ em cada componente)
    af = []
    for nome, smi in SOLUTOS.items():
        g = grupos(smi)
        linha = {"molecula": nome, "tem_grupos_unifac": bool(g)}
        if g:
            for k in SOLV:
                linha[k] = round(math.log(gama_inf(g, [G[k]], [1])), 2)
        af.append(linha)
    af = pd.DataFrame(af)
    af.to_csv(AQUI / "compat_afinidade.csv", index=False)
    pd.set_option("display.width", 250)
    print(af.to_string(index=False))


if __name__ == "__main__":
    main()
