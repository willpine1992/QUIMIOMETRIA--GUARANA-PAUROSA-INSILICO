# QUIMIOMETRIA · GUARANÁ & PAU-ROSA *IN SILICO*

**Por que só o extrato etanólico de guaraná neutraliza radicais livres?**

Estudo computacional (química quântica e quimioinformática) que procura, dentro das moléculas, a explicação para
os resultados experimentais de emulsões fitocosméticas feitas com extratos de sementes de guaraná
(*Paullinia cupana*) e óleo essencial de pau-rosa (*Aniba rosaeodora*) (FONSECA JÚNIOR, 2021).

## 🔗 Painel interativo

**https://willpine1992.github.io/QUIMIOMETRIA--GUARANA-PAUROSA-INSILICO/**

O painel conta a história em capítulos: cada um tem uma pergunta central, uma resposta curta gerada a partir dos
cálculos, a evidência (gráficos, tabelas, moléculas em 3D) e o "por que importa". Os botões **i**, **ii** e **iii**
de cada capítulo trazem a explicação, a metodologia detalhada para o artigo e os principais insights para a
discussão. A página tem uma tela de senha (a senha é fornecida pelo autor do estudo).

> A tela de senha só esconde a página; não criptografa os dados. O repositório é público.

## O enigma

| Material | Ensaio DPPH | Composição |
|---|---|---|
| Extrato **etanólico** de guaraná | ✓ ativo | fenólicos e flavonoides (taninos condensados) |
| Extrato **diclorometânico** de guaraná | ✗ inativo | rico em cafeína |
| **Óleo essencial** de pau-rosa | ✗ inativo | até 93,6 % de linalol |

**Hipótese:** as flavan-3-óis do extrato etanólico ((−)-epicatequina, unidade das procianidinas) têm um anel
catecol que cede o átomo de H com pouca energia; cafeína e linalol não têm esse grupo.

## Capítulos e módulos

| Cap. | Pergunta | Módulo | Método |
|---|---|---|---|
| 0 | O método reproduz um valor experimental conhecido? | M2 | calibração com o fenol |
| 1–4 | Quem cede H•, por qual mecanismo, em qual hidroxila, e isso bate com o DPPH? | M2 | GFN2-xTB/ALPB + M06-2X/6-311+G(d,p)/IEF-PCM (metanol): BDE, IP, PDE, PA, ETE |
| 5 | Onde estão os elétrons mais disponíveis, e que luz cada molécula absorve? | M3 | TD-CAM-B3LYP/6-31+G(d), PCM etanol; HOMO, LUMO e MEP em 3D |
| 6 | Depois de doar o H, o radical espalha o elétron? | M2b | densidade de spin (UKS M06-2X) |
| 7 | As moléculas explicam a instabilidade das emulsões? | M7 | partição e volatilidade (xTB/ALPB), UNIFAC, Ostwald |
| 8 | O que atravessa a pele e o que pede cuidado? | M5 | Potts–Guy, alertas estruturais, alérgenos (Reg. CE 1223/2009) |
| 9 | Esses constituintes se parecem com os ativos que a indústria usa? | M6 | ECFP4/Tanimoto, PCA, tamanho do cromóforo |
| 10 | A cafeína se prende à epicatequina e cobre o catecol? | M8 / M8b | busca de complexos GFN2-xTB/ALPB(água), qRRHO, índice NCI; refino ωB97X-V/def2-SVP + BSSE |

## Principais resultados até aqui

- **O catecol da epicatequina é o melhor doador de H• dos extratos:** BDE(3'-OH) = 81,7 kcal/mol, contra
  88,7 do fenol no mesmo nível de teoria. Cafeína e teobromina ficam acima de 109 kcal/mol.
- **Só os polifenóis completam as duas etapas do mecanismo SPLET** em metanol, o solvente do ensaio.
- **Os constituintes não se parecem com filtros UV** (maior Tanimoto 0,52; mediana 0,16); os flavan-3-óis são
  os mais próximos de antioxidantes cosméticos já usados.
- **Cafeína e epicatequina se associam por empilhamento π, mas de forma fraca:** anel A e anel B empatados
  (ΔE ≈ −13 kcal/mol com xTB), ΔG(1 M) de +1 a +2 kcal/mol; o catecol fica livre em ~35 % dos arranjos e a maior
  parte da epicatequina deve estar livre no extrato.

Os cálculos do TD-DFT (M3), do M2b e do refino M8b ainda estão em andamento; o painel mostra o estado de cada
capítulo ("respondida", "parcial" ou "aguardando").

## Estrutura

```
biblioteca/                  M1 · biblioteca curada de 76 constituintes (LOTUS + literatura)
modelagem/
  dft_antioxidantes/         M2/M2b · DFT dos antioxidantes; dashboard.py e template do painel; metodologia.py
  tddft_uv/                  M3 · TD-DFT, espectros UV, orbitais (.npz)
  estabilidade/              M7 · partição, volatilidade, compatibilidade com o emulsionante
  seguranca_espaco/          M5/M6 · segurança tópica e espaço químico
  interacoes/                M8/M8b · complexos cafeína–epicatequina, NCI, refino DFT
painel_assets/vol/           volumes 3D (orbitais, MEP, spin, NCI) carregados sob demanda pelo painel
dashboard_dft_antioxidantes.html   painel autocontido
rodar_*.sh                   filas de cálculo (agendar | status | parar)
proposta_quimioinformatica_edson.md   proposta e estado de cada módulo
```

## Reprodução

- Python 3.14, PySCF 2.14, RDKit 2026.03, tblite 0.8 (GFN2-xTB), geomeTRIC 1.1, NumPy, pandas, 3Dmol.js 2.4.
- Os scripts importam `quimiolib`, uma biblioteca de apoio compartilhada entre projetos do autor
  (`_compartilhado/quimiolib`, fora deste repositório): `qm.py` concentra as rotinas xTB, DFT e TD-DFT.
- Tudo roda em CPU (macOS). Cada módulo é retomável: cada cálculo pronto vira um `.json` e não é refeito.
- O painel é gerado por `modelagem/dft_antioxidantes/dashboard.py` a partir dos resultados.

## Referência experimental

FONSECA JÚNIOR, E. Q. **Desenvolvimento de emulsões com atividades antimicrobiana, antioxidante e fotoprotetora,
utilizando extratos de guaraná (*Paullinia cupana*) e óleo essencial de pau-rosa (*Aniba rosaeodora*)**. 2021.
Dissertação (Mestrado em Biotecnologia e Recursos Naturais da Amazônia) – Universidade do Estado do Amazonas,
Manaus, 2021.

Nenhum dado numérico da dissertação é usado como entrada dos cálculos; apenas as tendências qualitativas
(ativo/inativo, FPS, estabilidade) entram na comparação. O texto da dissertação não faz parte deste repositório.
