# Proposta de estudo quimioinformático complementar — guaraná e pau-rosa em fitocosméticos

**Estudo de referência:** FONSECA JÚNIOR, E. Q. *Desenvolvimento de emulsões com atividades antimicrobiana, antioxidante e fotoprotetora, utilizando extratos de guaraná (Paullinia cupana) e óleo essencial de pau-rosa (Aniba rosaeodora)*. 2021. Dissertação (Mestrado em Biotecnologia e Recursos Naturais da Amazônia) – Universidade do Estado do Amazonas, Manaus, 2021. Orientação: Dra. Patrícia M. Albuquerque; coorientação: Dr. Geverson F. da Silva.

**Versão:** 2 (03/11/2026). Substitui a versão 1, que propunha correções e reanálise dos dados.

---

## 1. Escopo

O trabalho de referência é de terceiros. Por isso, este estudo:

- **não corrige, não refaz e não reanalisa** nenhum dado de Fonseca Júnior (2021);
- é **100 % in silico**: modelagem molecular, química quântica e quimioinformática;
- usa os resultados de Fonseca Júnior (2021) apenas como **referência experimental qualitativa**, citados conforme a ABNT.

### 1.1 Resultados experimentais prévios usados como referência (só tendências)

Em estudos prévios, Fonseca Júnior (2021) observou as tendências abaixo, que servem de referência qualitativa.

| Atividade | Observação qualitativa usada |
|---|---|
| Antioxidante (DPPH, metanol) | Só o extrato **etanólico** de guaraná é ativo. Os extratos DCM e hexânico e o óleo de pau-rosa são praticamente inativos. A quercetina é o padrão mais ativo. |
| Fotoproteção (UV, Mansur) | Os extratos **etanólico e DCM** absorvem fortemente no UV; o óleo essencial absorve pouco. |
| Antimicrobiana (CIM) | O óleo de pau-rosa inibe gram-positivas e *C. tropicalis*, mas não *P. aeruginosa* nem *C. parapsilosis*. Os extratos DCM e EtOH inibem *Candida* e *E. faecalis*. O extrato hexânico é inativo. |
| Citotoxicidade | Nenhum extrato reduz a viabilidade de HCT116. Os extratos EtOH e Hex são tóxicos para *A. salina*. |
| Fitoquímica (CCD) | O Dragendorff indica cafeína nos extratos DCM e EtOH. FeCl₃ e AlCl₃ indicam fenólicos e flavonoides, sobretudo no EtOH. |

**Não serão usados como dado de entrada:** os valores absolutos de FPS acima de 200 a 0,1 %, os CE₅₀ extrapolados além da faixa testada e as CIMs para QSAR. Os motivos são o número pequeno de amostras e o fato de os ensaios serem feitos com misturas, não com moléculas isoladas. Esses números não são reinterpretados nem corrigidos; o estudo simplesmente não depende deles.

### 1.2 Perguntas que o estudo responde

1. **Antioxidante:** quais moléculas e qual mecanismo (HAT, SET-PT ou SPLET) explicam por que só o extrato etanólico sequestra radicais?
2. **Fotoproteção:** quais cromóforos dos extratos absorvem em UVB e UVA, e como se comparam à benzofenona-3?
3. **Antimicrobiana:** os constituintes majoritários têm afinidade plausível por alvos dos micro-organismos sensíveis, ou o efeito é mais compatível com ação sobre a membrana?
4. **Segurança tópica:** quais constituintes, incluindo os produtos de oxidação do linalol, merecem atenção quanto a permeação e sensibilização cutânea?
5. **Posicionamento:** onde os constituintes se situam no espaço químico em relação a filtros UV e antioxidantes cosméticos já usados?
6. **Estabilidade:** que mecanismos moleculares são compatíveis com as alterações vistas na estabilidade preliminar (perda de odor, mudança de cor, separação de fase)?

---

## 2. Módulos e estado

| # | Módulo | Pergunta | Estado |
|---|---|---|---|
| M1 | Biblioteca de constituintes | base de todas | ✅ **concluído** |
| M2 | DFT antioxidante (HAT/SET-PT/SPLET) | 1 | ⏳ **em execução** (geometrias prontas; DFT em fila) |
| M3 | TD-DFT, espectros UV | 2 | ⏳ **preparado e agendado** (inicia sozinho quando o M2 terminar; `rodar_tddft_uv.sh`) |
| M4 | Docking antimicrobiano | 3 | planejado |
| M5 | Segurança tópica in silico | 4 | ✅ **calculado** (`analise_m5_m6.py`; Cap. 8 do painel). Pendente: submeter `para_modelos_externos.smi` aos modelos públicos |
| M6 | Espaço químico e similaridade | 5 | ✅ **calculado** (biblioteca `cosmeticos_referencia` criada; PCA + Tanimoto; Cap. 9 do painel). UMAP não feito |
| M7 | Mecanismos de instabilidade das emulsões | 6 | ⏳ **em execução** (M7a ✅; M7c ✅; M7b na fila do M3, geometrias prontas) |

---

### M1 — Biblioteca de constituintes ✅

**Produto:** 76 estruturas únicas, com estrutura padronizada, proveniência e cerca de 30 descritores RDKit. Também têm OH fenólico, catecol, hidroperóxido e log *K*p de Potts–Guy. O arquivo 3D (SDF, ETKDGv3 + MMFF94) está pronto para modelagem.

| Biblioteca (compartilhada) | Estruturas | Fontes |
|---|---|---|
| `paullinia_cupana_semente` | 22 | LOTUS/Wikidata: metilxantinas, (±)-catequinas, procianidinas B1–B4, A2 e C1, lipídios e voláteis da semente |
| `aniba_rosaeodora_OE` | 48 | Xavier et al. 2021 (CG-EM de folhas e galhos), Teles et al. 2021 (OE), hidroperóxidos de linalol (Sköld et al. 2004) |
| `padroes_referencia` | 7 | Padrões usados por Fonseca Júnior (2021): quercetina, BP-3, levofloxacino, cetoconazol, ácido gálico, ácido cinâmico, doxorrubicina |

**Curadoria:**
- **Registros excluídos:** 12 registros do LOTUS atribuídos ao pau-rosa foram descartados. Oito eram de outras espécies de *Aniba* ou de cabreúva (erro de anotação no banco). Quatro eram constituintes da madeira, que não entram no OE. O khusinol foi excluído porque não tem estrutura no PubChem.
- **Isômeros corrigidos:** dois casos em que o número CAS levava ao isômero errado (*cis*-óxido de linalol e 14-hidroxi-9-*epi*-cariofileno).
- **Coluna `fracao_tese`:** indica em qual amostra de Fonseca Júnior (2021) a molécula provavelmente está (Hex, DCM, EtOH, OE ou formulação). É uma **inferência** por polaridade e pela CCD, não uma medida.
- **Estereoquímica no OE:** é a que o PubChem atribui, porque o CG-EM não determina configuração absoluta.

**Arquivos:**
- Bibliotecas: `QUIMIOINFORMÁTICA/_compartilhado/bibliotecas/`.
- Consolidação do projeto: `edson/biblioteca/montar_biblioteca_edson.py`, que gera `saida/biblioteca_edson.csv` e `biblioteca_edson_3d.sdf`.

---

### M2 — DFT antioxidante ⏳

**Hipótese:** o anel B catecol das flavan-3-óis (epicatequina e as unidades das procianidinas) tem BDE O–H baixa, comparável à da quercetina. Metilxantinas e linalol têm apenas C–H, N–H ou O–H alcoólico de energia alta. Isso explicaria, em nível molecular, a atividade exclusiva do extrato etanólico e a inatividade do DCM (rico em cafeína) e do OE (rico em linalol).

**Moléculas e sítios** (cada sítio é rotulado no SMILES):

| Molécula | Representa | Sítios |
|---|---|---|
| (−)-epicatequina | extrato EtOH; unidade das procianidinas; (+)-catequina é diastereoisômero | 3-OH, 5-OH, 7-OH, 3′-OH, 4′-OH |
| cafeína | extratos DCM/EtOH | C8–H |
| teobromina | extratos DCM/EtOH | N1–H, C8–H |
| teofilina | extratos DCM/EtOH | N7–H, C8–H |
| linalol | OE de pau-rosa (até 93,6 %) | O–H, C4–H alílico |
| quercetina | padrão antioxidante (FONSECA JÚNIOR, 2021) | 3-OH, 5-OH, 7-OH, 3′-OH, 4′-OH |
| ácido gálico | padrão Folin | 3-OH (= 5-OH), 4-OH |
| fenol | calibração | O–H |

**Protocolo de triagem: M06-2X/6-311+G(d,p)//GFN2-xTB, em metanol (solvente do ensaio DPPH).**

1. **Conformero:** ETKDGv3 (20 conformeros) + MMFF94, escolhendo o de menor energia GFN2-xTB/ALPB(metanol).
2. **Geometria e frequências:** GFN2-xTB/ALPB(metanol) com o otimizador geomeTRIC. A hessiana sai por diferenças finitas, e a correção térmica é RRHO a 298,15 K. Isso é feito para o neutro, cada radical, cada ânion (só nos sítios O–H e N–H) e o cátion-radical.
3. **Energia eletrônica:** single point M06-2X/6-311+G(d,p) com IEF-PCM (ε = 32,613), SCF direto, no PySCF.
4. **Descritores** (entalpias, kcal/mol):
   - HAT: **BDE** = H(R•) + H(H•) − H(RH)
   - SET-PT: **IP** = H(RH•⁺) + H(e⁻) − H(RH); **PDE** = H(R•) + H(H⁺) − H(RH•⁺)
   - SPLET: **PA** = H(R⁻) + H(H⁺) − H(RH); **ETE** = H(R•) + H(e⁻) − H(R⁻)
5. **Constantes:**
   - H(H•) é calculado no mesmo nível de teoria.
   - ΔsolvH(H⁺, MeOH) = −1038 kJ/mol (IEF-PCM; Rottmannová et al. 2013).
   - ΔsolvH(e⁻, MeOH) = −86 kJ/mol (Rimarčík et al. 2010; valor obtido por citação secundária).
   - IP, PDE, PA e ETE em valor absoluto dependem dessas constantes. Por isso também são reportados **Δ em relação ao fenol**, que não depende delas.
6. **Controle de qualidade:** número de frequências imaginárias por espécie, convergência do SCF e ⟨S²⟩ dos radicais (contaminação de spin).

**Calibração:** no teste inicial, a BDE do fenol deu **88,7 kcal/mol**; o valor experimental fica entre 87 e 88,7 kcal/mol. O teste usou *density fitting* e PCM com 302 pontos por átomo. A execução final usa SCF direto e 110 pontos por átomo, e recalcula o fenol nessas mesmas condições.

**Execução:**
- **Estado:** as 52 geometrias xTB estão prontas. As 53 espécies de DFT estão numa fila de 1 processo em prioridade baixa, molécula a molécula (das menores para as maiores), e levam de 12 a 15 h no total.
- **Comandos:** `edson/rodar_dft_antioxidantes.sh` com as opções `iniciar`, `status`, `parar` e `analise`. O cálculo retoma de onde parou se for interrompido.
- **Saídas:** `descritores_por_sitio.csv`, `resumo_por_molecula.csv` e `avisos.txt`, em `edson/modelagem/dft_antioxidantes/`.

#### Fila de execução do M2 (53 cálculos DFT)

A fila roda um cálculo por vez, em prioridade baixa (`nice`), para não atrapalhar o uso do computador. Cada molécula fecha completa antes da próxima começar, então os resultados vão saindo aos poucos. As maiores ficam no fim.

**Para que serve cada tipo de espécie:**

| Espécie | Carga / multiplicidade | Entra no cálculo de |
|---|---|---|
| neutro (RH) | 0 / singleto | todos os descritores |
| radical (R•), um por sítio | 0 / dubleto | BDE, PDE, ETE |
| ânion (R⁻), só sítios O–H e N–H | −1 / singleto | PA, ETE |
| cátion-radical (RH•⁺), um por molécula | +1 / dubleto | IP, PDE |
| átomo de H | 0 / dubleto | BDE (entalpia do H•) |

**Blocos da fila, com tempo estimado** (1 processo; o tempo da molécula pequena foi medido: fenol neutro, 13 átomos, 1,3 min):

| Ordem | Molécula | Espécies | Átomos | Tempo estimado | Resultado ao fechar o bloco |
|---|---|---|---|---|---|
| 1 | átomo de H | 1 | 1 | < 1 s | — |
| 2 | fenol (calibração) | 4 | 12–13 | ~6 min | confere o método (BDE exp. 87–88,7) |
| 3 | cafeína | 3 | 23–24 | ~20 min | C8–H: espera-se BDE alta → explica DCM inativo |
| 4 | teobromina | 5 | 20–21 | ~25 min | N1–H e C8–H |
| 5 | teofilina | 5 | 20–21 | ~25 min | N7–H e C8–H |
| 6 | ácido gálico | 5 | 17–18 | ~15 min | padrão Folin (referência fenólica) |
| 7 | linalol | 5 | 28–29 | ~1 h | O–H e C4–H alílico → explica OE inativo |
| 8 | epicatequina | 12 | 34–35 | ~5–7 h | catecol do anel B → explica EtOH ativo |
| 9 | quercetina | 12 | 31–32 | ~4–5 h | padrão positivo (F6) |
| | **total** | **52 + H** | | **~11–14 h** | |

**Fila detalhada** (estado em 03/11/2026, 21h40; atualizar com `./rodar_dft_antioxidantes.sh status`):

| # | Molécula | Espécie | Sítio | Carga | Multiplicidade | Átomos | Freq. imag. (xTB) | Estado |
|---|---|---|---|---|---|---|---|---|
| 1 | H | átomo H | — | 0 | 2 (dubleto) | 1 | — | ✅ |
| 2 | fenol | neutro | — | 0 | 1 (singleto) | 13 | 0 | ✅ 1,3 min |
| 3 | fenol | ânion | O–H | −1 | 1 (singleto) | 12 | 0 | pendente |
| 4 | fenol | cátion-radical | — | +1 | 2 (dubleto) | 13 | 0 | pendente |
| 5 | fenol | radical | O–H | 0 | 2 (dubleto) | 12 | 0 | pendente |
| 6 | cafeína | neutro | — | 0 | 1 (singleto) | 24 | 0 | pendente |
| 7 | cafeína | cátion-radical | — | +1 | 2 (dubleto) | 24 | **2** | pendente |
| 8 | cafeína | radical | C8–H | 0 | 2 (dubleto) | 23 | **1** | pendente |
| 9 | teobromina | neutro | — | 0 | 1 (singleto) | 21 | **1** | pendente |
| 10 | teobromina | ânion | N1–H | −1 | 1 (singleto) | 20 | **1** | pendente |
| 11 | teobromina | cátion-radical | — | +1 | 2 (dubleto) | 21 | **2** | pendente |
| 12 | teobromina | radical | N1–H | 0 | 2 (dubleto) | 20 | 0 | pendente |
| 13 | teobromina | radical | C8–H | 0 | 2 (dubleto) | 20 | **1** | pendente |
| 14 | teofilina | neutro | — | 0 | 1 (singleto) | 21 | 0 | pendente |
| 15 | teofilina | ânion | N7–H | −1 | 1 (singleto) | 20 | 0 | pendente |
| 16 | teofilina | cátion-radical | — | +1 | 2 (dubleto) | 21 | 0 | pendente |
| 17 | teofilina | radical | N7–H | 0 | 2 (dubleto) | 20 | 0 | pendente |
| 18 | teofilina | radical | C8–H | 0 | 2 (dubleto) | 20 | 0 | pendente |
| 19 | ácido gálico | neutro | — | 0 | 1 (singleto) | 18 | 0 | pendente |
| 20 | ácido gálico | ânion | 3-OH (= 5-OH) | −1 | 1 (singleto) | 17 | 0 | pendente |
| 21 | ácido gálico | ânion | 4-OH | −1 | 1 (singleto) | 17 | 0 | pendente |
| 22 | ácido gálico | cátion-radical | — | +1 | 2 (dubleto) | 18 | 0 | pendente |
| 23 | ácido gálico | radical | 3-OH (= 5-OH) | 0 | 2 (dubleto) | 17 | 0 | pendente |
| 24 | ácido gálico | radical | 4-OH | 0 | 2 (dubleto) | 17 | 0 | pendente |
| 25 | linalol | neutro | — | 0 | 1 (singleto) | 29 | 0 | pendente |
| 26 | linalol | ânion | O–H | −1 | 1 (singleto) | 28 | 0 | pendente |
| 27 | linalol | cátion-radical | — | +1 | 2 (dubleto) | 29 | 0 | pendente |
| 28 | linalol | radical | O–H | 0 | 2 (dubleto) | 28 | 0 | pendente |
| 29 | linalol | radical | C4–H (alílico) | 0 | 2 (dubleto) | 28 | 0 | pendente |
| 30 | epicatequina | neutro | — | 0 | 1 (singleto) | 35 | 0 | pendente |
| 31 | epicatequina | ânion | 3′-OH | −1 | 1 (singleto) | 34 | 0 | pendente |
| 32 | epicatequina | ânion | 4′-OH | −1 | 1 (singleto) | 34 | 0 | pendente |
| 33 | epicatequina | ânion | 3-OH | −1 | 1 (singleto) | 34 | 0 | pendente |
| 34 | epicatequina | ânion | 5-OH | −1 | 1 (singleto) | 34 | 0 | pendente |
| 35 | epicatequina | ânion | 7-OH | −1 | 1 (singleto) | 34 | 0 | pendente |
| 36 | epicatequina | cátion-radical | — | +1 | 2 (dubleto) | 35 | 0 | pendente |
| 37 | epicatequina | radical | 3′-OH | 0 | 2 (dubleto) | 34 | 0 | pendente |
| 38 | epicatequina | radical | 4′-OH | 0 | 2 (dubleto) | 34 | 0 | pendente |
| 39 | epicatequina | radical | 3-OH | 0 | 2 (dubleto) | 34 | 0 | pendente |
| 40 | epicatequina | radical | 5-OH | 0 | 2 (dubleto) | 34 | 0 | pendente |
| 41 | epicatequina | radical | 7-OH | 0 | 2 (dubleto) | 34 | 0 | pendente |
| 42 | quercetina | neutro | — | 0 | 1 (singleto) | 32 | 0 | pendente |
| 43 | quercetina | ânion | 3′-OH | −1 | 1 (singleto) | 31 | 0 | pendente |
| 44 | quercetina | ânion | 4′-OH | −1 | 1 (singleto) | 31 | 0 | pendente |
| 45 | quercetina | ânion | 3-OH | −1 | 1 (singleto) | 31 | 0 | pendente |
| 46 | quercetina | ânion | 5-OH | −1 | 1 (singleto) | 31 | 0 | pendente |
| 47 | quercetina | ânion | 7-OH | −1 | 1 (singleto) | 31 | 0 | pendente |
| 48 | quercetina | cátion-radical | — | +1 | 2 (dubleto) | 32 | 0 | pendente |
| 49 | quercetina | radical | 3′-OH | 0 | 2 (dubleto) | 31 | 0 | pendente |
| 50 | quercetina | radical | 4′-OH | 0 | 2 (dubleto) | 31 | 0 | pendente |
| 51 | quercetina | radical | 3-OH | 0 | 2 (dubleto) | 31 | 0 | pendente |
| 52 | quercetina | radical | 5-OH | 0 | 2 (dubleto) | 31 | 0 | pendente |
| 53 | quercetina | radical | 7-OH | 0 | 2 (dubleto) | 31 | 0 | pendente |

**Pendência de qualidade (geometrias xTB):** 6 espécies de cafeína e teobromina têm frequências imaginárias acima de 50 cm⁻¹ em módulo (em negrito na tabela). Isso indica que a geometria xTB não é um mínimo verdadeiro, provavelmente por causa da rotação de grupos metila. O efeito esperado é pequeno: só a correção térmica muda, em ≲ 0,5 kcal/mol, e a energia DFT não é afetada de modo relevante. Mesmo assim, antes da análise final essas 6 geometrias serão reotimizadas: deslocadas ao longo do modo imaginário, depois otimizadas com critério mais apertado e a correção térmica recalculada. O DFT delas será refeito se a geometria mudar. Por enquanto, os resultados dessas moléculas devem ser lidos como preliminares.

**Limitações:**
- É um protocolo de **triagem**. Para publicação, os sítios-chave (BDE mínima e PA mínima de cada molécula) devem ser refinados com otimização e frequências em DFT, por exemplo no ORCA.
- As procianidinas são representadas pelo monômero, porque um dímero ou trímero é inviável no hardware atual.
- O grupo COOH do ácido gálico não foi incluído no SPLET.
- O estudo é termodinâmico. A cinética (constantes de velocidade, protocolo QM-ORSA de Galano & Alvarez-Idaboy 2013) fica como extensão possível.

---

### M3 — TD-DFT: espectros UV (planejado)

**Objetivo:** identificar quais constituintes absorvem em UVB (290–320 nm) e UVA (320–400 nm). A comparação é qualitativa: os extratos EtOH e DCM absorvem e o OE quase não.

- **Moléculas:** epicatequina, cafeína, teobromina, teofilina, linalol, óxidos de linalol, quercetina e benzofenona-3. Opcionalmente, um dímero de procianidina se o custo permitir.
- **Método:** TD-DFT com CAM-B3LYP (ou ωB97X-D)/6-31+G(d) e IEF-PCM em etanol, solvente usado para os extratos no ensaio de Mansur. Serão calculados de 15 a 20 estados excitados sobre as geometrias xTB do M2. Antes, um teste com a BP-3, que tem λmax experimental conhecido.
- **Saídas:** espectros simulados com alargamento gaussiano, λmax, força do oscilador e área integrada em UVB e UVA. Também um λ crítico teórico, comparado ao da BP-3.
- **Limitação:** sem composição quantitativa dos extratos, não é possível reconstruir o espectro da mistura. A comparação é molécula a molécula.
- **Execução (preparada em 04/11/2026):**
  - Script `modelagem/tddft_uv/tddft_uv.py` e executor `rodar_tddft_uv.sh` (iniciar, agendar, status, parar, analise).
  - Nível: TD-CAM-B3LYP/6-31+G(d), IEF-PCM etanol (ε = 24,852), 15 singletos, PySCF.
  - Geometrias: as do M2 (GFN2-xTB/ALPB metanol) para linalol, xantinas, epicatequina e quercetina; novas otimizações (ALPB etanol) para benzofenona-3 e cis-óxido de linalol (furanoide).
  - Fila (1 por vez, nice 15): benzofenona-3 (calibração) → linalol → óxido de linalol → cafeína → teobromina → teofilina → epicatequina → quercetina.
  - Espectro: gaussianas em energia (FWHM 0,40 eV); saídas `resumo_uv.csv` (λmax, estado mais intenso, áreas UVC/UVB/UVA2/UVA1, λ crítico) e `espectros_uv.csv`.
  - A teobromina usa a geometria com 1 frequência imaginária; se a reotimização do M2 mudar a geometria, refazer o TD-DFT dela.

### M4 — Docking antimicrobiano (planejado)

**Objetivo:** avaliar a plausibilidade mecanística da atividade observada.

- **Alvos** (estruturas do PDB):
  - *S. aureus*: DNA girase B, FabI e sortase A.
  - *C. albicans*: CYP51 (lanosterol 14α-desmetilase).
  - *E. faecalis*: alvo escolhido pela disponibilidade de estrutura.
  - Controles positivos: levofloxacino na girase e cetoconazol na CYP51.
- **Método:**
  - AutoDock Vina sobre os SDF 3D da biblioteca.
  - Validação por re-docking do ligante cristalográfico (RMSD < 2 Å).
  - Dinâmica molecular curta e MM-GBSA só para os melhores pares, se o hardware permitir.
- **Interpretação cautelosa:** o linalol e outros monoterpenos agem principalmente sobre a membrana, e o docking é evidência fraca para eles. Isso será complementado com descritores de partição (logP, TPSA, volume) e com a discussão da literatura.

### M5 — Segurança tópica in silico (planejado)

- **Permeação cutânea:** log *K*p pela equação de Potts–Guy. Já está calculado para toda a biblioteca, usando logP estimado (Crippen).
- **Sensibilização e toxicidade:**
  - alertas estruturais por SMARTS (aceptores de Michael, hidroperóxidos, aldeídos);
  - modelos públicos: VEGA (sensibilização cutânea), Pred-Skin, ProTox-3 e Ames.
- **Foco:** os **hidroperóxidos de linalol**, alérgenos de contato reconhecidos que se formam por autoxidação, processo favorecido por calor e luz. Eles já estão na biblioteca. Também entram as metilxantinas e flavan-3-óis, para o perfil do extrato EtOH.
- **Comparação qualitativa:** com a toxicidade em *A. salina* e a ausência de citotoxicidade em HCT116 relatadas por Fonseca Júnior (2021).

### M6 — Espaço químico e similaridade (planejado)

- **Nova biblioteca compartilhada:** filtros UV permitidos (Anexo VI do Reg. CE 1223/2009) e antioxidantes cosméticos de referência.
- **Método:** fingerprints Morgan/ECFP4 e descritores, com PCA e UMAP. Similaridade de Tanimoto dos constituintes com os filtros e antioxidantes comerciais.
- **Produto:** figura de posicionamento para o artigo, com o argumento do "substituto natural" sustentado por distância química.


### M7 — Mecanismos de instabilidade das emulsões (em execução, desde 04/11/2026)

**Motivação:** o Edson está escrevendo um artigo sobre a estabilidade das emulsões. No estudo de estabilidade preliminar (estresse térmico, gela-degela, exposição e proteção à luz), Fonseca Júnior (2021) relatou três alterações qualitativas que o M7 procura explicar no nível molecular, sem reinterpretar os números das Tabelas 12–17:

| Alteração observada (Fonseca Júnior, 2021) | Formulações | Pergunta molecular | Parte |
|---|---|---|---|
| Perda do odor do OE de pau-rosa após estresse térmico e gela-degela | F2, F8, F9, F10 | O linalol escapa para o ar, oxida, ou as duas coisas? | M7a + M2 |
| Mudança de cor com o extrato etanólico | F4, F9, F10 (não F3, com extrato DCM) | Os produtos de oxidação do catecol absorvem no visível? | M7b + M2 |
| Separação de fase em superfície em todas após estresse térmico; ativos aumentam a viscosidade, exceto a quercetina (F6) | todas | Onde cada molécula fica na emulsão e como afeta a interface? | M7a + M7c |

Condições das formulações usadas como entrada: emulsão O/A não iônica (base Croda com Ceteareth-20; fase oleosa com óleo mineral, miristato de isopropila e dimeticona), pH 3,2–4,2. Nessa faixa os fenóis (pKa > 7) e as xantinas estão neutros, então log *D* ≈ log *P*.

**M7a — Partição e volatilidade (✅ concluído em 04/11/2026).**
- Script `modelagem/estabilidade/particao_volatilidade.py` → `particao_volatilidade.csv`.
- Energia livre de solvatação GFN2-xTB/ALPB em água, 1-octanol e n-hexadecano (modelo da fase oleosa), cada fase com confôrmero e geometria relaxados (Ehlert *et al.*, 2021).
- Saídas: log *P* (octanol/água), log *K* (óleo/água), log *K* (ar/água, constante de Henry adimensional) e log *K* (ar/óleo), além da localização provável na emulsão (gota de óleo, interface ou fase aquosa).
- Moléculas (15): linalol, cis-óxido de linalol, α-terpineol, geraniol, limoneno, benzoato de benzila, os dois hidroperóxidos de linalol, cafeína, teobromina, teofilina, epicatequina, catequina, o-quinona da epicatequina e quercetina.
- Validação: log *K*ow experimental de 8 moléculas e constante de Henry do linalol e do α-terpineol (PubChem, consultado em 04/11/2026). O erro médio vai junto com o resultado.
- **Resultado (04/11/2026, `particao_volatilidade.csv`):**
  - *Validação.* A constante de Henry do linalol é bem reproduzida (log *K*ar/água −3,15 calculado; −3,06 experimental); a do α-terpineol erra cerca de 1 unidade log (−2,97 contra −4,03). O log *P* octanol/água do ALPB **não serve**: erro médio de 1,05 e viés de +1,0, chegando a +2 nas xantinas (cafeína 2,01 contra −0,07). Para log *P* usa-se o valor experimental ou o Crippen (erro médio 0,62); do xTB ficam só as tendências entre fases (gás, água, hexadecano).
  - *Terpenoides do OE* (linalol, óxido, geraniol, α-terpineol, benzoato de benzila): ficam na gota de óleo (log *K*óleo/água de 1,7 a 4,0). O linalol escapa da água para o ar com log *K* ≈ −3, mas da gota de óleo com log *K* ≈ −5,7, cerca de 360 vezes menos. Ou seja, a gota de óleo retém o linalol; a perda de odor por evaporação exige que ele saia do óleo, o que o calor do estresse térmico favorece. O limoneno é muito mais volátil (log *K*ar/água −0,65).
  - *Hidroperóxidos de linalol:* log *K*óleo/água ≈ −0,75 e log *K*ar/água ≈ −7. A oxidação desloca a molécula da gota de óleo para a interface e a torna cerca de 10 000 vezes menos volátil que o linalol, ou seja, praticamente sem odor. Assim, **as duas vias (evaporação e autoxidação) são compatíveis com a perda de odor**; o M2 (BDE C–H alílica do linalol) ajuda a pesar a segunda.
  - *Polifenóis* (epicatequina, catequina, quercetina) e *xantinas*: ficam na fase aquosa ou na interface, nunca na gota de óleo. É na fase aquosa (pH 3,2–4,2, com O₂ dissolvido) que a oxidação do catecol e a mudança de cor ocorreriam. A o-quinona é um pouco menos hidrofílica que a epicatequina (log *K*óleo/água −3,8 contra −5,6). Ressalva: o hexadecano penaliza demais os compostos polares, e a fase oleosa real tem miristato de isopropila, mais polar. Por isso a quercetina (log *K*ow exp 2,13) provavelmente fica mais na interface do que o cálculo sugere. A coluna `local_provavel` é uma tendência, não uma localização.

**M7b — Origem da cor (na fila do M3).**
- A o-quinona da epicatequina e a desidrodicatequina A (dímero amarelo de oxidação da catequina; PubChem CID 182270; Guyot *et al.*, 1996) entram no fim da fila do TD-DFT, com o mesmo nível do M3. O espectro foi estendido até 600 nm, com área na faixa visível (400–500 nm).
- Leitura: se os produtos de oxidação absorverem em 400–500 nm e a epicatequina não, a mudança de cor de F4, F9 e F10 é compatível com a oxidação do catecol. O potencial de ionização e a BDE O–H do M2 dizem com que facilidade essa oxidação começa. As xantinas do extrato DCM (F3) não têm catecol, o que é coerente com a ausência de mudança de cor.
- O dímero tem 66 átomos e é o cálculo mais caro da fila (opcional).

**M7c — Compatibilidade com o sistema emulsionante (✅ concluído em 04/11/2026).**
- Script `modelagem/estabilidade/compatibilidade_emulsionante.py` → `compat_validacao.csv`, `compat_ostwald.csv`, `compat_afinidade.csv`.
- Método: EHL de Griffin da Ceteareth-20 calculado da estrutura; UNIFAC original (grupos pelo ugropy, modelo do thermo) a 25 °C, soluto a diluição infinita; volumes molares pelo modelo de Abdulelah-Gani. Os parâmetros de Hansen previstos no plano foram trocados pelo UNIFAC, que é implementado, citável e validável aqui (tabelas de Hansen por contribuição de grupos não estavam disponíveis numa fonte verificável). O EHL requerido da fase oleosa depende de valores tabelados de literatura e não foi estimado.
- **Concentração real dos ativos (leitura do método):** os extratos e o OE foram diluídos em propilenoglicol a 10 mg/mL e essa solução entrou a 10 % (m/m) (Fonseca Júnior, 2021, Tabela 10 e texto do preparo). O OE na emulsão fica, portanto, em torno de **0,1 %**, diluído numa fase lipofílica de cerca de 10–16 % (óleo mineral, miristato de isopropila, dimeticona e a parte oleosa da base). Essa leitura é a premissa das conclusões abaixo.
- **Resultados:**
  - *EHL da Ceteareth-20:* 15,3 (C18) a 15,7 (C16), na faixa de emulsionantes O/A.
  - *Validação do UNIFAC:* log *K*ow de 8 moléculas com erro médio de 0,31 e viés de +0,07, bem melhor que o ALPB do M7a. O UNIFAC Dortmund foi testado e descartado (erro maior nos terpenoides). A solubilidade do hexadecano em água sai cerca de 3 ordens de grandeza alta pelo UNIFAC; onde há dado experimental (PubChem), ele é usado.
  - *Amadurecimento de Ostwald:* os constituintes do OE são de 10⁴ a 10⁸ vezes mais solúveis em água (em fração de volume) que o óleo mineral e o miristato de isopropila; o linalol, 10⁷·⁸ vezes. Uma gota só de OE amadureceria depressa. Mas, com o OE a ~0,1 % dentro de gotas dominadas por óleos quase insolúveis, o efeito composicional de Kabalnov tende a frear esse processo. Junto com a observação de que **a separação de fase ocorreu em todas as formulações, inclusive na F1 (branco)** (Fonseca Júnior, 2021), isso indica que a separação após o estresse térmico é uma propriedade da base nessas condições, e não um efeito dos ativos.
  - *Onde os terpenoides se acomodam:* o linalol, o óxido de linalol, o α-terpineol e o geraniol têm afinidade muito maior pelo álcool cetílico (ln γ∞ ≈ 0–0,2) e pelo segmento de PEG (≈ 0,6–0,8) do que pelo óleo mineral (≈ 2,3) ou pela dimeticona (≈ 1,6–2,0). Ou seja, tendem a se concentrar na camada interfacial (álcool cetoestearílico + Ceteareth-20), e não no núcleo da gota. O limoneno é a exceção e fica no óleo.
  - *Polifenóis × PEG (hipótese do aumento de viscosidade):* **não sustentada pelo UNIFAC.** A epicatequina prefere água, propilenoglicol e glicerina (ln γ∞ negativos) ao segmento de PEG (+1,0). A ressalva é forte: o UNIFAC descreve mal moléculas com muitos OH vizinhos, não distingue estereoisômeros (catequina = epicatequina) e não tem grupos para quercetina, o-quinona, xantinas e hidroperóxidos. A complexação polifenol–PEG por ligação de H, conhecida para taninos, exigiria um cálculo de energia de interação (xTB/DFT ou dinâmica molecular). Por ora a hipótese fica **em aberto**.
- Dinâmica molecular da interface (GROMACS) continua fora do núcleo por custo.

**Limitação:** o estudo de referência é de estabilidade **preliminar**. As análises propõem mecanismos plausíveis, mas não demonstram degradação; isso exigiria CG-EM do óleo ou CLAE dos extratos após o estresse, que estão fora do escopo (seção 4).

---

## 3. Produto esperado

Um artigo complementar in silico sobre as **bases moleculares das atividades antioxidante, fotoprotetora e antimicrobiana de constituintes de *P. cupana* e *A. rosaeodora* para fitocosméticos**. O estudo de Fonseca Júnior (2021) seria citado como referência experimental e a autoria seria combinada com a orientadora.

- **Núcleo mínimo publicável:** M1 + M2 + M3.
- **Complementos:** M4, M5 e M6 entram como seções adicionais ou material suplementar.
- **Artigo de estabilidade do Edson:** o M7 (com M2 e M3) pode entrar como seção de discussão mecanística ou material suplementar, conforme a autoria combinada.

## 4. Fora do escopo

Não fazem parte deste estudo, por se tratar de trabalho de terceiro: correções no texto ou nas tabelas do estudo de referência, reanálise estatística ou reológica dos dados de formulação e novos ensaios de bancada (CG-EM, LC-MS, checkerboard). Se os autores originais quiserem, esses pontos podem ser discutidos à parte.

## 5. Reprodutibilidade

| Item | Local |
|---|---|
| Ambiente Python (RDKit, PySCF, tblite, geomeTRIC, ugropy, thermo) | `QUIMIOINFORMÁTICA/.venv` |
| Código compartilhado (`quimiolib`: LOTUS, PubChem, Europe PMC, descritores, bibliotecas, `qm`) | `QUIMIOINFORMÁTICA/_compartilhado/` |
| Bibliotecas curadas (com `build.py` reprodutível) | `_compartilhado/bibliotecas/` |
| Biblioteca do projeto | `edson/biblioteca/` |
| DFT antioxidante (script, espécies em JSON, tabelas) | `edson/modelagem/dft_antioxidantes/` e `edson/rodar_dft_antioxidantes.sh` |
| TD-DFT UV (inclui M7b) | `edson/modelagem/tddft_uv/` e `edson/rodar_tddft_uv.sh` |
| Partição e volatilidade (M7a) e compatibilidade com o emulsionante (M7c) | `edson/modelagem/estabilidade/` |

## 6. Referências usadas na metodologia

- Ehlert S., Stahn M., Spicher S., Grimme S. (2021). Robust and efficient implicit solvation model for fast semiempirical methods. *J. Chem. Theory Comput.* 17, 4250–4261.
- Fredenslund A., Jones R. L., Prausnitz J. M. (1975). Group-contribution estimation of activity coefficients in nonideal liquid mixtures. *AIChE J.* 21, 1086–1099.
- Galano A., Alvarez-Idaboy J. R. (2013). A computational methodology for accurate predictions of rate constants in solution: application to the assessment of primary antioxidant activity. *J. Comput. Chem.* 34, 2430–2445.
- Guyot S., Vercauteren J., Cheynier V. (1996). Structural determination of colourless and yellow dimers resulting from (+)-catechin coupling catalysed by grape polyphenoloxidase. *Phytochemistry* 42, 1279–1288.
- Kabalnov A. S., Shchukin E. D. (1992). Ostwald ripening theory: applications to fluorocarbon emulsion stability. *Adv. Colloid Interface Sci.* 38, 69–97.
- Potts R. O., Guy R. H. (1992). Predicting skin permeability. *Pharm. Res.* 9, 663–669.
- Rimarčík J., Lukeš V., Klein E., Ilčin M. (2010). *J. Mol. Struct. THEOCHEM* 952, 25. doi:10.1016/j.theochem.2010.04.002
- Rottmannová L., Škorňa P., Rimarčík J., Lukeš V., Klein E. (2013). Solvation enthalpies of the proton in polar and non-polar solvents: theoretical study. *Acta Chim. Slovaca* 6, 60–63. doi:10.2478/acs-2013-0011
- Rutz A. et al. (2022). The LOTUS initiative for open knowledge management in natural products research. *eLife* 11, e70780.
- Sköld M. et al. (2004). *Chem. Res. Toxicol.* 17, 1697. doi:10.1021/tx049831z
- Teles A. M. et al. (2021). *Aniba rosaeodora* (var. *amazonica* Ducke) essential oil: chemical composition, antibacterial, antioxidant and antitrypanosomal activity. *Antibiotics* 10, 24. doi:10.3390/antibiotics10010024
- Xavier J. K. A. M. et al. (2021). Essential oil composition and DNA barcode and identification of *Aniba* species (Lauraceae) growing in the Amazon region. *Molecules* 26, 1914. doi:10.3390/molecules26071914
