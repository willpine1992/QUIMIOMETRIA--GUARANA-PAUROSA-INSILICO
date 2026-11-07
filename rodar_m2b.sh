#!/bin/zsh
# M2b · refino das geometrias com freq. imaginária + funções de onda dos radicais-chave; depois o dímero do M3.
#
#   ./rodar_m2b.sh agendar    espera o DFT do M2 e o TD-DFT do M3 terminarem e então roda, em sequência:
#                               1) m2b_refino_densidades.py (etapas A e B)
#                               2) análise do M2 (BDE/PA/IP com as geometrias refinadas), superfícies 3D
#                                  (visual_assets.py: HOMO, LUMO, MEP do M3 e spin do M2b) e painel
#                               3) TD-DFT opcional do M3 (desidrodicatequina A, 66 átomos; o mais longo)
#   ./rodar_m2b.sh iniciar    o mesmo, sem esperar (recusa se o M2 ou o M3 estiverem rodando)
#   ./rodar_m2b.sh status     progresso (log)
#   ./rodar_m2b.sh parar      interrompe (o que terminou fica salvo; é retomável)

set -u
EDSON="${0:A:h}"
M2="$EDSON/modelagem/dft_antioxidantes"
M3="$EDSON/modelagem/tddft_uv"
PY="$EDSON/../../.venv/bin/python"
LOG="$M2/m2b.log"
PIDF="$M2/.pid_m2b"
ENV=(OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONUNBUFFERED=1 PYTHONWARNINGS=ignore)

vivo() { [[ -f "$1" ]] && kill -0 "$(cat "$1")" 2>/dev/null }
ocupado() { vivo "$M2/.pid" || vivo "$M3/.pid_agenda" || vivo "$M3/.pid" }

cadeia() {
  cd "$M2"
  env $ENV nice -n 15 "$PY" m2b_refino_densidades.py >> "$LOG" 2>&1
  env $ENV nice -n 15 "$PY" -c "import dft_antioxidantes as d; d.analisar()" >> "$LOG" 2>&1
  echo "[$(date '+%d/%m %H:%M')] gerando superfícies 3D (orbitais, MEP, spin) para o painel" >> "$LOG"
  env $ENV nice -n 15 "$PY" visual_assets.py >> "$LOG" 2>&1
  env $ENV nice -n 15 "$PY" dashboard.py >> "$LOG" 2>&1
  echo "[$(date '+%d/%m %H:%M')] iniciando o TD-DFT opcional (dímero)" >> "$LOG"
  cd "$M3"
  env $ENV nice -n 15 "$PY" tddft_uv.py --opcionais >> "$M3/run.log" 2>&1
  cd "$M2"
  env $ENV nice -n 15 "$PY" visual_assets.py >> "$LOG" 2>&1
  env $ENV nice -n 15 "$PY" dashboard.py >> "$LOG" 2>&1
  echo "[$(date '+%d/%m %H:%M')] fila noturna concluída" >> "$LOG"
}

case "${1:-status}" in
  agendar|iniciar)
    if vivo "$PIDF"; then echo "Já agendado/rodando (PID $(cat "$PIDF"))."; exit 0; fi
    if [[ "$1" == iniciar ]] && ocupado; then echo "⚠️  M2 ou M3 ainda rodando. Use: $0 agendar"; exit 1; fi
    nohup caffeinate -i "$0" _executar > /dev/null 2>&1 &
    echo $! > "$PIDF"
    echo "M2b agendado (PID $(cat "$PIDF")): começa sozinho depois do M2 e do M3."
    ;;
  _executar)
    while ocupado; do sleep 60; done
    sleep 120   # deixa o vigia e a análise do M3 fecharem
    echo "[$(date '+%d/%m %H:%M')] M2 e M3 terminaram; iniciando M2b" >> "$LOG"
    cadeia
    rm -f "$PIDF"
    ;;
  status)
    if vivo "$PIDF"; then
      if ocupado; then echo "Agendado: aguardando M2/M3."; else echo "Rodando (PID $(cat "$PIDF"))."; fi
    else echo "Não está agendado nem rodando."; fi
    [[ -f "$LOG" ]] && tail -6 "$LOG"
    ;;
  parar)
    if vivo "$PIDF"; then pkill -f "rodar_m2b.sh _executar" 2>/dev/null; kill "$(cat "$PIDF")" 2>/dev/null; rm -f "$PIDF"; echo "Parado."
    else echo "Nada rodando."; fi
    pkill -f m2b_refino_densidades.py 2>/dev/null; pkill -f "tddft_uv.py --opcionais" 2>/dev/null
    ;;
  *) echo "uso: $0 agendar|iniciar|status|parar"; exit 1 ;;
esac
