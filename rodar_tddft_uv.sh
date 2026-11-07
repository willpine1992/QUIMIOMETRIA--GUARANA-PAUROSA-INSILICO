#!/bin/zsh
# M3 · TD-DFT dos espectros UV (Edson) — iniciar, acompanhar, parar, analisar e agendar.
#
#   ./rodar_tddft_uv.sh iniciar     inicia já (1 molécula por vez, prioridade baixa)
#   ./rodar_tddft_uv.sh agendar     espera o DFT do M2 terminar e então inicia sozinho
#   ./rodar_tddft_uv.sh status      progresso e processos
#   ./rodar_tddft_uv.sh parar       interrompe a fila e/ou o agendamento (o que terminou fica salvo)
#   ./rodar_tddft_uv.sh analise     espectros + resumo (resumo_uv.csv, espectros_uv.csv)
#   ./rodar_tddft_uv.sh vigiar      regenera o painel a cada molécula pronta e a cada 15 min (iniciar já liga isto)
#
# Retomável: cada molécula pronta vira modelagem/tddft_uv/resultados/<id>.json.
# Ordem: benzofenona-3 (calibração) → linalol → óxido de linalol → cafeína → teobromina → teofilina
#        → epicatequina → quercetina → [M7] o-quinona da epicatequina. O dímero (desidrodicatequina A, 66 átomos) roda depois do M2b. Prioridade baixa (nice 15); o Mac não dorme (caffeinate).

set -u
EDSON="${0:A:h}"
DIR="$EDSON/modelagem/tddft_uv"
PY="$EDSON/../../.venv/bin/python"
LOG="$DIR/run.log"
PIDF="$DIR/.pid"
AGENDA_PIDF="$DIR/.pid_agenda"
DFT_PIDF="$EDSON/modelagem/dft_antioxidantes/.pid"
VIGIA_PIDF="$DIR/.pid_vigia"
TOTAL=9   # + 1 opcional (desidrodicatequina A), que roda depois do M2b via rodar_m2b.sh

vivo() { [[ -f "$1" ]] && kill -0 "$(cat "$1")" 2>/dev/null }

case "${1:-status}" in
  iniciar)
    if vivo "$PIDF"; then echo "Já está rodando (PID $(cat "$PIDF"))."; exit 0; fi
    if vivo "$DFT_PIDF"; then echo "⚠️  O DFT do M2 ainda está rodando. Use: $0 agendar"; exit 1; fi
    cd "$DIR"
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONUNBUFFERED=1 PYTHONWARNINGS=ignore \
      nohup caffeinate -i nice -n 15 "$PY" tddft_uv.py >> "$LOG" 2>&1 &
    echo $! > "$PIDF"
    echo "TD-DFT iniciado (PID $!). Log: $LOG"
    "$0" vigiar >/dev/null 2>&1
    ;;
  agendar)
    if vivo "$AGENDA_PIDF"; then echo "Já agendado (PID $(cat "$AGENDA_PIDF"))."; exit 0; fi
    nohup zsh -c '
      while [[ -f "'"$DFT_PIDF"'" ]] && kill -0 "$(cat "'"$DFT_PIDF"'")" 2>/dev/null; do sleep 60; done
      sleep 90   # deixa o vigia do M2 gerar o painel final
      echo "[$(date "+%d/%m %H:%M")] DFT do M2 terminou; iniciando TD-DFT" >> "'"$LOG"'"
      "'"$0"'" iniciar >> "'"$LOG"'" 2>&1
      rm -f "'"$AGENDA_PIDF"'"' >/dev/null 2>&1 &
    echo $! > "$AGENDA_PIDF"
    echo "Agendado (PID $!): o TD-DFT começa sozinho quando o DFT do M2 terminar."
    ;;
  status)
    arqs=("$DIR"/resultados/*.json(N)); n=${#arqs}
    echo "TD-DFT: $n / $TOTAL moléculas prontas"
    if vivo "$PIDF"; then
      echo "Rodando (PID $(cat "$PIDF")) · $(grep '^\[' "$LOG" | tail -1)"
    elif vivo "$AGENDA_PIDF"; then
      echo "Agendado: aguardando o DFT do M2 terminar."
    else
      echo "Não está rodando."
    fi
    grep -E "ok em" "$LOG" 2>/dev/null | tail -3
    grep -E "Traceback|Error" "$LOG" 2>/dev/null | tail -3
    ;;
  parar)
    vivo "$AGENDA_PIDF" && kill "$(cat "$AGENDA_PIDF")" && echo "Agendamento cancelado."
    rm -f "$AGENDA_PIDF"
    if vivo "$PIDF"; then kill "$(cat "$PIDF")"; pkill -f "tddft_uv.py" 2>/dev/null; echo "Fila parada."; fi
    rm -f "$PIDF"
    ;;
  vigiar)
    if vivo "$VIGIA_PIDF"; then exit 0; fi
    nohup zsh -c '
      gerar() { (cd "'"$EDSON"'/modelagem/dft_antioxidantes" && PYTHONWARNINGS=ignore nice -n 15 "'"$PY"'" dashboard.py >/dev/null 2>&1) }
      contar() { ls "'"$DIR"'"/resultados/*.json 2>/dev/null | wc -l | tr -d " " }
      ultimo=$(contar); t=0; gerar
      while [[ -f "'"$PIDF"'" ]] && kill -0 "$(cat "'"$PIDF"'")" 2>/dev/null; do
        sleep 60; t=$((t+1)); n=$(contar)
        if [[ "$n" != "$ultimo" || $t -ge 15 ]]; then gerar; ultimo=$n; t=0; fi
      done
      gerar; rm -f "'"$VIGIA_PIDF"'"' >/dev/null 2>&1 &
    echo $! > "$VIGIA_PIDF"
    ;;
  analise)
    cd "$DIR" && PYTHONWARNINGS=ignore "$PY" tddft_uv.py --analise
    ;;
  *)
    sed -n '2,8p' "$0"; exit 1 ;;
esac
