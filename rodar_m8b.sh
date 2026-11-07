#!/bin/zsh
# M8b · refino ωB97X-V/def2-SVP (PCM água + BSSE) dos melhores complexos cafeína–epicatequina.
#
#   ./rodar_m8b.sh agendar    espera a fila do TD-DFT (M3) terminar e então roda refino_dft.py (2 processos),
#                             em paralelo com o M2b (que usa 1 processo); no fim, refaz o painel
#   ./rodar_m8b.sh iniciar    o mesmo, sem esperar
#   ./rodar_m8b.sh status     progresso
#   ./rodar_m8b.sh parar      interrompe (o que terminou fica salvo; é retomável)

set -u
EDSON="${0:A:h}"
M8="$EDSON/modelagem/interacoes"
M2="$EDSON/modelagem/dft_antioxidantes"
M3="$EDSON/modelagem/tddft_uv"
PY="$EDSON/../../.venv/bin/python"
LOG="$M8/refino_dft.log"
PIDF="$M8/.pid_m8b"
ENV=(OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONUNBUFFERED=1 PYTHONWARNINGS=ignore)

vivo() { [[ -f "$1" ]] && kill -0 "$(cat "$1")" 2>/dev/null }

cadeia() {
  echo "[$(date '+%d/%m %H:%M')] iniciando o refino DFT (M8b)" >> "$LOG"
  cd "$M8"
  env $ENV nice -n 15 "$PY" refino_dft.py --workers 2 >> "$LOG" 2>&1
  cd "$M2"
  env $ENV nice -n 15 "$PY" dashboard.py >> "$LOG" 2>&1
  echo "[$(date '+%d/%m %H:%M')] M8b concluído" >> "$LOG"
}

case "${1:-status}" in
  agendar|iniciar)
    if vivo "$PIDF"; then echo "Já agendado/rodando (PID $(cat "$PIDF"))."; exit 0; fi
    nohup caffeinate -i "$0" "_${1}" > /dev/null 2>&1 &
    echo $! > "$PIDF"
    [[ "$1" == agendar ]] && echo "M8b agendado (PID $(cat "$PIDF")): começa quando o TD-DFT terminar." \
                          || echo "M8b iniciado (PID $(cat "$PIDF"))."
    ;;
  _agendar)
    while vivo "$M3/.pid" || vivo "$M3/.pid_agenda"; do sleep 60; done
    sleep 30
    cadeia; rm -f "$PIDF"
    ;;
  _iniciar)
    cadeia; rm -f "$PIDF"
    ;;
  status)
    if vivo "$PIDF"; then
      if vivo "$M3/.pid" || vivo "$M3/.pid_agenda"; then echo "Agendado: aguardando o TD-DFT."; else echo "Rodando (PID $(cat "$PIDF"))."; fi
    else echo "Não está agendado nem rodando."; fi
    echo "Cálculos prontos: $(ls "$M8/resultados/dft" 2>/dev/null | wc -l | tr -d ' ')/23"
    [[ -f "$LOG" ]] && tail -6 "$LOG"; true
    ;;
  parar)
    if vivo "$PIDF"; then pkill -f "rodar_m8b.sh _" 2>/dev/null; kill "$(cat "$PIDF")" 2>/dev/null; rm -f "$PIDF"; echo "Parado."
    else echo "Nada rodando."; fi
    pkill -f refino_dft.py 2>/dev/null
    ;;
  *) echo "uso: $0 agendar|iniciar|status|parar"; exit 1 ;;
esac
