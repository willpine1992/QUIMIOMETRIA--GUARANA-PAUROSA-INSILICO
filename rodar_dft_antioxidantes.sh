#!/bin/zsh
# DFT dos antioxidantes (Edson) — iniciar, acompanhar, parar e analisar.
#
#   ./rodar_dft_antioxidantes.sh iniciar [workers]   inicia em segundo plano (padrão: 3 workers)
#   ./rodar_dft_antioxidantes.sh status              progresso (espécies prontas / total) e processos
#   ./rodar_dft_antioxidantes.sh parar               interrompe (o que já terminou fica salvo)
#   ./rodar_dft_antioxidantes.sh analise             gera as tabelas a partir do que já foi calculado
#   ./rodar_dft_antioxidantes.sh dashboard           regenera o painel HTML (edson/dashboard_dft_antioxidantes.html)
#   ./rodar_dft_antioxidantes.sh vigiar              regenera o painel a cada cálculo concluído (iniciar já liga isto)
#
# O cálculo é RETOMÁVEL: cada espécie pronta vira um JSON em modelagem/dft_antioxidantes/especies/.
# Rodar "iniciar" de novo continua de onde parou.
# Leva de 5 a 8 h num Mac M1 de 8 GB. Ideal: deixar rodando à noite, com o Docker Desktop FECHADO
# (a VM do Docker reserva 4 GB e faz o cálculo travar por falta de memória).
# O Mac não dorme enquanto roda (caffeinate); pode fechar o terminal.
# Roda com prioridade baixa (nice 15): outros programas (fatiador 3D etc.) têm preferência na CPU.

set -u
EDSON="${0:A:h}"
DIR="$EDSON/modelagem/dft_antioxidantes"
PY="$EDSON/../../.venv/bin/python"
LOG="$DIR/run.log"
PIDF="$DIR/.pid"
TOTAL=53   # 52 espécies + átomo de H
DASH="$EDSON/dashboard_dft_antioxidantes.html"
VIGIA_PIDF="$DIR/.pid_vigia"

gerar_dashboard() { (cd "$DIR" && PYTHONWARNINGS=ignore nice -n 15 "$PY" dashboard.py >/dev/null 2>&1) }
contar() { ls "$DIR/especies/dft" 2>/dev/null | wc -l | tr -d ' ' }

rodando() { [[ -f "$PIDF" ]] && kill -0 "$(cat "$PIDF")" 2>/dev/null }

case "${1:-status}" in
  iniciar)
    if rodando; then echo "Já está rodando (PID $(cat "$PIDF")). Use: $0 status"; exit 0; fi
    W="${2:-3}"
    if pgrep -f com.docker.virtualization >/dev/null; then
      echo "⚠️  O Docker Desktop está aberto e reserva 4 GB de RAM; o DFT tende a travar por falta de memória."
      read "r?Continuar mesmo assim? [s/N] "; [[ "$r" == [sS] ]] || exit 1
    fi
    echo "Memória livre: $(memory_pressure | tail -1 | awk '{print $NF}')"
    cd "$DIR"
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONUNBUFFERED=1 PYTHONWARNINGS=ignore \
      nohup caffeinate -i nice -n 15 "$PY" dft_antioxidantes.py --workers "$W" >> "$LOG" 2>&1 &
    echo $! > "$PIDF"
    echo "Iniciado com $W workers (PID $!). Log: $LOG"
    "$0" vigiar >/dev/null 2>&1
    echo "Painel (atualiza sozinho): $DASH"
    echo "Acompanhe com: $0 status"
    ;;
  status)
    n=$(ls "$DIR/especies/dft" 2>/dev/null | wc -l | tr -d ' ')
    echo "DFT: $n / $TOTAL espécies prontas"
    if rodando; then
      echo "Rodando (PID $(cat "$PIDF")):"
      ps -Ao %cpu,etime,command | grep "[P]ython.app" | awk '$1>20 {print "  worker " $1 "% CPU, há " $2}'
      echo "Memória livre: $(memory_pressure | tail -1 | awk '{print $NF}')"
    else
      echo "Não está rodando."
    fi
    grep -E "ERRO|Traceback" "$LOG" 2>/dev/null | tail -3
    if [[ "$n" -ge "$TOTAL" ]]; then echo "✅ Completo — rode: $0 analise"; fi
    ;;
  parar)
    if rodando; then
      kill "$(cat "$PIDF")" 2>/dev/null
      pkill -f "dft_antioxidantes.py --workers" 2>/dev/null
      rm -f "$PIDF"; echo "Parado. O que já terminou está salvo em especies/."
      gerar_dashboard
    else
      pkill -f "dft_antioxidantes.py --workers" 2>/dev/null; echo "Nada rodando."
    fi
    ;;
  dashboard)
    gerar_dashboard && echo "Painel: $DASH"
    ;;
  vigiar)
    if [[ -f "$VIGIA_PIDF" ]] && kill -0 "$(cat "$VIGIA_PIDF")" 2>/dev/null; then exit 0; fi
    nohup zsh -c '
      n0=-1
      while true; do
        n=$(ls "'"$DIR"'/especies/dft" 2>/dev/null | wc -l | tr -d " ")
        if [[ "$n" != "$n0" ]]; then
          (cd "'"$DIR"'" && PYTHONWARNINGS=ignore nice -n 15 "'"$PY"'" dashboard.py >/dev/null 2>&1); n0=$n
        fi
        if ! { [[ -f "'"$PIDF"'" ]] && kill -0 "$(cat "'"$PIDF"'")" 2>/dev/null; }; then
          (cd "'"$DIR"'" && PYTHONWARNINGS=ignore "'"$PY"'" dashboard.py >/dev/null 2>&1); exit 0
        fi
        sleep 30
      done' >/dev/null 2>&1 &
    echo $! > "$VIGIA_PIDF"
    ;;
  analise)
    cd "$DIR" && PYTHONWARNINGS=ignore "$PY" dft_antioxidantes.py --so-analise
    echo
    echo "Tabelas: $DIR/descritores_por_sitio.csv  e  $DIR/resumo_por_molecula.csv  (avisos em avisos.txt)"
    ;;
  *)
    sed -n '2,9p' "$0"; exit 1 ;;
esac
