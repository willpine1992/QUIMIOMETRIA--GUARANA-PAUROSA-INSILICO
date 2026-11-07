#!/bin/zsh
# Troca única (05/11/2026): espera o óxido de linalol terminar no modo direto, para a fila, guarda a
# benzofenona-3 antiga (não convergiu em todos os estados) e reinicia com ajuste de densidade.
DIR="${0:A:h}"; EDSON="$DIR/../.."
while [[ ! -f "$DIR/resultados/oxido_linalol.json" ]]; do
  kill -0 "$(cat "$DIR/.pid" 2>/dev/null)" 2>/dev/null || { echo "[$(date '+%d/%m %H:%M')] fila parou antes do óxido terminar; troca feita mesmo assim" >> "$DIR/run.log"; break; }
  sleep 30
done
"$EDSON/rodar_tddft_uv.sh" parar >> "$DIR/run.log" 2>&1
sleep 5
mkdir -p "$DIR/resultados_direto_antigo"
mv "$DIR/resultados/benzofenona_3.json" "$DIR/resultados_direto_antigo/" 2>/dev/null
mv "$DIR/resultados/orbitais/benzofenona_3.npz" "$DIR/resultados_direto_antigo/" 2>/dev/null
echo "[$(date '+%d/%m %H:%M')] troca para ajuste de densidade; benzofenona-3 será refeita" >> "$DIR/run.log"
"$EDSON/rodar_tddft_uv.sh" iniciar >> "$DIR/run.log" 2>&1
rm -f "$DIR/.pid_troca"
