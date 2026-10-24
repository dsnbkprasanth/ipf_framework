#!/usr/bin/env bash
# Launch the IPF framework GUI. Activate the env first (conda activate ipffx
# OR source .venv/bin/activate), then run: ./run_gui.sh
set -e
cd "$(dirname "$0")"
exec streamlit run app/streamlit_app.py "$@"
