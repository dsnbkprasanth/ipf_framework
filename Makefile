# Convenience targets. Activate the environment first (conda activate ipffx).
.PHONY: demo download prepare atlas lincs real gui clean

demo:
	python scripts/run_pipeline.py --demo --out results_demo

download:
	bash scripts/download_data.sh

prepare:
	python scripts/prepare_bulk.py

atlas:
	python scripts/build_atlas_h5ad.py

lincs:
	python scripts/build_lincs.py --out data/lincs

real:
	bash scripts/run_real.sh

gui:
	./run_gui.sh

clean:
	rm -rf results_* __pycache__ ipffx/__pycache__
