.PHONY: all clean

all: paper/paper.pdf

# ── Preprocessing ──────────────────────────────────────────────────────────────
temp/panel.csv: input/california_prop99.csv code/preprocess.py
	python code/preprocess.py

# ── Analysis: produces every artifact consumed by the paper ───────────────────
output/tables/main_result.tex output/figures/sdid_trends.pdf: \
		temp/panel.csv code/analysis.py
	python code/analysis.py

# ── Paper compilation (run pdflatex twice for cross-references) ───────────────
paper/paper.pdf: paper/paper.tex paper/references.bib \
		output/tables/main_result.tex output/figures/sdid_trends.pdf
	cd paper && pdflatex -interaction=nonstopmode paper.tex \
	         && bibtex paper \
	         && pdflatex -interaction=nonstopmode paper.tex \
	         && pdflatex -interaction=nonstopmode paper.tex

# ── Clean ─────────────────────────────────────────────────────────────────────
clean:
	rm -f temp/panel.csv
	rm -f output/tables/main_result.tex
	rm -f output/figures/sdid_trends.pdf output/figures/sdid_trends.png
	rm -f paper/paper.pdf paper/paper.aux paper/paper.log \
	      paper/paper.bbl paper/paper.blg paper/paper.out paper/paper.toc
