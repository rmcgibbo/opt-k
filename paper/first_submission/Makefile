all: paper2.pdf README.md

paper2.pdf: paper2.tex bibliography.bib
	pdflatex paper2.tex
	bibtex paper
	pdflatex paper2.tex

README.md: paper2.tex postprocess_markdown.py
	pandoc --no-wrap --to=markdown -o README.md.tmp paper.tex
	python postprocess_markdown.py README.md.tmp > README.md
	rm README.md.tmp

clean:
	rm -f paper.bbl paper.blg paper.out paperNotes.bib texput.log paper.aux paper.log
	rm -f paper2.bbl paper2.blg paper2.out paperNotes2.bib texput2.log paper2.aux paper2.log

watch:
	echo "Watching .tex files for changes..."
	watchmedo shell-command --command='make all' -w -p='*.tex'
