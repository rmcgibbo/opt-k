all: latex markdown

latex: paper.tex
	pdflatex paper.tex
	bibtex paper
	pdflatex paper.tex

markdown: paper.tex
	pandoc --to=markdown -o README.md.tmp paper.tex
	python postprocess_markdown.py README.md.tmp > README.md
	rm README.md.tmp

clean:
	rm -f paper.bbl paper.blg paper.out paperNotes.bib texput.log paper.aux paper.log