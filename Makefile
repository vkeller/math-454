#!/bin/bash

pdf:
	pdflatex cours_vince.tex

clean:
	rm *.log *.aux *.nav *.out *.pdf *.out *.snm *.toc *.vrb
