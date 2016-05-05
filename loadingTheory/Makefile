#Simple Makefile to build docs from latex


DVIFLAGS=-output-format=dvi -src-specials -interaction=errorstopmode 

PDFFLAGS=-output-format=pdf -interaction=nonstopmode

MASTER=Theoryloading

DVIOUT=${MASTER}.dvi
PDFOUT=${MASTER}.pdf

#all: dvi
all:pdf 
#all: dvi pdf

dvi:${DVIOUT}

pdf:${PDFOUT}

${DVIOUT}:${MASTER}.tex
	latex ${DVIFLAGS} ${MASTER}.tex

${PDFOUT}:${MASTER}.tex
	latex ${PDFFLAGS} ${MASTER}.tex

clean:
	rm -f *aux ${DVIOUT} ${PDFOUT} ${MASTER}.log #${MASTER}.toc ${MASTER}.lof ${MASTER}.bbl ${MASTER}.blg
