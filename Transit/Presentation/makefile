## Sweave Document MAKE Compiler ##

# for Mac OSX and Linux

# Josephine D. Kressner
# Georgia Institute of Technology
# josiekressner@gatech.edu

# Edits and annotation by Greg Macfarlane
# gregmacfarlane@gatech.edu

# Enter the name of the final .pdf document here:
MASTER = Presentation.pdf
# Your principal .rnw or .tex file should have
# the same stem. If your principal file is
# foo.rnw, set MASTER = foo.pdf

# Enter the name of the source data file here:
#DATAFILE = /Volumes/CEE_LGarrow_Group\$/common/Equifax/Studies/HomePrices.Rdata

# If the source data is in another directory, 
# change the VPATH variable below to the correct
# directory.
# VPATH = /Volumes/CEE_LGarrow_Group\$/common/Equifax/Studies

# This file is designed for UNIX- type
# operating systems and successfully compiles
# on Mac OS X and Red Hat Linux at least. Compiling on
# Windows will require different LaTeX output
# suppressors (>/dev/null will not be recognized) at 
# a minumum.



RNWFILES = $(wildcard *.rnw) 
TEXFILES = $(wildcard *.tex)
DEPENDS = $(patsubst %.rnw,%.tex,$(RNWFILES)) $(TEXFILES) 

all: $(MASTER) 
	@ make clean

#$(MASTER): $(DEPENDS) 

%.tex: %.rnw
	@echo + Sweaving $@ from $< ... 
	knit $< 	

%.pdf: %.tex 
	@ echo + Writing $@ from $< ...
	texi2pdf $< >/dev/null
 
clean:
	@ echo + Cleaning ...
	@ rm -f *.aux *.lof *.log *.lot *.toc Rplots.pdf 
	@ rm -f *.bbl *.blg *.dvi *.Rout
	@ rm -f *.nav *.out *.snm 
#	@ rm -f $(patsubst %.rnw,%.tex,$(RNWFILES)) 

realclean: clean
	@echo + Really cleaning ...
	rm -f $(MASTER)


Stangle:$(RNWFILES)
	R64 CMD STANGLE '$<'
	@ open $< -a RStudio 

menu:
	@ echo + ==============================
	@ echo + .......GNU Make menu..........
	@ echo + data: ............ create data
	@ echo + all: ......... create document
	@ echo + clean: ...... delete aux files
	@ echo + realclean: . delete all output
	@ echo + Stangle: extract R code into R
	@ echo +          
	@ echo + Georgia  | Civil
	@ echo +    Tech /| Engineering
	@ echo + ==============================
