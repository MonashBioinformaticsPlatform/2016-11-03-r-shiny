
# Compile .Rmd files to create these html files
HTMLS=index.html tutorial.html slides.html

# Create stripped down versions of corresponding .Rmd files
RS=tutorial.R

all : $(RS) $(HTMLS)

%.html : %.Rmd
	Rscript -e 'rmarkdown::render("$<", "all")'

%.R : %.Rmd purify.py
	python purify.py <$< >$@

clean :
	rm $(HTMLS) $(RS)
