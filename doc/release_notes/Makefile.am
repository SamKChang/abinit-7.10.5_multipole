#                                                          -*- Automake -*-
# Makefile for the Abinit Documentation
# Generated by make-makefiles-doc on 2015/10/19 06:29:11 +0000

#
# IMPORTANT NOTE
#
# Any manual change to this file will systematically be overwritten.
# Please modify the make-makefiles-doc script or its config file instead.
#

#
# Suffixes
#

SUFFIXES = .html .pdf .tex .gv .md

# --------------------------------------------------------------------------- #

#
# Rules to build Graphviz documents
#

if DO_BUILD_GV
.gv.pdf:
	$(GRAPHVIZ_DOT) $(GVDFLAGS) -T pdf -o $@ $<
else
.gv.pdf:
	@echo "Please install or upgrade Graphviz to build $@"
endif

# --------------------------------------------------------------------------- #

#
# Rules to build LaTeX documents
#

if DO_BUILD_TEX
.tex.pdf:
	TEXINPUTS="$(abs_builddir):$(abs_srcdir):" \
	  $(LATEX) $(TEXFLAGS) $(TEXFLAGS_EXTRA) $<
	TEXINPUTS="$(abs_builddir):$(abs_srcdir):" \
	  $(LATEX) $(TEXFLAGS) $(TEXFLAGS_EXTRA) $<
else
.tex.pdf:
	@echo "Please install LaTeX to build $@"
endif

# --------------------------------------------------------------------------- #

#
# Rules to build MarkDown documents
#

.md.html:
	$(MARKDOWN) $(MDFLAGS) $< >$@

# --------------------------------------------------------------------------- #


all-local: doc-build-stamp

if DO_BUILD_TEX
doc-build-stamp:
	touch doc-build-stamp
else
doc-build-stamp:
	cd .. && make built-docs-stamp @SET_MAKE@
endif

abidocdir = $(docdir)/release_notes
nobase_dist_abidoc_DATA = \
  ./release_notes_v2.3.htm \
  ./release_notes_v3.0.htm \
  ./release_notes_v3.1.htm \
  ./release_notes_v3.2.htm \
  ./release_notes_v3.3.htm \
  ./release_notes_v3.4.htm \
  ./release_notes_v4.0.htm \
  ./release_notes_v4.1.htm \
  ./release_notes_v4.2.htm \
  ./release_notes_v4.3.html \
  ./release_notes_v4.4.html \
  ./release_notes_v4.5.html \
  ./release_notes_v4.6.html \
  ./release_notes_v4.7.html \
  ./release_notes_v5.0.html \
  ./release_notes_v5.1.html \
  ./release_notes_v5.2.html \
  ./release_notes_v5.3.html \
  ./release_notes_v5.4.html \
  ./release_notes_v5.5.html \
  ./release_notes_v5.6.html \
  ./release_notes_v5.7.html \
  ./release_notes_v5.8.html \
  ./release_notes_v6.0.html \
  ./release_notes_v6.10.html \
  ./release_notes_v6.11.html \
  ./release_notes_v6.12.html \
  ./release_notes_v6.13.html \
  ./release_notes_v6.14.html \
  ./release_notes_v6.15.html \
  ./release_notes_v6.2.html \
  ./release_notes_v6.4.html \
  ./release_notes_v6.5.html \
  ./release_notes_v6.6.html \
  ./release_notes_v6.7.html \
  ./release_notes_v6.8.html \
  ./release_notes_v6.9.html \
  ./release_notes_v7.0.html \
  ./release_notes_v7.1.html \
  ./release_notes_v7.10.html \
  ./release_notes_v7.11.html \
  ./release_notes_v7.2.html \
  ./release_notes_v7.3.html \
  ./release_notes_v7.4.html \
  ./release_notes_v7.5.html \
  ./release_notes_v7.6.html \
  ./release_notes_v7.7.html \
  ./release_notes_v7.8.html \
  ./release_notes_v7.9.html

clean-local:
	rm -f doc-build-stamp
