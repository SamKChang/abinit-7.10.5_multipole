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

abidocdir = $(docdir)/developers
nobase_dist_abidoc_DATA = \
  ./FFT_in_parallel.txt \
  ./HM2.document.txt \
  ./HM5.document.txt \
  ./Markdown-HOWTO \
  ./Patches-HOWTO \
  ./README \
  ./adding_input_variables.txt \
  ./checklist.txt \
  ./contributing.html \
  ./contributors.txt \
  ./debug_make_abiauty \
  ./debug_make_parents \
  ./formabinit.css \
  ./preprocessing-options.txt \
  ./programmer_guide.txt \
  ./robodoc.doc.txt \
  ./rules_OMP.txt \
  ./rules_coding.txt \
  ./rules_paral.txt \
  ./use_cpp.txt \
  ./strain_nonlocal/README \
  ./strain_nonlocal/d2term1.fort \
  ./strain_nonlocal/d2term1.out \
  ./strain_nonlocal/format_code.c \
  ./strain_nonlocal/fort_2fix.scr \
  ./strain_nonlocal/fort_fix3.scr \
  ./strain_nonlocal/fort_fixis.scr \
  ./strain_nonlocal/make_met2str.m \
  ./strain_nonlocal/make_metintstr.m \
  ./strain_nonlocal/make_metintstrv.m \
  ./strain_nonlocal/make_metstr3.m

clean-local:
	rm -f doc-build-stamp
