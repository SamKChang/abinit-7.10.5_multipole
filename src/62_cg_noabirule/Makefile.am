#
# Makefile for ABINIT                                      -*- Automake -*-
# Generated by make-makefiles-corelibs on 2015/10/19 06:29:18 +0000

#
# IMPORTANT NOTE
#
# Any manual change to this file will systematically be overwritten.
# Please modify the make-makefiles-corelibs script or its config file instead.
#

AM_CFLAGS = @ABI_CPPFLAGS@
AM_FCFLAGS = @FPPFLAGS@ @FCFLAGS_FREEFORM@ @FCFLAGS_MODDIR@ @fcflags_opt_62_cg_noabirule@

AM_CPPFLAGS = \
	@src_56_xc_incs@ \
	@src_16_hideleave_incs@ \
	@src_32_util_incs@ \
	@src_44_abitypes_defs_incs@ \
	@src_53_spacepar_incs@ \
	@src_11_memory_mpi_incs@ \
	@src_10_defs_incs@ \
	@src_56_recipspace_incs@ \
	@src_53_ffts_incs@ \
	@src_incs_incs@ \
	@fallbacks_incs@ \
	@fc_mod_incs@

# Regular source files
lib62_cg_noabirule_srcs = \
	bracketing.F90 \
	brent.F90 \
	cgpr.F90 \
	dotproduct.F90 \
	frskerker1.F90 \
	frskerker2.F90 \
	ftfvw1.F90 \
	ftfvw2.F90 \
	linmin.F90 \
	interfaces_62_cg_noabirule.F90

# Library description
noinst_LIBRARIES = lib62_cg_noabirule.a

lib62_cg_noabirule_a_SOURCES= $(lib62_cg_noabirule_srcs)

CLEANFILES = \
	bracketing_cpp.f90 \
	brent_cpp.f90 \
	cgpr_cpp.f90 \
	dotproduct_cpp.f90 \
	frskerker1_cpp.f90 \
	frskerker2_cpp.f90 \
	ftfvw1_cpp.f90 \
	ftfvw2_cpp.f90 \
	linmin_cpp.f90 \
	interfaces_62_cg_noabirule_cpp.f90

EXTRA_DIST = abinit.src

EXTRA_DIST += abinit.dep

#Dependencies (inside the directory) of the directory ./src/62_cg_noabirule
#
#COPYRIGHT
#Copyright (C) 2010-2014 ABINIT group
#This file is distributed under the terms of the
#GNU General Public License, see ~abinit/COPYING
#or http://www.gnu.org/copyleft/gpl.txt .
#
#THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
#To do that: config/scripts/abilint --dependencies . .

CLEANFILES += \
	frskerker1.$(MODEXT) \
	frskerker2.$(MODEXT) \
	ftfvw1.$(MODEXT) \
	ftfvw2.$(MODEXT) \
	interfaces_62_cg_noabirule.$(MODEXT) 

brent.$(OBJEXT): interfaces_62_cg_noabirule.$(OBJEXT) 

cgpr.$(OBJEXT): interfaces_62_cg_noabirule.$(OBJEXT) 

frskerker1.$(OBJEXT): interfaces_62_cg_noabirule.$(OBJEXT) 

frskerker2.$(OBJEXT): interfaces_62_cg_noabirule.$(OBJEXT) 

ftfvw2.$(OBJEXT): interfaces_62_cg_noabirule.$(OBJEXT) 

linmin.$(OBJEXT): interfaces_62_cg_noabirule.$(OBJEXT) 

EXTRA_DIST += abinit.amf

EXTRA_DIST += dummy_functions.inc

EXTRA_DIST += _62_cg_noabirule_
