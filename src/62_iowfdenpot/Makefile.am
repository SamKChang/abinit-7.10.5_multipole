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
AM_FCFLAGS = @FPPFLAGS@ @FCFLAGS_FREEFORM@ @FCFLAGS_MODDIR@ @fcflags_opt_62_iowfdenpot@

AM_CPPFLAGS = \
	@src_32_util_incs@ \
	@src_28_numeric_noabirule_incs@ \
	@src_44_abitypes_defs_incs@ \
	@src_42_libpaw_incs@ \
	@src_51_manage_mpi_incs@ \
	@src_43_wvl_wrappers_incs@ \
	@src_16_hideleave_incs@ \
	@src_56_recipspace_incs@ \
	@src_52_fft_mpi_noabirule_incs@ \
	@src_49_gw_toolbox_oop_incs@ \
	@src_56_io_mpi_incs@ \
	@src_18_timing_incs@ \
	@src_12_hide_mpi_incs@ \
	@src_11_memory_mpi_incs@ \
	@src_10_defs_incs@ \
	@src_14_hidewrite_incs@ \
	@src_41_geometry_incs@ \
	@src_27_toolbox_oop_incs@ \
	@src_incs_incs@ \
	@lib_bigdft_incs@ \
	@lib_etsf_io_incs@ \
	@lib_mpi_incs@ \
	@lib_netcdf_incs@ \
	@fallbacks_incs@ \
	@fc_mod_incs@

# Regular source files
lib62_iowfdenpot_srcs = \
	calcdensph.F90 \
	initwf.F90 \
	ioarr.F90 \
	m_io_screening.F90 \
	m_wfk.F90 \
	out1dm.F90 \
	outwant.F90 \
	randac.F90 \
	rdnpw.F90 \
	WffReadEigK.F90 \
	WffReadSkipK.F90 \
	interfaces_62_iowfdenpot.F90

# Library description
noinst_LIBRARIES = lib62_iowfdenpot.a

lib62_iowfdenpot_a_SOURCES= $(lib62_iowfdenpot_srcs)

CLEANFILES = \
	calcdensph_cpp.f90 \
	initwf_cpp.f90 \
	ioarr_cpp.f90 \
	m_io_screening_cpp.f90 \
	m_wfk_cpp.f90 \
	out1dm_cpp.f90 \
	outwant_cpp.f90 \
	randac_cpp.f90 \
	rdnpw_cpp.f90 \
	WffReadEigK_cpp.f90 \
	WffReadSkipK_cpp.f90 \
	interfaces_62_iowfdenpot_cpp.f90

EXTRA_DIST = abinit.src

EXTRA_DIST += abinit.dep

#Dependencies (inside the directory) of the directory ./src/62_iowfdenpot
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
	interfaces_62_iowfdenpot.$(MODEXT) \
	m_io_screening.$(MODEXT) \
	m_wfk.$(MODEXT) 

initwf.$(OBJEXT): interfaces_62_iowfdenpot.$(OBJEXT) m_wfk.$(OBJEXT) 

EXTRA_DIST += _62_iowfdenpot_
