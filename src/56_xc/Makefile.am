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
AM_FCFLAGS = @FPPFLAGS@ @FCFLAGS_FREEFORM@ @FCFLAGS_MODDIR@ @fcflags_opt_56_xc@

AM_CPPFLAGS = \
	@src_54_abiutil_incs@ \
	@src_32_util_incs@ \
	@src_28_numeric_noabirule_incs@ \
	@src_44_abitypes_defs_incs@ \
	@src_42_libpaw_incs@ \
	@src_51_manage_mpi_incs@ \
	@src_16_hideleave_incs@ \
	@src_52_fft_mpi_noabirule_incs@ \
	@src_41_xc_lowlevel_incs@ \
	@src_41_geometry_incs@ \
	@src_01_linalg_ext_incs@ \
	@src_18_timing_incs@ \
	@src_12_hide_mpi_incs@ \
	@src_53_spacepar_incs@ \
	@src_11_memory_mpi_incs@ \
	@src_10_defs_incs@ \
	@src_14_hidewrite_incs@ \
	@src_53_ffts_incs@ \
	@src_incs_incs@ \
	@lib_fft_incs@ \
	@lib_libxc_incs@ \
	@fallbacks_incs@ \
	@fc_mod_incs@

# Regular source files
lib56_xc_srcs = \
	calc_smeared_density.F90 \
	gammapositron.F90 \
	gammapositron_fft.F90 \
	hartre.F90 \
	interfaces_56_xc.F90 \
	m_electronpositron.F90 \
	m_kxc.F90 \
	m_xc_vdw.F90 \
	mkcore.F90 \
	mkvxc3.F90 \
	mkvxcgga3.F90 \
	phase.F90 \
	rhohxc.F90 \
	rhohxcpositron.F90 \
	xcpot.F90 \
	xcden.F90 \
	interfaces_56_xc.F90

# Library description
noinst_LIBRARIES = lib56_xc.a

lib56_xc_a_SOURCES= $(lib56_xc_srcs)

CLEANFILES = \
	calc_smeared_density_cpp.f90 \
	gammapositron_cpp.f90 \
	gammapositron_fft_cpp.f90 \
	hartre_cpp.f90 \
	interfaces_56_xc_cpp.f90 \
	m_electronpositron_cpp.f90 \
	m_kxc_cpp.f90 \
	m_xc_vdw_cpp.f90 \
	mkcore_cpp.f90 \
	mkvxc3_cpp.f90 \
	mkvxcgga3_cpp.f90 \
	phase_cpp.f90 \
	rhohxc_cpp.f90 \
	rhohxcpositron_cpp.f90 \
	xcpot_cpp.f90 \
	xcden_cpp.f90 \
	interfaces_56_xc_cpp.f90

EXTRA_DIST = abinit.src

EXTRA_DIST += abinit.dep

#Dependencies (inside the directory) of the directory ./src/56_xc
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
	interfaces_56_xc.$(MODEXT) \
	m_electronpositron.$(MODEXT) \
	m_kxc.$(MODEXT) \
	m_xc_vdw.$(MODEXT) 

gammapositron_fft.$(OBJEXT): interfaces_56_xc.$(OBJEXT) m_electronpositron.$(OBJEXT) 

interfaces_56_xc.$(OBJEXT): m_electronpositron.$(OBJEXT) 

m_kxc.$(OBJEXT): interfaces_56_xc.$(OBJEXT) 

mkvxc3.$(OBJEXT): interfaces_56_xc.$(OBJEXT) 

mkvxcgga3.$(OBJEXT): interfaces_56_xc.$(OBJEXT) 

rhohxc.$(OBJEXT): interfaces_56_xc.$(OBJEXT) m_electronpositron.$(OBJEXT) m_xc_vdw.$(OBJEXT) 

rhohxcpositron.$(OBJEXT): interfaces_56_xc.$(OBJEXT) m_electronpositron.$(OBJEXT) 

xcden.$(OBJEXT): interfaces_56_xc.$(OBJEXT) 

xcpot.$(OBJEXT): interfaces_56_xc.$(OBJEXT) 

EXTRA_DIST += _56_xc_
