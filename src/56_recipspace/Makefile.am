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
AM_FCFLAGS = @FPPFLAGS@ @FCFLAGS_FREEFORM@ @FCFLAGS_MODDIR@ @fcflags_opt_56_recipspace@

AM_CPPFLAGS = \
	@src_32_util_incs@ \
	@src_28_numeric_noabirule_incs@ \
	@src_44_abitypes_defs_incs@ \
	@src_42_libpaw_incs@ \
	@src_51_manage_mpi_incs@ \
	@src_16_hideleave_incs@ \
	@src_52_fft_mpi_noabirule_incs@ \
	@src_27_toolbox_oop_incs@ \
	@src_41_geometry_incs@ \
	@src_18_timing_incs@ \
	@src_12_hide_mpi_incs@ \
	@src_11_memory_mpi_incs@ \
	@src_10_defs_incs@ \
	@src_14_hidewrite_incs@ \
	@src_53_ffts_incs@ \
	@src_incs_incs@ \
	@fallbacks_incs@ \
	@fc_mod_incs@

# Regular source files
lib56_recipspace_srcs = \
	getcut.F90 \
	getkgrid.F90 \
	getkpgnorm.F90 \
	getmpw.F90 \
	getph.F90 \
	get_full_kgrid.F90 \
	initylmg.F90 \
	irrzg.F90 \
	kpgio.F90 \
	laplacian.F90 \
	listkk.F90 \
	m_ab7_kpoints.F90 \
	m_bz_mesh.F90 \
	m_double_grid.F90 \
	m_gsphere.F90 \
	m_nesting.F90 \
	mkkin.F90 \
	setmqgrid.F90 \
	setsym.F90 \
	smpbz.F90 \
	symkchk.F90 \
	symg.F90 \
	testkgrid.F90 \
	interfaces_56_recipspace.F90

# Library description
noinst_LIBRARIES = lib56_recipspace.a

lib56_recipspace_a_SOURCES= $(lib56_recipspace_srcs)

CLEANFILES = \
	getcut_cpp.f90 \
	getkgrid_cpp.f90 \
	getkpgnorm_cpp.f90 \
	getmpw_cpp.f90 \
	getph_cpp.f90 \
	get_full_kgrid_cpp.f90 \
	initylmg_cpp.f90 \
	irrzg_cpp.f90 \
	kpgio_cpp.f90 \
	laplacian_cpp.f90 \
	listkk_cpp.f90 \
	m_ab7_kpoints_cpp.f90 \
	m_bz_mesh_cpp.f90 \
	m_double_grid_cpp.f90 \
	m_gsphere_cpp.f90 \
	m_nesting_cpp.f90 \
	mkkin_cpp.f90 \
	setmqgrid_cpp.f90 \
	setsym_cpp.f90 \
	smpbz_cpp.f90 \
	symkchk_cpp.f90 \
	symg_cpp.f90 \
	testkgrid_cpp.f90 \
	interfaces_56_recipspace_cpp.f90

EXTRA_DIST = abinit.src

EXTRA_DIST += abinit.dep

#Dependencies (inside the directory) of the directory ./src/56_recipspace
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
	interfaces_56_recipspace.$(MODEXT) \
	m_ab7_kpoints.$(MODEXT) \
	m_bz_mesh.$(MODEXT) \
	m_double_grid.$(MODEXT) \
	m_gsphere.$(MODEXT) \
	m_nesting.$(MODEXT) 

getkgrid.$(OBJEXT): interfaces_56_recipspace.$(OBJEXT) 

m_ab7_kpoints.$(OBJEXT): interfaces_56_recipspace.$(OBJEXT) 

m_bz_mesh.$(OBJEXT): interfaces_56_recipspace.$(OBJEXT) 

m_double_grid.$(OBJEXT): m_bz_mesh.$(OBJEXT) 

m_nesting.$(OBJEXT): m_bz_mesh.$(OBJEXT) 

setsym.$(OBJEXT): interfaces_56_recipspace.$(OBJEXT) 

testkgrid.$(OBJEXT): interfaces_56_recipspace.$(OBJEXT) 

EXTRA_DIST += abinit.amf

EXTRA_DIST += ab7_kpoints.fortran.h

EXTRA_DIST += _56_recipspace_
