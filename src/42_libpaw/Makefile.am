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
AM_FCFLAGS = @FPPFLAGS@ @FCFLAGS_FREEFORM@ @FCFLAGS_MODDIR@ @fcflags_opt_42_libpaw@

AM_CPPFLAGS = \
	@src_41_xc_lowlevel_incs@ \
	@src_12_hide_mpi_incs@ \
	@src_10_defs_incs@ \
	@src_11_memory_mpi_incs@ \
	@src_16_hideleave_incs@ \
	@src_14_hidewrite_incs@ \
	@src_incs_incs@ \
	@lib_fox_incs@ \
	@fallbacks_incs@ \
	@fc_mod_incs@

# Regular source files
lib42_libpaw_srcs = \
	m_paral_atom.F90 \
	m_sphharm.F90 \
	m_pawio.F90 \
	m_pawang.F90 \
	m_pawrad.F90 \
	m_pawtab.F90 \
	m_paw_an.F90 \
	m_paw_ij.F90 \
	m_pawfgrtab.F90 \
	m_paw_numeric.F90 \
	m_paw_finegrid.F90 \
	m_pawcprj.F90 \
	m_pawrhoij.F90 \
	m_pawdij.F90 \
	m_pawxc.F90 \
	m_atompaw.F90 \
	m_gaussfit.F90 \
	m_pawxmlps.F90 \
	m_pawpsp.F90 \
	interfaces_42_libpaw.F90

# Library description
noinst_LIBRARIES = lib42_libpaw.a

lib42_libpaw_a_SOURCES= $(lib42_libpaw_srcs)

CLEANFILES = \
	m_paral_atom_cpp.f90 \
	m_sphharm_cpp.f90 \
	m_pawio_cpp.f90 \
	m_pawang_cpp.f90 \
	m_pawrad_cpp.f90 \
	m_pawtab_cpp.f90 \
	m_paw_an_cpp.f90 \
	m_paw_ij_cpp.f90 \
	m_pawfgrtab_cpp.f90 \
	m_paw_numeric_cpp.f90 \
	m_paw_finegrid_cpp.f90 \
	m_pawcprj_cpp.f90 \
	m_pawrhoij_cpp.f90 \
	m_pawdij_cpp.f90 \
	m_pawxc_cpp.f90 \
	m_atompaw_cpp.f90 \
	m_gaussfit_cpp.f90 \
	m_pawxmlps_cpp.f90 \
	m_pawpsp_cpp.f90 \
	interfaces_42_libpaw_cpp.f90

EXTRA_DIST = abinit.src

EXTRA_DIST += abinit.dep

#Dependencies (inside the directory) of the directory ./src/42_libpaw
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
	interfaces_42_libpaw.$(MODEXT) \
	m_atompaw.$(MODEXT) \
	m_gaussfit.$(MODEXT) \
	m_paral_atom.$(MODEXT) \
	m_paw_an.$(MODEXT) \
	m_paw_finegrid.$(MODEXT) \
	m_paw_ij.$(MODEXT) \
	m_paw_numeric.$(MODEXT) \
	m_pawang.$(MODEXT) \
	m_pawcprj.$(MODEXT) \
	m_pawdij.$(MODEXT) \
	m_pawfgrtab.$(MODEXT) \
	m_pawio.$(MODEXT) \
	m_pawpsp.$(MODEXT) \
	m_pawrad.$(MODEXT) \
	m_pawrhoij.$(MODEXT) \
	m_pawtab.$(MODEXT) \
	m_pawxc.$(MODEXT) \
	m_pawxmlps.$(MODEXT) \
	m_sphharm.$(MODEXT) 

m_atompaw.$(OBJEXT): m_paw_numeric.$(OBJEXT) m_pawrad.$(OBJEXT) m_pawtab.$(OBJEXT) 

m_gaussfit.$(OBJEXT): m_paw_numeric.$(OBJEXT) m_pawrad.$(OBJEXT) 

m_paw_an.$(OBJEXT): m_paral_atom.$(OBJEXT) m_pawang.$(OBJEXT) m_pawtab.$(OBJEXT) 

m_paw_finegrid.$(OBJEXT): m_paw_numeric.$(OBJEXT) m_pawtab.$(OBJEXT) m_sphharm.$(OBJEXT) 

m_paw_ij.$(OBJEXT): m_paral_atom.$(OBJEXT) m_pawio.$(OBJEXT) m_pawtab.$(OBJEXT) 

m_pawang.$(OBJEXT): m_sphharm.$(OBJEXT) 

m_pawcprj.$(OBJEXT): m_pawtab.$(OBJEXT) 

m_pawdij.$(OBJEXT): m_paral_atom.$(OBJEXT) m_paw_an.$(OBJEXT) m_paw_finegrid.$(OBJEXT) m_paw_ij.$(OBJEXT) m_pawang.$(OBJEXT) m_pawfgrtab.$(OBJEXT) m_pawio.$(OBJEXT) m_pawrad.$(OBJEXT) m_pawrhoij.$(OBJEXT) m_pawtab.$(OBJEXT) 

m_pawfgrtab.$(OBJEXT): m_paral_atom.$(OBJEXT) 

m_pawpsp.$(OBJEXT): m_atompaw.$(OBJEXT) m_gaussfit.$(OBJEXT) m_paw_numeric.$(OBJEXT) m_pawang.$(OBJEXT) m_pawrad.$(OBJEXT) m_pawtab.$(OBJEXT) m_pawxc.$(OBJEXT) m_pawxmlps.$(OBJEXT) 

m_pawrhoij.$(OBJEXT): m_paral_atom.$(OBJEXT) m_pawang.$(OBJEXT) m_pawio.$(OBJEXT) m_pawtab.$(OBJEXT) 

m_pawxc.$(OBJEXT): m_pawang.$(OBJEXT) m_pawrad.$(OBJEXT) 

m_pawxmlps.$(OBJEXT): m_paw_numeric.$(OBJEXT) m_pawrad.$(OBJEXT) 

EXTRA_DIST += _42_libpaw_
