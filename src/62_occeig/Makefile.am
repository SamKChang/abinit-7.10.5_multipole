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
AM_FCFLAGS = @FPPFLAGS@ @FCFLAGS_FREEFORM@ @FCFLAGS_MODDIR@ @fcflags_opt_62_occeig@

AM_CPPFLAGS = \
	@src_32_util_incs@ \
	@src_28_numeric_noabirule_incs@ \
	@src_44_abitypes_defs_incs@ \
	@src_16_hideleave_incs@ \
	@src_56_io_mpi_incs@ \
	@src_56_recipspace_incs@ \
	@src_41_geometry_incs@ \
	@src_18_timing_incs@ \
	@src_12_hide_mpi_incs@ \
	@src_11_memory_mpi_incs@ \
	@src_10_defs_incs@ \
	@src_14_hidewrite_incs@ \
	@src_27_toolbox_oop_incs@ \
	@src_incs_incs@ \
	@lib_bigdft_incs@ \
	@lib_etsf_io_incs@ \
	@lib_netcdf_incs@ \
	@fallbacks_incs@ \
	@fc_mod_incs@

# Regular source files
lib62_occeig_srcs = \
	dos_hdr_write.F90 \
	dos_degeneratewfs.F90 \
	gaus_dos.F90 \
	getnel.F90 \
	get_dos_1band.F90 \
	get_dos_1band_m.F90 \
	init_bess_spl.F90 \
	init_occ_ent.F90 \
	m_ebands.F90 \
	newocc.F90 \
	occeig.F90 \
	pareigocc.F90 \
	printbxsf.F90 \
	prtbltztrp_out.F90 \
	prtbltztrp_tau_out.F90 \
	recip_ylm.F90 \
	sphericaldens.F90 \
	tetrahedron.F90 \
	interfaces_62_occeig.F90

# Library description
noinst_LIBRARIES = lib62_occeig.a

lib62_occeig_a_SOURCES= $(lib62_occeig_srcs)

CLEANFILES = \
	dos_hdr_write_cpp.f90 \
	dos_degeneratewfs_cpp.f90 \
	gaus_dos_cpp.f90 \
	getnel_cpp.f90 \
	get_dos_1band_cpp.f90 \
	get_dos_1band_m_cpp.f90 \
	init_bess_spl_cpp.f90 \
	init_occ_ent_cpp.f90 \
	m_ebands_cpp.f90 \
	newocc_cpp.f90 \
	occeig_cpp.f90 \
	pareigocc_cpp.f90 \
	printbxsf_cpp.f90 \
	prtbltztrp_out_cpp.f90 \
	prtbltztrp_tau_out_cpp.f90 \
	recip_ylm_cpp.f90 \
	sphericaldens_cpp.f90 \
	tetrahedron_cpp.f90 \
	interfaces_62_occeig_cpp.f90

EXTRA_DIST = abinit.src

EXTRA_DIST += abinit.dep

#Dependencies (inside the directory) of the directory ./src/62_occeig
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
	interfaces_62_occeig.$(MODEXT) \
	m_ebands.$(MODEXT) 

gaus_dos.$(OBJEXT): interfaces_62_occeig.$(OBJEXT) 

getnel.$(OBJEXT): interfaces_62_occeig.$(OBJEXT) 

m_ebands.$(OBJEXT): interfaces_62_occeig.$(OBJEXT) 

newocc.$(OBJEXT): interfaces_62_occeig.$(OBJEXT) 

tetrahedron.$(OBJEXT): interfaces_62_occeig.$(OBJEXT) 

EXTRA_DIST += _62_occeig_
