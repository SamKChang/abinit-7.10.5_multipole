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
AM_FCFLAGS = @FPPFLAGS@ @FCFLAGS_FREEFORM@ @FCFLAGS_MODDIR@ @fcflags_opt_12_hide_mpi@

AM_CPPFLAGS = \
	@src_11_memory_mpi_incs@ \
	@src_10_defs_incs@ \
	@src_incs_incs@ \
	@lib_mpi_incs@ \
	@fallbacks_incs@ \
	@fc_mod_incs@

# Regular source files
lib12_hide_mpi_srcs = \
	m_xmpi.F90 \
	m_xomp.F90 \
	interfaces_12_hide_mpi.F90

# Library description
noinst_LIBRARIES = lib12_hide_mpi.a

lib12_hide_mpi_a_SOURCES= $(lib12_hide_mpi_srcs)

CLEANFILES = \
	m_xmpi_cpp.f90 \
	m_xomp_cpp.f90 \
	interfaces_12_hide_mpi_cpp.f90

EXTRA_DIST = abinit.src

EXTRA_DIST += abinit.dep

#Dependencies (inside the directory) of the directory ./src/12_hide_mpi
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
	interfaces_12_hide_mpi.$(MODEXT) \
	m_xmpi.$(MODEXT) \
	m_xomp.$(MODEXT) 

EXTRA_DIST += abinit.amf

EXTRA_DIST += \
	xmpi_allgather.finc \
        xmpi_allgatherv.finc \
        xmpi_gather.finc \
        xmpi_gatherv.finc \
        xmpi_bcast.finc \
        xmpi_min.finc \
        xmpi_sum_master.finc \
        xmpi_exch.finc \
        xmpi_recv.finc \
        xmpi_irecv.finc \
        xmpi_sum.finc \
        xmpi_alltoall.finc \
        xmpi_ialltoall.finc \
        xmpi_alltoallv.finc \
        xmpi_ialltoallv.finc \
        xmpi_max.finc \
        xmpi_send.finc \
        xmpi_isend.finc \
        xmpi_scatterv.finc

EXTRA_DIST += _12_hide_mpi_
