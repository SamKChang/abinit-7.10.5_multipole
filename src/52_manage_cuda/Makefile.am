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
AM_FCFLAGS = @FPPFLAGS@ @FCFLAGS_FREEFORM@ @FCFLAGS_MODDIR@ @fcflags_opt_52_manage_cuda@

AM_CPPFLAGS = \
	@src_12_hide_mpi_incs@ \
	@src_11_memory_mpi_incs@ \
	@src_10_defs_incs@ \
	@src_44_abitypes_defs_incs@ \
	@src_15_gpu_toolbox_incs@ \
	@src_incs_incs@ \
	@lib_gpu_incs@ \
	@fallbacks_incs@ \
	@fc_mod_incs@

# Regular source files
lib52_manage_cuda_srcs = \
	alloc_hamilt_gpu.F90 \
	cuda_rec.cu \
	cuda_rec_cut.cu \
	density.cu \
	gpu_apply_local_potential.cu \
	gpu_compute_nl_hamiltonian.cu \
	gpu_compute_nl_projections.cu \
	gpu_density_accumulation.cu \
	gpu_fourwf.cu \
	gpu_mkkpg.cu \
	gpu_nonlop.cu \
	gpu_sphere.cu \
	prt_utils_rec.cu \
	rec_dens_calc.cu \
	rec_kernels.cu \
	recursion_bth.cu \
	recursion_cut_bth.cu \
	recursion_cut_no_bth.cu \
	recursion_no_bth.cu \
	interfaces_52_manage_cuda.F90

# Library description
noinst_LIBRARIES = lib52_manage_cuda.a

lib52_manage_cuda_a_SOURCES= $(lib52_manage_cuda_srcs)

CLEANFILES = \
	alloc_hamilt_gpu_cpp.f90 \
	cuda_rec_cpp.cu \
	cuda_rec_cut_cpp.cu \
	density_cpp.cu \
	gpu_apply_local_potential_cpp.cu \
	gpu_compute_nl_hamiltonian_cpp.cu \
	gpu_compute_nl_projections_cpp.cu \
	gpu_density_accumulation_cpp.cu \
	gpu_fourwf_cpp.cu \
	gpu_mkkpg_cpp.cu \
	gpu_nonlop_cpp.cu \
	gpu_sphere_cpp.cu \
	prt_utils_rec_cpp.cu \
	rec_dens_calc_cpp.cu \
	rec_kernels_cpp.cu \
	recursion_bth_cpp.cu \
	recursion_cut_bth_cpp.cu \
	recursion_cut_no_bth_cpp.cu \
	recursion_no_bth_cpp.cu \
	interfaces_52_manage_cuda_cpp.f90

EXTRA_DIST = abinit.src

EXTRA_DIST += abinit.dep

#Dependencies (inside the directory) of the directory ./src/52_manage_cuda
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
	interfaces_52_manage_cuda.$(MODEXT) 

SUFFIXES = .o .cu

.cu.o:
	$(NVCC) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(NVCC_CPPFLAGS) $(NVCC_CFLAGS) -c $<

cuda_rec.$(OBJEXT): cuda_rec.cu

cuda_rec_cut.$(OBJEXT): cuda_rec_cut.cu

density.$(OBJEXT): density.cu

gpu_apply_local_potential.$(OBJEXT): gpu_apply_local_potential.cu

gpu_compute_nl_hamiltonian.$(OBJEXT): gpu_compute_nl_hamiltonian.cu

gpu_compute_nl_projections.$(OBJEXT): gpu_compute_nl_projections.cu

gpu_density_accumulation.$(OBJEXT): gpu_density_accumulation.cu

gpu_fourwf.$(OBJEXT): gpu_fourwf.cu

gpu_mkkpg.$(OBJEXT): gpu_mkkpg.cu

gpu_nonlop.$(OBJEXT): gpu_nonlop.cu

gpu_sphere.$(OBJEXT): gpu_sphere.cu

prt_utils_rec.$(OBJEXT): prt_utils_rec.cu

rec_dens_calc.$(OBJEXT): rec_dens_calc.cu

rec_kernels.$(OBJEXT): rec_kernels.cu

recursion_bth.$(OBJEXT): recursion_bth.cu

recursion_cut_bth.$(OBJEXT): recursion_cut_bth.cu

recursion_cut_no_bth.$(OBJEXT): recursion_cut_no_bth.cu

recursion_no_bth.$(OBJEXT): recursion_no_bth.cu

EXTRA_DIST += abinit.amf

# FIXME: the following file belongs to nVidia Corp.
#        IT MUST BE REMOVED AS SOON AS POSSIBLE!!!
#nodist_lib52_manage_cuda_a_SOURCES = fortran.c

# FIXME: the following file belongs to nVidia Corp.
#        IT MUST BE REMOVED AS SOON AS POSSIBLE!!!
#EXTRA_DIST = deviceQuery.cu

#noinst_HEADERS = cuda_rec_head.h


EXTRA_DIST += _52_manage_cuda_
