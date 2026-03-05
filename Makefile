# This file is a part of DUDI, the Fortran-95 implementation
# of the two-body model for dust dynamics
# Version 1.2.2
# License: GPL-3.0

# =========================
#  DUDI explicit targets
# =========================
# Usage:
#   make flyby_profile
#   make run-flyby       # run flyby_profile with params 5, 17, 7.1, 7.2, 7.3, 21.1, 21.2, 21.3
#   make list | make clean | make distclean

# -------- Compiler settings --------
FC      ?= gfortran
FFLAGS  ?= -O3 -fimplicit-none -Wall -Wno-tabs -Wno-unused-variable
# strict flags
STRICT_WARNINGS = -O0 -g -fimplicit-none -Wall -Wextra -Wconversion -Wsurprising \
                  -Warray-temporaries -Wcharacter-truncation -Wreal-q-constant \
                  -Wtarget-lifetime -Wimplicit-interface
LDFLAGS ?=
PYTHON  ?= python3

# Refuse/override Fortran 77 compilers
NEEDS_FC_GUARD := $(filter-out clean distclean help list,$(MAKECMDGOALS))
ifneq ($(NEEDS_FC_GUARD),)
  ifneq (,$(findstring f77,$(FC)))
    $(info FC='$(FC)' looks like Fortran 77; switching to gfortran)
    override FC := gfortran
  endif
endif

# Avoid built-in implicit rules that might pick f77
.SUFFIXES:
MAKEFLAGS += --no-builtin-rules

# -------- Layout --------
MODDIR := build
BINDIR := bin
SRCDIR := src
EXDIR  := examples
RESDIR := results

# -------- Core sources (ordered: providers before users) --------
CORE_SOURCES := \
  $(SRCDIR)/const.f90 \
  $(SRCDIR)/comparison_utils.f90 \
  $(SRCDIR)/define_types.f90 \
  $(SRCDIR)/help.f90 \
  $(SRCDIR)/bgmod.f90 \
  $(SRCDIR)/distributions_fun.f90 \
  $(SRCDIR)/gu.f90 \
  $(SRCDIR)/twobody_fun.f90 \
  $(SRCDIR)/inputdata.f90 \
  $(SRCDIR)/dataoutmod.f90 \
  $(SRCDIR)/integrator.f90 \
  $(SRCDIR)/image_construction.f90

# -------- Program sources --------
FLYBY_PROFILE_SRC        ?= $(EXDIR)/flyby_profile.f90
VERTICAL_STRUCTURE_SRC   ?= $(EXDIR)/plume_vert_slice.f90
HORIZONTAL_STRUCTURE_SRC ?= $(EXDIR)/plume_horizontal_structure.f90

# -------- Phony targets --------
.PHONY: all help list clean distclean flyby_profile run-flyby vertical_structure horizontal_structure

# Default: build flyby_profile
all: flyby_profile

help:
	@echo "Build:"
	@echo "  make flyby_profile"
	@echo "  make vertical_structure"
	@echo "  make horizontal_structure"
	@echo ""
	@echo "Run:"
	@echo "  make run-flyby   # run flyby_profile with params 5, 5.2, 17, 17.17, 7.1, 7.2, 7.3, 21.1, 21.2, 21.3"
	@echo ""
	@echo "Other:"
	@echo "  make list | make clean | make distclean"

list:
	@echo "Core src:               $(CORE_SOURCES)"
	@echo "FLYBY_PROFILE_SRC:      $(FLYBY_PROFILE_SRC)"
	@echo "VERTICAL_STRUCTURE_SRC: $(VERTICAL_STRUCTURE_SRC)"
	@echo "HORIZONTAL_STRUCTURE_SRC: $(HORIZONTAL_STRUCTURE_SRC)"

# Ensure dirs exist
$(MODDIR) $(BINDIR) $(RESDIR):
	mkdir -p $@

# -----------------------------
#  Build rules (explicit only)
# -----------------------------

$(BINDIR)/flyby_profile: $(CORE_SOURCES) $(FLYBY_PROFILE_SRC) | $(BINDIR) $(MODDIR)
	$(FC) $(FFLAGS) -fopenmp -J$(MODDIR) -I$(MODDIR) $(CORE_SOURCES) $(FLYBY_PROFILE_SRC) -o $@ $(LDFLAGS) -fopenmp
flyby_profile: $(BINDIR)/flyby_profile

$(BINDIR)/vertical_structure: $(CORE_SOURCES) $(VERTICAL_STRUCTURE_SRC) | $(BINDIR) $(MODDIR)
	$(FC) $(FFLAGS) -fopenmp -J$(MODDIR) -I$(MODDIR) $(CORE_SOURCES) $(VERTICAL_STRUCTURE_SRC) -o $@ $(LDFLAGS) -fopenmp
vertical_structure: $(BINDIR)/vertical_structure

$(BINDIR)/horizontal_structure: $(CORE_SOURCES) $(HORIZONTAL_STRUCTURE_SRC) | $(BINDIR) $(MODDIR)
	$(FC) $(FFLAGS) -fopenmp -J$(MODDIR) -I$(MODDIR) $(CORE_SOURCES) $(HORIZONTAL_STRUCTURE_SRC) -o $@ $(LDFLAGS) -fopenmp
horizontal_structure: $(BINDIR)/horizontal_structure

run-flyby: $(BINDIR)/flyby_profile | $(RESDIR)
	$(BINDIR)/flyby_profile 5
	$(BINDIR)/flyby_profile 5.2
	$(BINDIR)/flyby_profile 17
	$(BINDIR)/flyby_profile 17.17
	$(BINDIR)/flyby_profile 7.1
	$(BINDIR)/flyby_profile 7.2
	$(BINDIR)/flyby_profile 7.3
	$(BINDIR)/flyby_profile 21.1
	$(BINDIR)/flyby_profile 21.2
	$(BINDIR)/flyby_profile 21.3

# -----------------------------
#  Cleaning
# -----------------------------
clean:
	@rm -f $(MODDIR)/*.o $(MODDIR)/*.mod

distclean: clean
	@rm -rf $(BINDIR) $(RESDIR)/*

# Strict warnings sweep: clean + rebuild (compile only)
clean-warnings:
	$(MAKE) clean
	$(MAKE) -B flyby_profile FFLAGS='$(STRICT_WARNINGS)'
