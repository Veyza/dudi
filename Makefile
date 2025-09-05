# This file is a part of DUDI, the Fortran-95 implementation 
# of the two-body model for dust dynamics
# Version 1.2.1
# License: GPL-3.0

# =========================
#  DUDI explicit targets
# =========================
# Usage:
#   make [dudi]
#   make enceladus | europa | io        # build *_dudi, run it, then plot
#   make enceladus_dudi | europa_dudi | io_dudi   # build-only
#   make run-dudi | run-enceladus_dudi | run-europa_dudi | run-io_dudi
#   make list | make clean | make distclean

# -------- Compiler settings --------
FC      ?= gfortran
#FFLAGS  ?= -O3 -fimplicit-none -Wall -Wno-tabs -Wno-unused-variable
FFLAGS ?= -O3 -fimplicit-none -Wno-tabs -Wno-unused-variable
LDFLAGS ?=
PYTHON  ?= python3

# -------- Layout --------
MODDIR := build
BINDIR := bin
SRCDIR := src
EXDIR  := examples
RESDIR := results

# -------- Core sources --------
CORE_SOURCES := $(wildcard $(SRCDIR)/*.f90)

# -------- Example/main program sources --------
MAIN_SRC            ?= $(EXDIR)/main_program.f90
ENCELADUS_SRC       ?= $(EXDIR)/enceladus_example.f90
EUROPA_SRC          ?= $(EXDIR)/europa_example.f90
IO_SRC              ?= $(EXDIR)/io_example.f90

# ---- Plot scripts ----
ENCELADUS_PYSCRIPT ?= scripts/e2plot.py
EUROPA_PYSCRIPT    ?= scripts/deposition.py
IO_PYSCRIPT        ?= scripts/volcano_image.py

# -------- Phony targets --------
.PHONY: all help list clean distclean \
        dudi enceladus europa io enceladus_dudi europa_dudi io_dudi \
        run-dudi run-enceladus_dudi run-europa_dudi run-io_dudi

# Default: build the library demo binary
all: dudi

help:
	@echo "Build-only binaries:"
	@echo "  make dudi | enceladus_dudi | europa_dudi | io_dudi"
	@echo ""
	@echo "Pipelines (build -> run -> plot):"
	@echo "  make enceladus | europa | io"
	@echo ""
	@echo "Run-only helpers:"
	@echo "  make run-dudi | run-enceladus_dudi | run-europa_dudi | run-io_dudi"
	@echo ""
	@echo "Override plot scripts/Python if needed:"
	@echo "  ENCELADUS_PYSCRIPT=$(ENCELADUS_PYSCRIPT)"
	@echo "  EUROPA_PYSCRIPT=$(EUROPA_PYSCRIPT)"
	@echo "  IO_PYSCRIPT=$(IO_PYSCRIPT)"
	@echo "  PYTHON=$(PYTHON)"

list:
	@echo "Core src:          $(CORE_SOURCES)"
	@echo "MAIN_SRC:          $(MAIN_SRC)"
	@echo "ENCELADUS_SRC:     $(ENCELADUS_SRC)"
	@echo "EUROPA_SRC:        $(EUROPA_SRC)"
	@echo "IO_SRC:            $(IO_SRC)"
	@echo "Plot scripts: enc=$(ENCELADUS_PYSCRIPT)  eur=$(EUROPA_PYSCRIPT)  io=$(IO_PYSCRIPT)"

# Ensure dirs exist
$(MODDIR) $(BINDIR) $(RESDIR):
	mkdir -p $@

# -----------------------------
#  Build rules (explicit only)
# -----------------------------

# dudi (library demo)
$(BINDIR)/dudi: $(CORE_SOURCES) $(MAIN_SRC) | $(BINDIR) $(MODDIR)
	$(FC) $(FFLAGS) -J$(MODDIR) -I$(MODDIR) $(CORE_SOURCES) $(MAIN_SRC) -o $@ $(LDFLAGS)
dudi: $(BINDIR)/dudi

# *_dudi binaries (build-only)
$(BINDIR)/enceladus_dudi: $(CORE_SOURCES) $(ENCELADUS_SRC) | $(BINDIR) $(MODDIR)
	$(FC) $(FFLAGS) -J$(MODDIR) -I$(MODDIR) $(CORE_SOURCES) $(ENCELADUS_SRC) -o $@ $(LDFLAGS)
enceladus_dudi: $(BINDIR)/enceladus_dudi

$(BINDIR)/europa_dudi: $(CORE_SOURCES) $(EUROPA_SRC) | $(BINDIR) $(MODDIR)
	$(FC) $(FFLAGS) -J$(MODDIR) -I$(MODDIR) $(CORE_SOURCES) $(EUROPA_SRC) -o $@ $(LDFLAGS)
europa_dudi: $(BINDIR)/europa_dudi

$(BINDIR)/io_dudi: $(CORE_SOURCES) $(IO_SRC) | $(BINDIR) $(MODDIR)
	$(FC) $(FFLAGS) -J$(MODDIR) -I$(MODDIR) $(CORE_SOURCES) $(IO_SRC) -o $@ $(LDFLAGS)
io_dudi: $(BINDIR)/io_dudi

# -----------------------------
#  Pipelines: build -> run -> plot
# -----------------------------
enceladus: $(BINDIR)/enceladus_dudi | $(RESDIR)
	@echo ">>> Running enceladus_dudi"
	$(BINDIR)/enceladus_dudi
	@echo ">>> Plotting with: $(PYTHON) $(ENCELADUS_PYSCRIPT)"
	$(PYTHON) $(ENCELADUS_PYSCRIPT)

europa: $(BINDIR)/europa_dudi | $(RESDIR)
	@echo ">>> Running europa_dudi"
	$(BINDIR)/europa_dudi
	@echo ">>> Plotting with: $(PYTHON) $(EUROPA_PYSCRIPT)"
	$(PYTHON) $(EUROPA_PYSCRIPT)

io: $(BINDIR)/io_dudi | $(RESDIR)
	@echo ">>> Running io_dudi"
	$(BINDIR)/io_dudi
	@echo ">>> Plotting with: $(PYTHON) $(IO_PYSCRIPT)"
	$(PYTHON) $(IO_PYSCRIPT)

# -----------------------------
#  Run helpers (binary only)
# -----------------------------
run-dudi:           $(BINDIR)/dudi           ; $(BINDIR)/dudi
run-enceladus_dudi: $(BINDIR)/enceladus_dudi ; $(BINDIR)/enceladus_dudi
run-europa_dudi:    $(BINDIR)/europa_dudi    ; $(BINDIR)/europa_dudi
run-io_dudi:        $(BINDIR)/io_dudi        ; $(BINDIR)/io_dudi

# -----------------------------
#  Cleaning
# -----------------------------
clean:
	@rm -f $(MODDIR)/*.o $(MODDIR)/*.mod

distclean: clean
	@rm -rf $(BINDIR) $(RESDIR)/*

