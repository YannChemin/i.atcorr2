MODULE_TOPDIR = $(HOME)/dev/grass

PGM = i.atcorr2

# ── grass_sixsv dependency ──────────────────────────────────────────────────
# grass_sixsv is the standalone 6SV2.1 library built and installed from
# ../libsixsv/ (https://github.com/YannChemin/libsixsv).
# The lazy = assignment ensures GRASS_LIB_VERSION_NUMBER is resolved after
# the include chain (Vars.make → Grass.make) has run.
SIXSV_LIB_NAME = grass_sixsv.$(GRASS_LIB_VERSION_NUMBER)
SIXSVLIB        = -l$(SIXSV_LIB_NAME)
SIXSVDEP        = $(ARCH_LIBDIR)/$(SHLIB_PREFIX)$(SIXSV_LIB_NAME)$(SHLIB_SUFFIX)

# ── Sources ─────────────────────────────────────────────────────────────────
MOD_OBJS = main.o

# ── Compiler / linker options ────────────────────────────────────────────────
# atcorr.h lives in ../libsixsv/include/ in the build tree; after
# install it is under $GISBASE/include/grass/ via the standard INC path.
EXTRA_INC     = -I$(MODULE_TOPDIR)/../libsixsv/include
EXTRA_CFLAGS  = -O3 -ffast-math -fopenmp -std=c11 \
                -Wall -Wextra -Wno-unused-parameter
EXTRA_LDFLAGS = -fopenmp
LIBES         = $(SIXSVLIB) $(RASTERLIB) $(GISLIB) $(MATHLIB)
DEPENDENCIES  = $(SIXSVDEP)

include $(MODULE_TOPDIR)/include/Make/Module.make

default: cmd

# ── Installation ─────────────────────────────────────────────────────────────
install:
	$(INSTALL) $(ARCH_DISTDIR)/bin/$(PGM)$(EXE) $(INST_DIR)/bin/
	$(INSTALL_DATA) i.atcorr2.html $(INST_DIR)/docs/html/
	$(INSTALL_DATA) $(ARCH_DISTDIR)/docs/man/man1/$(PGM).1 \
	    $(INST_DIR)/docs/man/man1/ 2>/dev/null || true

# ── Python module installation ────────────────────────────────────────────────
# Installs the Python front-end, the sensors helper, and the SRF CSV data.
# sensors.py locates sensors_csv/ relative to __file__, so both must land
# in $(INST_DIR)/scripts/.
install_py:
	$(INSTALL) i.atcorr2.py $(INST_DIR)/scripts/$(PGM)
	$(INSTALL_DATA) sensors.py $(INST_DIR)/scripts/sensors.py
	$(INSTALL) -d $(INST_DIR)/scripts/sensors_csv
	$(INSTALL_DATA) sensors_csv/*.csv $(INST_DIR)/scripts/sensors_csv/
	$(INSTALL_DATA) sensors_csv/README $(INST_DIR)/scripts/sensors_csv/
	$(INSTALL_DATA) i.atcorr2.html $(INST_DIR)/docs/html/

.PHONY: install install_py
