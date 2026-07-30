###############################################################################
# Atmospheric Chemistry module build fragment.
# Included by framework/app.mk automatically when atmospheric_chemistry is the
# current APPLICATION_NAME, sharing KPP and SUNDIALS config across the module
# build and any external app that enables ATMOSPHERIC_CHEMISTRY.
###############################################################################

ATMCHEM_APP_DIR := $(MOOSE_DIR)/modules/atmospheric_chemistry

# SUNDIALS direct ODE solver support for box-mode integration.
# When SUNDIALS is available, SundialsBoxIntegrator in BoxIntegrator.C uses
# SUNDIALS CVODE/ARKODE instead of PETSc TS for the solver_type=sundials path.
# No-op when SUNDIALS is not detected at build time.
SUNDIALS_DIR ?= $(CONDA_PREFIX)
ifneq ($(wildcard $(SUNDIALS_DIR)/include/cvode/cvode.h),)
  ADDITIONAL_CPPFLAGS += -I$(SUNDIALS_DIR)/include
  ADDITIONAL_LIBS      += -L$(SUNDIALS_DIR)/lib \
                          -Wl,-rpath,$(SUNDIALS_DIR)/lib \
                          -lsundials_cvode \
                          -lsundials_arkode \
                          -lsundials_nvecserial \
                          -lsundials_sunlinsoldense \
                          -lsundials_sunmatrixdense
  ADDITIONAL_CPPFLAGS  += -DHAVE_SUNDIALS
endif

# KPP shared-library backend support.
ADDITIONAL_CPPFLAGS += -I$(ATMCHEM_APP_DIR)/kpp/adapter -I$(ATMCHEM_APP_DIR)/kpp/runtime
ADDITIONAL_CPPFLAGS += -DKPP_ENABLED
