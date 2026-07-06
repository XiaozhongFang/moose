###############################################################################
################### MOOSE Application Standard Makefile #######################
###############################################################################
#
# Optional Environment variables
# MOOSE_DIR        - Root directory of the MOOSE project
# MODULE_DIR       - Location of the MOOSE modules directory
# FRAMEWORK_DIR    - Location of the MOOSE framework
#
###############################################################################
MODULE_DIR         ?= $(shell dirname `pwd`)
MOOSE_DIR          ?= $(shell dirname $(MODULE_DIR))
FRAMEWORK_DIR      ?= $(MOOSE_DIR)/framework
ADDITIONAL_CPPFLAGS += -Wall -Wextra
###############################################################################

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

# framework
include $(FRAMEWORK_DIR)/build.mk
include $(FRAMEWORK_DIR)/moose.mk

# Module
MODULE_NAME := atmospheric_chemistry
include $(MODULE_DIR)/modules.mk

###############################################################################
# Additional special case targets should be added here
