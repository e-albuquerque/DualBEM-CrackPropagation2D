# ==============================================================================
# Makefile - Fracture 2D (Dual-BEM com Elementos Quadráticos Descontínuos)
# Compilador: GNU Fortran (gfortran) + Netlib LAPACK/BLAS + OpenMP
# ==============================================================================

FC      := gfortran
OPENMP  ?= 1

ifeq ($(OS),Windows_NT)
  STATIC ?= 1
else
  STATIC ?= 0
endif

FFLAGS  := -O3 -Wall -Wextra -Wno-unused-dummy-argument -Wno-maybe-uninitialized
LDFLAGS :=
LIBS    := -llapack -lblas

ifeq ($(OPENMP),1)
  FFLAGS  += -fopenmp
  LDFLAGS += -fopenmp
endif

ifeq ($(STATIC),1)
  LDFLAGS += -static
endif

SRC_DIR   := src
BUILD_DIR := build
BIN_DIR   := bin
TARGET    := $(BIN_DIR)/fracture_2d

ifeq ($(OS),Windows_NT)
  TARGET      := $(BIN_DIR)/fracture_2d.exe
  MKDIR_BUILD = if not exist "$(subst /,\,$(BUILD_DIR))" mkdir "$(subst /,\,$(BUILD_DIR))"
  MKDIR_BIN   = if not exist "$(subst /,\,$(BIN_DIR))" mkdir "$(subst /,\,$(BIN_DIR))"
  RM_RF       = if exist "$(subst /,\,$(BUILD_DIR))" rmdir /s /q "$(subst /,\,$(BUILD_DIR))"
  RM_BIN      = if exist "$(subst /,\,$(TARGET))" del /q "$(subst /,\,$(TARGET))"
else
  MKDIR_BUILD = mkdir -p $(BUILD_DIR)
  MKDIR_BIN   = mkdir -p $(BIN_DIR)
  RM_RF       = rm -rf $(BUILD_DIR)
  RM_BIN      = rm -f $(TARGET)
endif

OBJS := \
	$(BUILD_DIR)/constants_mod.o \
	$(BUILD_DIR)/quadrature_mod.o \
	$(BUILD_DIR)/model_mod.o \
	$(BUILD_DIR)/kelvin_sol_mod.o \
	$(BUILD_DIR)/integration_mod.o \
	$(BUILD_DIR)/assembler_mod.o \
	$(BUILD_DIR)/sif_mod.o \
	$(BUILD_DIR)/fatigue_prop_mod.o \
	$(BUILD_DIR)/io_mod.o \
	$(BUILD_DIR)/main.o

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	@$(MKDIR_BIN)
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90
	@$(MKDIR_BUILD)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

# Dependências explícitas entre módulos
$(BUILD_DIR)/quadrature_mod.o:   $(BUILD_DIR)/constants_mod.o
$(BUILD_DIR)/model_mod.o:        $(BUILD_DIR)/constants_mod.o $(BUILD_DIR)/quadrature_mod.o
$(BUILD_DIR)/kelvin_sol_mod.o:   $(BUILD_DIR)/constants_mod.o
$(BUILD_DIR)/integration_mod.o:  $(BUILD_DIR)/constants_mod.o $(BUILD_DIR)/quadrature_mod.o $(BUILD_DIR)/model_mod.o $(BUILD_DIR)/kelvin_sol_mod.o
$(BUILD_DIR)/assembler_mod.o:    $(BUILD_DIR)/constants_mod.o $(BUILD_DIR)/quadrature_mod.o $(BUILD_DIR)/model_mod.o $(BUILD_DIR)/kelvin_sol_mod.o $(BUILD_DIR)/integration_mod.o
$(BUILD_DIR)/sif_mod.o:          $(BUILD_DIR)/constants_mod.o $(BUILD_DIR)/quadrature_mod.o $(BUILD_DIR)/model_mod.o $(BUILD_DIR)/kelvin_sol_mod.o $(BUILD_DIR)/integration_mod.o
$(BUILD_DIR)/fatigue_prop_mod.o: $(BUILD_DIR)/constants_mod.o $(BUILD_DIR)/quadrature_mod.o $(BUILD_DIR)/model_mod.o $(BUILD_DIR)/assembler_mod.o $(BUILD_DIR)/sif_mod.o
$(BUILD_DIR)/io_mod.o:           $(BUILD_DIR)/constants_mod.o $(BUILD_DIR)/quadrature_mod.o $(BUILD_DIR)/model_mod.o $(BUILD_DIR)/integration_mod.o $(BUILD_DIR)/sif_mod.o $(BUILD_DIR)/fatigue_prop_mod.o
$(BUILD_DIR)/main.o:             $(BUILD_DIR)/constants_mod.o $(BUILD_DIR)/model_mod.o $(BUILD_DIR)/assembler_mod.o $(BUILD_DIR)/sif_mod.o $(BUILD_DIR)/fatigue_prop_mod.o $(BUILD_DIR)/io_mod.o

clean:
	$(RM_RF)
	$(RM_BIN)
