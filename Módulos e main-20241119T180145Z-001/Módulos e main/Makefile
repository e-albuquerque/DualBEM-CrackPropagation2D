FC      := gfortran
FFLAGS  := -O2 -Wall -Wextra
STATIC  ?= 1
LIBS    := -llapack -lblas

ifeq ($(STATIC),1)
  LDFLAGS := -static
else
  LDFLAGS :=
endif

BUILD_DIR := build
BIN_DIR   := ../../bin
TARGET    := $(BIN_DIR)/fracture_2d

ifeq ($(OS),Windows_NT)
  TARGET := $(BIN_DIR)/fracture_2d.exe
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

MOD_DIR := modules

OBJS := \
  $(BUILD_DIR)/error_messages.o \
  $(BUILD_DIR)/geometry.o \
  $(BUILD_DIR)/kelvin_fund_sol_2d.o \
  $(BUILD_DIR)/numerical_integration_points.o \
  $(BUILD_DIR)/vector_algebra.o \
  $(BUILD_DIR)/matrix_algebra.o \
  $(BUILD_DIR)/interpolation_functions.o \
  $(BUILD_DIR)/fracture_2d_model.o \
  $(BUILD_DIR)/fracture_2d_create_internal_points.o \
  $(BUILD_DIR)/fracture_2d_displ_correlation.o \
  $(BUILD_DIR)/fracture_2d_regular_integration.o \
  $(BUILD_DIR)/fracture_2d_weakly_sing_integration.o \
  $(BUILD_DIR)/fracture_2d_rad_exp_intgr.o \
  $(BUILD_DIR)/fracture_2d_assembler.o \
  $(BUILD_DIR)/fracture_2d_set_bound_results.o \
  $(BUILD_DIR)/fracture_2d_internal_points.o \
  $(BUILD_DIR)/fracture_2d_j_integral.o \
  $(BUILD_DIR)/fracture_2d_stress_intens_factor.o \
  $(BUILD_DIR)/fracture_2d_crack_propagation.o \
  $(BUILD_DIR)/fracture_2d_output.o \
  $(BUILD_DIR)/fracture_2d_rainflow.o \
  $(BUILD_DIR)/fracture_2d_input.o \
  $(BUILD_DIR)/Fracture_2D_v3-0.o

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(MKDIR_BIN)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

$(BUILD_DIR)/%.o: $(MOD_DIR)/%.f95
	$(MKDIR_BUILD)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/Fracture_2D_v3-0.o: Fracture_2D_v3-0.f95 \
  $(BUILD_DIR)/fracture_2d_input.o \
  $(BUILD_DIR)/fracture_2d_assembler.o \
  $(BUILD_DIR)/fracture_2d_set_bound_results.o \
  $(BUILD_DIR)/fracture_2d_internal_points.o \
  $(BUILD_DIR)/fracture_2d_output.o \
  $(BUILD_DIR)/fracture_2d_model.o \
  $(BUILD_DIR)/fracture_2d_stress_intens_factor.o \
  $(BUILD_DIR)/fracture_2d_crack_propagation.o \
  $(BUILD_DIR)/fracture_2d_create_internal_points.o \
  $(BUILD_DIR)/matrix_algebra.o
	$(MKDIR_BUILD)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR)/vector_algebra.o: $(BUILD_DIR)/error_messages.o
$(BUILD_DIR)/matrix_algebra.o: $(BUILD_DIR)/error_messages.o
$(BUILD_DIR)/interpolation_functions.o: $(BUILD_DIR)/vector_algebra.o $(BUILD_DIR)/error_messages.o $(BUILD_DIR)/geometry.o
$(BUILD_DIR)/fracture_2d_model.o: $(BUILD_DIR)/geometry.o $(BUILD_DIR)/vector_algebra.o $(BUILD_DIR)/interpolation_functions.o $(BUILD_DIR)/error_messages.o
$(BUILD_DIR)/fracture_2d_create_internal_points.o: $(BUILD_DIR)/fracture_2d_model.o $(BUILD_DIR)/numerical_integration_points.o
$(BUILD_DIR)/fracture_2d_displ_correlation.o: $(BUILD_DIR)/fracture_2d_model.o
$(BUILD_DIR)/fracture_2d_regular_integration.o: $(BUILD_DIR)/fracture_2d_model.o $(BUILD_DIR)/error_messages.o $(BUILD_DIR)/geometry.o $(BUILD_DIR)/kelvin_fund_sol_2d.o $(BUILD_DIR)/interpolation_functions.o $(BUILD_DIR)/numerical_integration_points.o
$(BUILD_DIR)/fracture_2d_weakly_sing_integration.o: $(BUILD_DIR)/interpolation_functions.o $(BUILD_DIR)/numerical_integration_points.o $(BUILD_DIR)/fracture_2d_model.o
$(BUILD_DIR)/fracture_2d_rad_exp_intgr.o: $(BUILD_DIR)/error_messages.o $(BUILD_DIR)/geometry.o $(BUILD_DIR)/matrix_algebra.o $(BUILD_DIR)/fracture_2d_model.o $(BUILD_DIR)/numerical_integration_points.o $(BUILD_DIR)/kelvin_fund_sol_2d.o
$(BUILD_DIR)/fracture_2d_assembler.o: $(BUILD_DIR)/fracture_2d_model.o $(BUILD_DIR)/fracture_2d_regular_integration.o $(BUILD_DIR)/fracture_2d_weakly_sing_integration.o $(BUILD_DIR)/fracture_2d_rad_exp_intgr.o
$(BUILD_DIR)/fracture_2d_set_bound_results.o: $(BUILD_DIR)/error_messages.o $(BUILD_DIR)/interpolation_functions.o $(BUILD_DIR)/fracture_2d_model.o $(BUILD_DIR)/fracture_2d_rad_exp_intgr.o $(BUILD_DIR)/fracture_2d_regular_integration.o
$(BUILD_DIR)/fracture_2d_internal_points.o: $(BUILD_DIR)/fracture_2d_model.o $(BUILD_DIR)/fracture_2d_regular_integration.o $(BUILD_DIR)/kelvin_fund_sol_2d.o $(BUILD_DIR)/numerical_integration_points.o
$(BUILD_DIR)/fracture_2d_j_integral.o: $(BUILD_DIR)/fracture_2d_model.o $(BUILD_DIR)/fracture_2d_regular_integration.o $(BUILD_DIR)/fracture_2d_rad_exp_intgr.o $(BUILD_DIR)/kelvin_fund_sol_2d.o $(BUILD_DIR)/numerical_integration_points.o
$(BUILD_DIR)/fracture_2d_stress_intens_factor.o: $(BUILD_DIR)/fracture_2d_displ_correlation.o $(BUILD_DIR)/fracture_2d_j_integral.o
$(BUILD_DIR)/fracture_2d_crack_propagation.o: $(BUILD_DIR)/fracture_2d_model.o $(BUILD_DIR)/matrix_algebra.o $(BUILD_DIR)/fracture_2d_j_integral.o $(BUILD_DIR)/fracture_2d_assembler.o $(BUILD_DIR)/fracture_2d_set_bound_results.o $(BUILD_DIR)/fracture_2d_stress_intens_factor.o
$(BUILD_DIR)/fracture_2d_output.o: $(BUILD_DIR)/error_messages.o $(BUILD_DIR)/interpolation_functions.o $(BUILD_DIR)/fracture_2d_model.o $(BUILD_DIR)/fracture_2d_stress_intens_factor.o $(BUILD_DIR)/fracture_2d_crack_propagation.o
$(BUILD_DIR)/fracture_2d_rainflow.o: $(BUILD_DIR)/fracture_2d_model.o $(BUILD_DIR)/fracture_2d_crack_propagation.o $(BUILD_DIR)/fracture_2d_output.o
$(BUILD_DIR)/fracture_2d_input.o: $(BUILD_DIR)/error_messages.o $(BUILD_DIR)/geometry.o $(BUILD_DIR)/interpolation_functions.o $(BUILD_DIR)/fracture_2d_model.o $(BUILD_DIR)/fracture_2d_rad_exp_intgr.o $(BUILD_DIR)/fracture_2d_stress_intens_factor.o $(BUILD_DIR)/fracture_2d_crack_propagation.o $(BUILD_DIR)/fracture_2d_rainflow.o

clean:
	$(RM_RF)
	$(RM_BIN)
