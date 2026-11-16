ifeq ($(OS),Windows_NT)
  MAKE_CMD := mingw32-make
  RUN_IN_SRC = cmd /c "pushd M*dulos*-* && pushd M*dulos* && $(MAKE_CMD) $(1) && popd && popd"
else
  MAKE_CMD := make
  RUN_IN_SRC = $(MAKE_CMD) -C "Módulos e main-20241119T180145Z-001/Módulos e main" $(1)
endif

.PHONY: all clean

all:
	$(call RUN_IN_SRC,all)

clean:
	$(call RUN_IN_SRC,clean)
