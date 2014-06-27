# Copyright InfoTeCS. All rights reserved.

include ../Makefile.cf.linux

PROJECT=signature_contest
TARGET_BINARY=$(BIN_DIR)/$(PROJECT)

OBJ= \
	src/linux/load_dl.o \
	src/main.o \


OBJ_IN_RESULT=$(addprefix $(OBJ_DIR)/$(PROJECT)/, $(OBJ))
TEMP_OBJ_DIRS_IN_RESULT:= $(sort $(dir $(OBJ_IN_RESULT) ))

CFLAGS +=
CXXFLAGS += -std=c++0x
LDFLAGS += -lstdc++ -ldl -s

all: $(TARGET_BINARY)

$(OBJ_DIR)/$(PROJECT)/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<
$(OBJ_DIR)/$(PROJECT)/%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

$(TARGET_BINARY): $(OBJ_IN_RESULT)
	$(CC) -o $(TARGET_BINARY) $(OBJ_IN_RESULT) $(LDFLAGS)

prepare: make_dirs_structure

clean:
	rm -rf $(RESULT_DIR)

make_dirs_structure:
	mkdir -p $(TEMP_OBJ_DIRS_IN_RESULT)
	mkdir -p $(BIN_DIR)

install:
	mkdir -p $(INCLUDE_DIR)/$(COMPONENT)
	cp -f ../sign_engine.h $(INCLUDE_DIR)/$(COMPONENT)
	cp -f $(TARGET_BINARY) $(INSTALL_BIN_DIR)/

.PHONY: all clean make_dirs_structure install

