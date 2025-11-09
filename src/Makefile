# ===============================
# Makefile for CS6135 HW3
# ===============================

CXX := g++
CXXFLAGS := -O3 -march=native -DNDEBUG -flto -funroll-loops -fomit-frame-pointer -ffast-math -std=c++17

BIN_DIR := ../bin
OUT_DIR := ../output
SRC := hw3.cpp
TARGET := $(BIN_DIR)/hw3

.PHONY: all clean run

# -------------------------------
# Default build rule
# -------------------------------
all: $(TARGET)

$(TARGET): $(SRC)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@
	@echo "Build complete: $(TARGET)"

clean:
	@rm -f $(TARGET)
	@echo "Cleaned executable"

# -------------------------------
# Run all testcases
# 使用方式：
#   make run
# -------------------------------
run1: $(TARGET)
	@$(TARGET) ../testcase/public1.lef ../testcase/public1.def $(OUT_DIR)/public1_out.def

run: $(TARGET)
	@echo "▶ Running all DEF/LEF testcases..."
	@$(TARGET) ../testcase/public1.lef ../testcase/public1.def $(OUT_DIR)/public1.def
	@$(TARGET) ../testcase/public2.lef ../testcase/public2.def $(OUT_DIR)/public2.def
	@$(TARGET) ../testcase/public3.lef ../testcase/public3.def $(OUT_DIR)/public3.def
	@$(TARGET) ../testcase/public4.lef ../testcase/public4.def $(OUT_DIR)/public4.def
	
	@echo "Logs saved to ../log_public*.txt"
