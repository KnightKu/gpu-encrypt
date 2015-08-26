CUDA_PATH       := /usr/local/cuda-5.0
CUDA_INC_PATH   := $(CUDA_PATH)/include
CUDA_BIN_PATH   := $(CUDA_PATH)/bin
CUDA_LIB_PATH   := $(CUDA_PATH)/lib64
NVCC            := $(CUDA_BIN_PATH)/nvcc
GCC             := g++

LDFLAGS   := -L$(CUDA_LIB_PATH) -lcudart
CCFLAGS   := -m64
NVCCFLAGS := -m64

hill: main.cu
	$(NVCC) $(NVCCFLAGS) $(INCLUDS) $(LDFLAGS) -O3 -o $@ $< matInverse.cpp
