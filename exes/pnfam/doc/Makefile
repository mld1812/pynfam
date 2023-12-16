docGen=doxygen
TARGET=Doxyfile
outputDir=html latex
DOXYHOME=$(shell pwd)

all:
	@echo "Compiling doxygen documentation..."
	@$(docGen) $(TARGET)

clean:
	@rm -rf $(outputDir)
