# Minimal Makefile for Sphinx documentation

# Configurable variables, can be overridden from the command line
# SPHINXOPTS: Options passed to the sphinx-build command
# Default value: '-W --keep-going' to treat warnings as errors and continue on errors
SPHINXOPTS    ?= -W --keep-going

# SPHINXBUILD: Command to run Sphinx's build tool
# Default value: Automatically detects the 'sphinx-build' executable
SPHINXBUILD   ?= $(shell which sphinx-build)

# SOURCEDIR: Directory containing the source (.rst) files
SOURCEDIR     = .

# BUILDDIR: Directory where generated documentation will be stored
# This is customizable, default is '_build'
BUILDDIR      ?= _build

# The 'help' target provides guidance on available make commands
# It is executed when running 'make' without additional arguments
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# The 'clean' target removes all build artifacts from the build directory
# Use: 'make clean' to clean the output directory
clean:
	rm -rf $(BUILDDIR)/*
	@echo "Build directory cleaned."

# This is a generic "catch-all" target that routes all unknown
# make commands to Sphinx. It uses the '-M' "make mode" option.
# $(O) is passed as a shortcut for additional Sphinx options.
.PHONY: help clean Makefile

%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)