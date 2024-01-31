# Compile HGSpy
#   Compile with g++ or Intel C++ Compiler (/opt/intel/bin/icpc)
#   Compile with the most aggressive optimization setting (O3)
#   Use the most pedantic compiler settings: must compile with no warnings at all
#
# The user may override any desired internal variable by redefining it via command-line:
#   make CXX=g++ [...]
#   make OPTL=-O2 [...]
#   make FLAGS="-Wall -g" [...]
#
# Caleb Fuster, Arnau Miro, Manel Soria
# UPC - ESEIAAT

# Compilers
#
# Automatically detect if the intel compilers are installed and use
# them, otherwise default to the GNU compilers
PYTHON = python
PIP    = pip

# One rule to compile them all, one rule to find them,
# One rule to bring them all and in the compiler link them.
all:  requirements install
	@echo ""
	@echo "HGSpy deployed successfully"


# Python
#
requirements: requirements.txt
	@${PIP} install -r $<

install: 
	@${PIP} install .

install_dev: 
	@${PIP} install -e .

package-build:
	@${PYTHON} -m build


# Clean
#
clean:
	-@cd HGSpy; rm -rf __pycache__
uninstall: clean
	-@${PIP} uninstall HGSpy
	-@rm -rf HGSpy.egg-info