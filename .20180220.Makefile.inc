# PYTHIA configuration file.
# Generated on Fri Oct 6 17:34:24 BST 2017 with the user supplied options:
# --with-root=/cvmfs/lhcb.cern.ch/lib/lcg/releases/LCG_84/ROOT/6.06.02/x86_64-slc6-gcc49-opt/
# --with-fastjet3=/home/nkw/gs2015/fastjet
# --with-evtgen-include=/home/nkw/gs2015/evtgen/1.5.0/
# --with-evtgen-lib=/home/nkw/gs2015/evtgen/1.5.0/lib
# --with-hepmc2=/home/nkw/gs2015/HepMC-2.06.08

# Install directory prefixes.
PREFIX_BIN=/home/nkw/pythia8230/bin
PREFIX_INCLUDE=/home/nkw/pythia8230/include
PREFIX_LIB=/home/nkw/pythia8230/lib
PREFIX_SHARE=/home/nkw/pythia8230/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
ENABLE_SHARED=false
CXX=g++
CXX_COMMON=-O2  -pedantic -W -Wall -Wshadow -fPIC
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so

# EVTGEN configuration.
EVTGEN_USE=true
EVTGEN_BIN=
EVTGEN_INCLUDE=/home/nkw/gs2015/evtgen/1.5.0/
EVTGEN_LIB=/home/nkw/gs2015/evtgen/1.5.0/lib

# FASTJET3 configuration.
FASTJET3_USE=true
FASTJET3_BIN=/home/nkw/gs2015/fastjet/bin/
FASTJET3_INCLUDE=/home/nkw/gs2015/fastjet/include
FASTJET3_LIB=/home/nkw/gs2015/fastjet/lib

# HEPMC2 configuration.
HEPMC2_USE=true
HEPMC2_BIN=/home/nkw/gs2015/HepMC-2.06.08/bin/
HEPMC2_INCLUDE=/home/nkw/gs2015/HepMC-2.06.08/include
HEPMC2_LIB=/home/nkw/gs2015/HepMC-2.06.08/lib

# HEPMC3 configuration.
HEPMC3_USE=false
HEPMC3_BIN=
HEPMC3_INCLUDE=./
HEPMC3_LIB=./

# LHAPDF5 configuration.
LHAPDF5_USE=false
LHAPDF5_BIN=
LHAPDF5_INCLUDE=./
LHAPDF5_LIB=./

# LHAPDF6 configuration.
LHAPDF6_USE=false
LHAPDF6_BIN=
LHAPDF6_INCLUDE=./
LHAPDF6_LIB=./

# POWHEG configuration.
POWHEG_USE=false
POWHEG_BIN=
POWHEG_INCLUDE=./
POWHEG_LIB=./

# PROMC configuration.
PROMC_USE=false
PROMC_BIN=
PROMC_INCLUDE=./
PROMC_LIB=./

# ROOT configuration.
ROOT_USE=true
ROOT_BIN=/cvmfs/lhcb.cern.ch/lib/lcg/releases/LCG_84/ROOT/6.06.02/x86_64-slc6-gcc49-opt//bin/
ROOT_INCLUDE=/cvmfs/lhcb.cern.ch/lib/lcg/releases/LCG_84/ROOT/6.06.02/x86_64-slc6-gcc49-opt//include
ROOT_LIB=/cvmfs/lhcb.cern.ch/lib/lcg/releases/LCG_84/ROOT/6.06.02/x86_64-slc6-gcc49-opt//lib

# GZIP configuration.
GZIP_USE=false
GZIP_BIN=
GZIP_INCLUDE=./
GZIP_LIB=./

# BOOST configuration.
BOOST_USE=false
BOOST_BIN=
BOOST_INCLUDE=./
BOOST_LIB=./

# PYTHON configuration.
PYTHON_USE=false
PYTHON_BIN=
PYTHON_INCLUDE=./
PYTHON_LIB=./
