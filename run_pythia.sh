#!/bin/bash
#PBS -j oe
#PBS -l walltime=120:00:00

#set -x

shopt -s expand_aliases

# CHange nkw to your own userid, e.g. user27, on the next line.
export USER=user72

cd /home/$USER


# Change this to your own directory, e.g.
#export MYPYTHIA=/home/user74/year4/Project/pythia_example_newest

export PATH=/bin:/usr/local/bin:/usr/bin:$PATH
export HOME=/home/$USER

export MYPYTHIA=$HOME/Project/Dark_Photons_Project


cd $MYPYTHIA
# Try to set this up explicitly here, alias does not seem to have been set up for batch workers.
source /cvmfs/lhcb.cern.ch/group_login.sh
#source ./setup_env.sh
SetupProject root


# Create a separate output file for each job submitted
#export workdir=pythia_$$

#mkdir -p $workdir
#cd $workdir

# Run the executable, assume compiled already before submitting to the batch system

#../gs_ex300 > run.out
cd $MYPYTHIA
./gs_project36 > run.out




