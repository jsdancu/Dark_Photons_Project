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
#SetupProject root
lb-run ROOT/6.12.04 ./gs_project26 > run26.out$$
lb-run ROOT/6.12.04 ./gs_project27 > run27.out$$
lb-run ROOT/6.12.04 ./gs_project28 > run28.out$$
lb-run ROOT/6.12.04 ./gs_project29 > run29.out$$
lb-run ROOT/6.12.04 ./gs_project210 > run210.out$$
lb-run ROOT/6.12.04 ./gs_project211 > run211.out$$
lb-run ROOT/6.12.04 ./gs_project65 > run65.out$$
lb-run ROOT/6.12.04 ./gs_project66 > run66.out$$
lb-run ROOT/6.12.04 ./gs_project67 > run67.out$$
lb-run ROOT/6.12.04 ./gs_project68 > run68.out$$
lb-run ROOT/6.12.04 ./gs_project69 > run69.out$$
lb-run ROOT/6.12.04 ./gs_project610 > run610.out$$

# Create a separate output file for each job submitted
#export workdir=pythia_$$

#mkdir -p $workdir
#cd $workdir

# Run the executable, assume compiled already before submitting to the batch system

#../gs_ex300 > run.out
cd $MYPYTHIA
#./gs_project65 > run.out




