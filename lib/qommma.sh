#!/bin/bash
Indir=/home/neil/quick_jobs/crambin/
cd $Indir

export TMPDIR=/temp0/neil/qmmm
export GAUSS_SCRDIR=$TMPDIR
export LD_LIBRARY_PATH=/usr/local/gcc-6.3.0/lib64

if [ ! -d $GAUSS_SCRDIR ]
then
mkdir -p $GAUSS_SCRDIR
fi

export g09root=/usr/local/chem/g09D01
bash $g09root/bsd/g09.profile
export GAUSS_EXEDIR=$g09root

/home/neil/software_and_scripts/qommma_8.07/qommma $Indir/qommma.in


