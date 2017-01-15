#!/bin/sh
#------------------------------------------------------
# Script for setting up the AIPS environment and
# starting ParselTongue.
# Created on 20050929 by H.T. Intema, Sterrewacht Leiden
# 20061122 HTI - Adapted script to work with global
#                ParselTongue installation.
# 20080728 HTI - Changed to sh to prevent CTRL-C kill
#------------------------------------------------------

OBIT_PYTHONPATH=${SPAM_PATH}/Obit/python
OBIT_PATH=${SPAM_PATH}/Obit/bin
PT_PYTHONPATH=${SPAM_PATH}/ParselTongue/share/parseltongue/python
PT_PATH=${SPAM_PATH}/ParselTongue/bin

PYTHONPATH=${PT_PYTHONPATH}:${OBIT_PYTHONPATH}:${PYTHONPATH}
export PYTHONPATH
PATH=${PT_PATH}:${OBIT_PATH}:${PATH}
export PATH
. ${SPAM_PATH}/AIPS/LOGIN.SH


# get name and path of this script
scriptname=`basename $0`
scriptpath=`dirname $0`

# print welcome message
if [ $# -le 2 ]
then
  echo ${scriptname} ": Startup script for ParselTongue"
fi

# check if MS directory name is given
if [ $# -lt 2 ]
then
  echo "Usage : " ${scriptname} " working_directory aips_user_id [script.py]"
  echo "  working_directory: directory that contains the following AIPS directories"
  echo "    /workx : AIPS working directories where x = 01 to 32"
  echo "    /runfil: directory that contains AIPS runfiles"
  echo "    /prtfil: directory for AIPS printfiles"
  echo "    /fits  : directory that contains fits files"
  echo "  aips_user_id: AIPS user ID"
  echo "  script.py: name of python script (path relative to working dir)"
  exit 1
fi

# check max. number of arguments
if [ $# -gt 10 ]
then
  echo "  ERROR: too many arguments"
  exit 1
fi

# check if working directory exists
if ( ! [ $1 = "." ] )
then
  if ( ! [ -d $1 ] )
  then
#    echo "  ERROR: Cannot find working_directory" $1
#    exit 1
    mkdir -p $1
  fi
fi

# set AIPS userno in system variable
AIPS_USERID=$2
export AIPS_USERID

# store current directory position and enter working directory
startdir=`pwd -P`
cd $1
currentdir=`pwd -P`

# convert user ID to ehex format
ehex=`${scriptpath}/ehex.py $2`
cadname="CAD000000."${ehex}";"

# check DA00 directory
da00_dir="da00"
if ( ! [ -d ${da00_dir} ] )
then
  mkdir -p ${da00_dir}
  cp -f ${AIPS_ROOT}/DA00/${SPAM_HOST}/* ${da00_dir}
fi

# check AIPS working directories
x=1
x_max=35
while [ $x -le $x_max ]
do
  if [ $x -ge 10 ]
  then
    workdir="work"$x
  else
    workdir="work0"$x
  fi
  if ( ! [ -d ${workdir} ] )
  then
    mkdir ${workdir}
    touch ${workdir}/SPACE
    cp ${scriptpath}/'CAD000000.001;' ${workdir}/
    cp ${scriptpath}/'CAD000000.001;' ${workdir}/${cadname}
  else
    if ( ! [ -f ${workdir}/'CAD000000.001;' ] )
    then
      cp ${scriptpath}/'CAD000000.001;' ${workdir}/
    fi
    if ( ! [ -f ${workdir}/${cadname} ] )
    then
      cp ${scriptpath}/'CAD000000.001;' ${workdir}/${cadname}
    fi
  fi
  x=`expr $x + 1`
done
if ( ! [ -f 'work01/MSD001000.001;' ] )
then
  cp ${scriptpath}/'MSD001000.001;' work01
fi

# create .dadevs file in home directory
if [ -f ${HOME}/.dadevs ]
then
  mv -f ${HOME}/.dadevs ${HOME}/.dadevs_last
fi
touch ${HOME}/.dadevs
ls -1 -d work* > wworkdirs.txt
workcount=`wc wworkdirs.txt | awk '{print $1}'`
cat wworkdirs.txt | while read workdir
do
  echo "+  "${currentdir}"/"${workdir} >> ${HOME}/.dadevs
done
rm -f wworkdirs.txt

# create runfile directory
runfildir='runfil'
if ( ! [ -d ${runfildir} ] )
then
  mkdir ${runfildir}
fi

# create prtfile directory
prtfildir='prtfil'
if ( ! [ -d ${prtfildir} ] )
then
  mkdir ${prtfildir}
fi

# create fits directory
fitsdir='fits'
if ( ! [ -d $fitsdir ] ) then
  mkdir $fitsdir
fi

# create data directory
datfildir='datfil'
if ( ! [ -d ${datfildir} ] )
then
  mkdir ${datfildir}
fi

# create python directory
pythondir=${currentdir}/python
if ( ! [ -d ${pythondir} ] )
then
  mkdir ${pythondir}
fi

# set AIPS environment
AIPS_TEK_EMULATOR=none; export AIPS_TEK_EMULATOR
AIPSDIR=${currentdir}; export AIPSDIR
FIT=${AIPSDIR}/fits; export FIT
FITS=${AIPSDIR}/fits; export FITS
RUN=${AIPSDIR}/runfil; export RUN
PRT=${AIPSDIR}/prtfil; export PRT
DAT=${AIPSDIR}/datfil; export DAT
DA00=${AIPSDIR}/da00; export DA00

DA01=${currentdir}/work01; export DA01
DA02=${currentdir}/work02; export DA02
DA03=${currentdir}/work03; export DA03
DA04=${currentdir}/work04; export DA04
DA05=${currentdir}/work05; export DA05
DA06=${currentdir}/work06; export DA06
DA07=${currentdir}/work07; export DA07
DA08=${currentdir}/work08; export DA08
DA09=${currentdir}/work09; export DA09
DA0A=${currentdir}/work10; export DA0A
DA0B=${currentdir}/work11; export DA0B
DA0C=${currentdir}/work12; export DA0C
DA0D=${currentdir}/work13; export DA0D
DA0E=${currentdir}/work14; export DA0E
DA0F=${currentdir}/work15; export DA0F
DA0G=${currentdir}/work16; export DA0G
DA0H=${currentdir}/work17; export DA0H
DA0I=${currentdir}/work18; export DA0I
DA0J=${currentdir}/work19; export DA0J
DA0K=${currentdir}/work20; export DA0K
DA0L=${currentdir}/work21; export DA0L
DA0M=${currentdir}/work22; export DA0M
DA0N=${currentdir}/work23; export DA0N
DA0O=${currentdir}/work24; export DA0O
DA0P=${currentdir}/work25; export DA0P
DA0Q=${currentdir}/work26; export DA0Q
DA0R=${currentdir}/work27; export DA0R
DA0S=${currentdir}/work28; export DA0S
DA0T=${currentdir}/work29; export DA0T
DA0U=${currentdir}/work30; export DA0U
DA0V=${currentdir}/work31; export DA0V
DA0W=${currentdir}/work32; export DA0W
DA0X=${currentdir}/work33; export DA0X
DA0Y=${currentdir}/work34; export DA0Y
DA0Z=${currentdir}/work35; export DA0Z
DA10=${currentdir}/work10; export DA10
DA11=${currentdir}/work11; export DA11
DA12=${currentdir}/work12; export DA12
DA13=${currentdir}/work13; export DA13
DA14=${currentdir}/work14; export DA14
DA15=${currentdir}/work15; export DA15
DA16=${currentdir}/work16; export DA16
DA17=${currentdir}/work17; export DA17
DA18=${currentdir}/work18; export DA18
DA19=${currentdir}/work19; export DA19
DA20=${currentdir}/work20; export DA20
DA21=${currentdir}/work21; export DA21
DA22=${currentdir}/work22; export DA22
DA23=${currentdir}/work23; export DA23
DA24=${currentdir}/work24; export DA24
DA25=${currentdir}/work25; export DA25
DA26=${currentdir}/work26; export DA26
DA27=${currentdir}/work27; export DA27
DA28=${currentdir}/work28; export DA28
DA29=${currentdir}/work29; export DA29
DA30=${currentdir}/work30; export DA30
DA31=${currentdir}/work31; export DA31
DA32=${currentdir}/work32; export DA32
DA33=${currentdir}/work33; export DA33
DA34=${currentdir}/work34; export DA34
DA35=${currentdir}/work35; export DA35

# add directories to python search path
PYTHONPATH=${pythondir}:${PYTHONPATH}; export PYTHONPATH

# start parseltongue
if [ $# -ge 3 ]
then
  if [ $3 == "debug" ]
  then
    . ${AIPS_VERSION}/SYSTEM/UNIX/DADEVS.SH
#    . ${AIPS_VERSION}/SYSTEM/UNIX/PRDEVS.SH
#    OBIT_PYTHONPATH=${SPAM_PATH}/Obit/python
#    PYTHONPATH=$PYTHONPATH:${SPAM_PATH}/Obit/python
#    export PYTHONPATH
    PYTHONSTARTUP=${SPAM_PATH}/ParselTongue/share/parseltongue/python/ParselTongue.py
    export PYTHONSTARTUP
    gdb ${PYTHON}
#  else
#    if [ $3 == "test" ]
#    then
#      . ${AIPS_VERSION}/SYSTEM/UNIX/DADEVS.SH
#      . ${AIPS_VERSION}/SYSTEM/UNIX/PRDEVS.SH
##      OBIT_PYTHONPATH=${SPAM_PATH}/Obit/python
##      PYTHONPATH=$PYTHONPATH:${OBIT_PYTHONPATH}
##      export PYTHONPATH
#      PYTHONSTARTUP=${PT_PYTHONPATH}/ParselTongue.py
#      export PYTHONSTARTUP
#      ${PYTHON}
  else
#    ParselTongue -V $3 $4 $5 $6 $7 $8 $9 $10
  PYTHONSTARTUP=${PT_PYTHONPATH}/ParselTongue.py
  export PYTHONSTARTUP
#  echo ${PYTHON} ${PYTHONSTARTUP} 
  ${PYTHON} $3 $4 $5 $6 $7 $8 $9 $10
  fi
else
#  ParselTongue
  . ${AIPS_VERSION}/SYSTEM/UNIX/DADEVS.SH
  PYTHONSTARTUP=${PT_PYTHONPATH}/ParselTongue.py
  export PYTHONSTARTUP
  echo ${PYTHON} ${PYTHONSTARTUP} 
	${PYTHON} -c 'print(2)'
fi

# cleanup and exit
if [ -d ${HOME}/.matplotlib/tex-cache ]
then
  rm -f ${HOME}/.matplotlib/tex-cache/*
fi

cd ${startdir}

if [ $# -le 2 ]
then
  echo "Done"
fi

