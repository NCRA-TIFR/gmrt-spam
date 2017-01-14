#!/bin/sh
#------------------------------------------------------
# Csh script for setting up the AIPS environment and
# starting AIPS.
# Created on 20050621 by H.T. Intema, Sterrewacht Leiden
# 20080801 HTI - Changed csh script to sh script
# 20120810 HTI - Changed to work at NRAO SOC
# 20130214 HTI - Changed to work with different TVs
# 20141004 HTI - Changed to work on NRAO DSOC cluster
#------------------------------------------------------

# AIPS startup
#if [ -f /home/AIPS/HOSTS.LIST ]
#then
#  if [ -f ${SPAM_PATH}/AIPS/HOSTS.LIST ]
#  then
#    if [ ! -f ${SPAM_PATH}/AIPS/HOSTS.LIST.OLD ]
#    then
#      mv ${SPAM_PATH}/AIP/HOSTS.LIST ${SPAM_PATH}/AIPS/HOSTS.LIST.OLD
#    fi
#  fi
#  cp /home/AIPS/HOSTS.LIST ${SPAM_PATH}/AIPS/
#  echo "+  LOCALHOST   LNX64  SPAM     NONE        cluster node running SPAM" >> /lustre/hintema/spam/AIPS/HOSTS.LIST
#fi

. ${SPAM_PATH}/AIPS/LOGIN.SH
aips_version="TST"
aips_printer=1

# get name and path of this script
scriptname=`basename $0`
scriptpath=`dirname $0`

# print welcome message
echo ${scriptname} ": Startup script for AIPS"

# check if MS directory name is given
if [ $# -lt 1 ]
then
  echo "Usage : " ${scriptname} " working_directory [display_id]/[user_id runfile [arguments]]"
  echo "  working_directory: directory that contains the following AIPS directories"
  echo "    /workx : AIPS working directories where x = 1 to 32"
  echo "    /runfil: directory that contains AIPS runfiles"
  echo "    /prtfil: directory for AIPS printfiles"
  echo "    /fits  : directory that contains fits files"
  echo "  display_id: AIPS display ID to use; 0 = new display, omit for default 1"
  echo "  user_id: AIPS user ID"
  echo "  runfile: name of runfile to run (UPPER CASE, no extension), located in /runfil directory"
  echo "    runfile layout: init, main(user_id,argument1,...), clean"
  echo "  arguments: up to 10 runfile arguments, space separated"
  echo "    put strings between double quotes"
  exit 1
fi


# check max. number of arguments
if [ $# -gt 11 ]
then
  echo "  ERROR: more than 8 arguments for runfile"
  exit 1
fi

# check if working directory exists
if ( ! [ $1 = "." ] )
then
  if ( ! [ -d $1 ] )
  then
    echo "  ERROR: Cannot find working_directory" $1
    exit 1
  fi
fi

# store current directory position and enter working directory
startdir=`pwd -P`
cd $1
currentdir=`pwd -P`

# get display id
display_id=-1
if [ $# = 2 ]
then
  display_id=$2
fi

# get user id and user file extension
# check runfile
if [ $# -gt 2 ]
then
  user_id=$2
  runbase=`basename $3`
  runfile="runfil/"${runbase}".001"
  if ( ! [ -f ${runfile} ] )
  then
    echo "  ERROR: Cannot find runfile" ${currentdir}"/"${runfile}
    cd ${startdir}
    exit 1
  fi
fi

# check DA00 directory and create it if necessary
if ( ! [ -d ${DA00} ] )
then
  mkdir -p ${DA00}
  cp ${AIPS_ROOT}/DA00/${SPAM_HOST}/* ${DA00}
fi

# check and create AIPS working directories
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
  fi
  x=`expr $x + 1`
done
da01_source_dir=${AIPS_ROOT}'/DATA/'${HOST}'_1'
if ( ! [ -f 'work01/MSD001000.001;' ] )
then
  cp ${scriptpath}/'MSD001000.001;' work01
fi

#HOME=.
# create .dadevs file in home directory
if [ -f ${HOME}/.dadevs ]
then
  mv -f ${HOME}/.dadevs ${HOME}/.dadevs_last
fi
touch ${HOME}/.dadevs
ls -1 -d work* > wworkdirs.txt
#workcount=`wc wworkdirs.txt | awk '{print $1}'`
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
if ( ! [ -d $prtfildir ] )
then
  mkdir $prtfildir
fi

# create fits directory
fitsdir='fits'
if ( ! [ -d $fitsdir ] ) then
  mkdir $fitsdir
fi

# create data directory
datfildir='datfil'
if ( ! [ -d $datfildir ] )
then
  mkdir $datfildir
fi

# create python directory
pythondir=${currentdir}/python
if ( ! [ -d $pythondir ] )
then
  mkdir $pythondir
fi

# set AIPS environment
AIPSDIR=${currentdir}; export AIPSDIR
FIT=${AIPSDIR}/fits; export FIT
RUN=${AIPSDIR}/runfil; export RUN
PRT=${AIPSDIR}/prtfil; export PRT
DAT=${AIPSDIR}/datfil; export DAT
AIPS_TEK_EMULATOR=none; export AIPS_TEK_EMULATOR
AIPS_TV_BUFFERED="YES"; export AIPS_TV_BUFFERED

# read in AIPS TV window settings
xrdb -merge ${SPAM_PATH}/bin/Xdefaults_AIPS

#echo "DA00 = " $DA00

# create AIPS input file and start AIPS
if [ $# -gt 2 ]
then
  echo $user_id > input.txt
  echo "run" ${runbase} >> input.txt
  echo "init" >> input.txt
  if [ $# = 4 ]
  then
    echo "main("$user_id","$4")" >> input.txt
  else
    if [ $# = 5 ]
    then
      echo "main("$user_id","$4","$5")" >> input.txt
    else
      if [ $# = 6 ]
      then
        echo "main("$user_id","$4","$5","$6")" >> input.txt
      else
        if [ $# = 7 ]
        then
          echo "main("$user_id","$4","$5","$6","$7")" >> input.txt
        else
          if [ $# = 8 ]
          then
            echo "main("$user_id","$4","$5","$6","$7","$8")" >> input.txt
          else
            if [ $# = 9 ]
            then
              echo "main("$user_id","$4","$5","$6","$7","$8","$9")" >> input.txt
            else
              if [ $# = 10 ]
              then
                echo "main("$user_id","$4","$5","$6","$7","$8","$9","$10")" >> input.txt
              else
                if [ $# = 11 ]
                then
                  echo "main("$user_id","$4","$5","$6","$7","$8","$9","$10","$11")" >> input.txt
                fi
              fi
            fi
          fi
        fi
      fi
    fi
  fi
  echo "clean" >> input.txt
  echo "exit" >> input.txt
  date
  ${AIPS_ROOT}/START_AIPS ${aips_version} tv=local pr=${aips_printer} tpok < input.txt
  date
  rm -f input.txt
else
  if [ ${display_id} -ge 0 ]; then
    ${AIPS_ROOT}/START_AIPS ${aips_version} tv=local:${display_id} pr=${aips_printer} tpok
  else
    ${AIPS_ROOT}/START_AIPS ${aips_version} tv=local pr=${aips_printer} tpok
  fi
fi

# cleanup and exit
cd ${startdir}
echo "Done"

