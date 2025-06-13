#!/bin/bash

######################################################
## modXNA.sh                                        ##
## Script to generate modified nucleotides.         ##
######################################################
VERSION='1.8'

# Check for required programs
if [ -z "$CPPTRAJ" ] ; then
  CPPTRAJ=`which cpptraj`
  if [ -z "$CPPTRAJ" ] ; then
    echo -e "  \e[31mCPPTRAJ not found.\e[39m"
    exit 1
  fi
elif [ ! -f "$CPPTRAJ" ] ; then
  echo -e "  \e[31mCPPTRAJ $CPPTRAJ not found.\e[39m"
  exit 1
fi

# CPPTRAJ version check
cc_version=`$CPPTRAJ --version`
if [ -z "$cc_version" ] ; then
  echo -e "  \e[31mError: Could not get version from cpptraj $CPPTRAJ.\e[39m"
  exit 1
fi
#echo "DEBUG: $cc_version"
cc_version=`echo $cc_version | awk 'BEGIN{line=0;}{
  if (line == 0 && $2 == "Version") {
    print substr($3,2);
    exit 0;
  } else {
    # Very old version.
    line++;
  }
  if (line == 1 && $2 == "Trajectory") {
    print substr($4,2);
    exit 0;
  }
}'`
#cc_version=`echo $cc_version | sed -e 's/-/./'`
cc_version_major=`echo $cc_version | cut -d'.' -f1`
cc_version_minor=`echo $cc_version | cut -d'.' -f2`
cc_version_patch=`echo $cc_version | cut -d'.' -f3`
#echo "DEBUG: CPPTRAJ version $cc_version major $cc_version_major minor $cc_version_minor revision $cc_version_patch"
if [ -z "$cc_version_patch" ] ; then
   echo -e "  \e[31mError: CPPTRAJ version from AmberTools $cc_version_major.$cc_version_minor is too old.\e[39m"
   echo -e "  \e[31mError: Requires a more modern CPPTRAJ version, at least 6.26.0\e[39m"
   exit 1
fi
cc_version_ok=1
if [ $cc_version_major -lt 6 ] ; then
  cc_version_ok=0
elif [ $cc_version_major -eq 6 -a $cc_version_minor -lt 26 ] ; then
  cc_version_ok=0
fi
if [ $cc_version_ok -eq 0 ] ; then
  echo -e "  \e[31mError: CPPTRAJ version is too old. Require at least 6.26.0\e[39m"
  exit 1
fi

# Check LEAP
LEAP=`which tleap`
if [ -z "$LEAP" ] ; then
  echo -e "  \e[31mtleap not found.\e[39m"
  exit 1
fi

LEAP=`which sander`
if [ -z "$LEAP" ] ; then
  echo -e "  \e[31mSANDER not found.\e[39m"
  exit 1
fi

# Variables
SUMMARY_CONTENTS=""
TEST=0
LENGTH=1
IS_5CAP=0
IS_3CAP=0

# Library path for the templates
export DAT='/home/ros/work/modxna/git/modXNA'
export DATBB=$DAT'/dat/lib_backbone'
export DATSU=$DAT'/dat/lib_sugar'
export DATBA=$DAT'/dat/lib_base'

# ==============================================================================
Help() {
  echo "Command line options"
  echo " Required:"
  echo "  -i <file>         : Input file (required)"
  echo "  Input file format : <backbone> <sugar> <base>"
  echo ""
  echo " Optional:"
  echo "  -m <name>    : Optional name of generated residue instead of random name"
  echo "  --5cap       : Create a 5'-terminal residue"
  echo "  --3cap       : Create a 3'-terminal residue"
  exit 1 # Exit script after printing help
}

on_exit() {
  echo ""
  echo "   Cleaning up...(remove tmp files, etc)"
  if [ -f tmp.* ] ; then
    rm tmp.*
  fi
}

# ==============================================================================
echo "-------------------------------------------------"
echo "-----               modXNA                  -----"
echo "----- Modified nucleotide residue Generator -----"
echo "-------------------------------------------------"
echo "modXNA.sh Version $VERSION"
echo "CPPTRAJ version $cc_version_major.$cc_version_minor.$cc_version_patch detected."
# ==============================================================================

# Parse command line options
while [ ! -z "$1" ] ; do
  case "$1" in
    '-i'            ) shift ; INPUT=$1 ;;
    '-m'            ) shift ; RESNAME=$1 ;;
    '--5cap'        ) IS_5CAP=1 ;;
    '--3cap'        ) IS_3CAP=1 ;;
    '-h' | '--help' ) Help ; exit 0 ;;
    *               ) echo "Unrecognized command line option: $1" ; exit 1 ;;
  esac
  shift
done

echo "  INPUT FILE: $INPUT"
if [ ! -f "$INPUT" ] ; then
    echo -e "  \e[31mInput file $INPUT not found.\e[39m"
    echo "  Run -h or --help for more information."
    echo "** WARNING - Input file needs to have a RET after last line as of this version **"
  exit 1
fi

# ==============================================================================
# Remove empty lines from input file
sed -i '/^$/d' $INPUT

# Read input file
echo "  Reading input from file:$INPUT"

while read OPTLINE ; do

    ## Ignore comments
    [[ $OPTLINE = \#* ]] && continue

    if [ -z "$RESNAME" ] ; then
        ## Generate random new residue name
        RESNAME0=$(cat /dev/urandom | tr -dc 'a-zA-Z' | head -c 3)
        ## Change residue name to UPPERCASE
        RESNAME=${RESNAME0^^}

    fi

    [[ $OPTLINE = \#* || -z "$OPTLINE" ]] && continue
    SUMMARY_CONTENTS+="  $OPTLINE"$'\n'
    
    echo "  ==========================================================="
    echo "  New residue number :" $LENGTH
    echo "  New residue name is: "$RESNAME
    BACKBONE=`echo "$OPTLINE" | awk '{print $1;}'`
    echo -e "  + \e[32mBackbone:\e[39m" $BACKBONE
    SUGAR=`echo "$OPTLINE" | awk '{print $2;}'`
    echo -e "  + \e[32mSugar:\e[39m" $SUGAR
    BASE=`echo "$OPTLINE" | awk '{print $3;}'`
    echo -e "  + \e[32mBase:\e[39m" $BASE

    ## Check if all fragments file exist
    if [ ! -f "$DATBB/$BACKBONE.mol2" ]; then
	echo -e "  \e[31mBackbone:\e[39m" $BACKBONE " not found."
	exit 1
    elif [ ! -f "$DATSU/$SUGAR.mol2" ]; then
	echo -e "  \e[31mSugar:\e[39m" $SUGAR " not found."
	exit 1
    elif [ ! -f "$DATBA/$BASE.mol2" ]; then
	echo -e "  \e[31mBase:\e[39m" $BASE " not found."
	exit 1
    fi

    #==================================================
    # Reading masks for BASE
    HEAD01BASE=$(grep modXNA $DATBA/$BASE.mol2 |awk -F: '{print $4}')
    HEAD01BASESTRIP=$(grep modXNA $DATBA/$BASE.mol2 |awk -F: '{print $6}')

    # Reading masks for BACKBONE
    HEAD01BACKBONE=$(grep modXNA $DATBB/$BACKBONE.mol2 |awk -F: '{print $4}')
    HEAD01BACKBONESTRIP=$(grep modXNA $DATBB/$BACKBONE.mol2 |awk -F: '{print $6}')
    TAIL01BACKBONE=$(grep modXNA $DATBB/$BACKBONE.mol2 |awk -F: '{print $8}')
    TAIL01BACKBONESTRIP=$(grep modXNA $DATBB/$BACKBONE.mol2 |awk -F: '{print $10}')

    # Reading masks for SUGAR
    HEAD01SUGAR=$(grep modXNA $DATSU/$SUGAR.mol2 |awk -F: '{print $4}')
    HEAD01SUGARSTRIP=$(grep modXNA $DATSU/$SUGAR.mol2 |awk -F: '{print $6}')
    TAIL01SUGAR=$(grep modXNA $DATSU/$SUGAR.mol2 |awk -F: '{print $8}')
    TAIL01SUGARSTRIP=$(grep modXNA $DATSU/$SUGAR.mol2 |awk -F: '{print $10}')
    ANCHOR03SUGAR=$(grep modXNA $DATSU/$SUGAR.mol2 |awk -F: '{print $12}')
    ANCHOR03SUGARSTRIP=$(grep modXNA $DATSU/$SUGAR.mol2 |awk -F: '{print $14}')

    echo "  ==========================================================="
    echo "  + Checking masks from fragments... "    
    echo -e "  + \e[32mBackbone head atom to link:\e[39m" $BACKBONE $HEAD01BACKBONE
    echo -e "  + \e[32mBackbone head atoms to strip:\e[39m" $BACKBONE $HEAD01BACKBONESTRIP    
    echo -e "  + \e[32mBackbone tail atom to link:\e[39m" $BACKBONE $TAIL01BACKBONE
    echo -e "  + \e[32mBackbone tail atoms to strip:\e[39m" $BACKBONE $TAIL01BACKBONESTRIP    
    echo ""
    echo -e "  + \e[32mBase head atom to link:\e[39m" $BASE $HEAD01BASE
    echo -e "  + \e[32mBase head atoms to strip:\e[39m" $BASE $HEAD01BASESTRIP
    echo ""
    echo -e "  + \e[32mSugar head atom to link:\e[39m" $SUGAR $HEAD01SUGAR
    echo -e "  + \e[32mSugar head atoms to strip:\e[39m" $SUGAR $HEAD01SUGARSTRIP
    echo -e "  + \e[32mSugar tail atom to link:\e[39m" $SUGAR $TAIL01SUGAR
    echo -e "  + \e[32mSugar tail atoms to strip:\e[39m" $SUGAR $TAIL01SUGARSTRIP
    echo -e "  + \e[32mSugar anchor atom to link:\e[39m" $SUGAR $ANCHOR03SUGAR
    echo -e "  + \e[32mSugar anchor atoms to strip:\e[39m" $SUGAR $ANCHOR03SUGARSTRIP
    
    sleep 1

    # Clean temp files
    for TF in tmp.base.mol2 tmp.sugar.mol2 tmp.bb.mol2 tmp.o3.sugar.mol2 tmp.o3.c1.sugar.mol2 tmp.o3.c1.h1.sugar.mol2 tmp2.sugar.mol2 tmp2.base.mol2 ; do
      if [ -f "$TF" ] ; then
        rm $TF
      fi
    done

    # Stage each fragment to a temporary file
    cp $DATBA/$BASE.mol2 tmp.base.mol2
    cp $DATSU/$SUGAR.mol2 tmp.sugar.mol2
    cp $DATBB/$BACKBONE.mol2 tmp.bb.mol2

    # Copy O3' charge from backbone to sugar
    if [ $IS_5CAP -eq 1 ]; then
      # 5'-fragment does not have an O3' atom
      VALBO=$(grep "O5' " tmp.bb.mol2 | awk '{print $9}')
    else
      VALBO=$(grep "O3' " tmp.bb.mol2 | awk '{print $9}')
    fi    
    VALSO=$(grep "O3' " tmp.sugar.mol2 | awk '{print $9}')
    echo "Replacing sugar O3' $VALSO with backbone OP3 $VALBO"
    sed "s/$VALSO/$VALBO/" tmp.sugar.mol2 > tmp.o3.sugar.mol2

    # Copy C1' H1' from base to sugar
    VALC1base=`grep "C1' " tmp.base.mol2 | awk '{print $9}'`
    VALC1sugar=`grep "C1' " tmp.o3.sugar.mol2 | awk '{print $9}'`
    sed "s/$VALC1sugar/$VALC1base/" tmp.o3.sugar.mol2 > tmp.o3.c1.sugar.mol2
    echo "Replacing sugar C1' $VALC1sugar with base C1' $VALC1base"
    VALH1base=`grep "H1' " tmp.base.mol2 | awk '{print $9}'`
    VALH1sugar=`grep "H1' " tmp.o3.c1.sugar.mol2 | awk '{print $9}'`
    sed "s/$VALH1sugar/$VALH1base/" tmp.o3.c1.sugar.mol2 > tmp.o3.c1.h1.sugar.mol2
    echo "Replacing sugar H1' $VALH1sugar with base H1' $VALH1base"
    cp tmp.o3.c1.h1.sugar.mol2 tmp.sugar.mol2

    # Determine if any component is modified
    has_modifications=0
    # Is th sugar modified?
    sugar_has_modifications=1
    for sugarname in DC2 RC3 ; do
      if [ "$sugarname" = "$SUGAR" ] ; then
        sugar_has_modifications=0
	break
      fi
    done
    # Is the base modified?
    base_has_modifications=1
    for basename in DAA DCC DGG DTT RAA RCC RUU RGG ; do
      if [ "$basename" = "$BASE" ] ; then
        base_has_modifications=0
        break
      fi
    done
    # If any are modified, all are modified
    if [ $sugar_has_modifications -eq 1 -o $base_has_modifications -eq 1 ] ; then
      has_modifications=1
    fi

    echo "Has mods: $has_modifications  Base has mods: $base_has_modifications  Sugar has mods: $sugar_has_modifications"
    ## ADD correction factors to the sugar
    if [ $sugar_has_modifications -eq 1 -o $base_has_modifications -eq 1 ] ; then
      cat > tmp.cpptraj.sugarmod.in <<EOF
parm tmp.sugar.mol2
loadcrd tmp.sugar.mol2 name sugarfragment
EOF
      if [ $base_has_modifications -eq 1 ] ; then
        echo "Applying chi corrections to sugar"
        cat >> tmp.cpptraj.sugarmod.in <<EOF
## Joint correction terms ###
change crdset sugarfragment charge by -0.4000 of @C1'
change crdset sugarfragment charge by  0.2000 of @H1'
EOF
      fi
      if [ $sugar_has_modifications -eq 1 ] ; then
        echo "Applying JCC to sugar"
        cat >> tmp.cpptraj.sugarmod.in <<EOF
# Avg value that 'works' for both DNA/RNA
change crdset sugarfragment charge by -0.1750 of @C3'
change crdset sugarfragment charge by  0.1750 of @H3'
change crdset sugarfragment charge by -0.2000 of @C5'
change crdset sugarfragment charge by  0.100 of @H5'
change crdset sugarfragment charge by  0.100 of @H5''
EOF
      fi
      cat >> tmp.cpptraj.sugarmod.in <<EOF
crdout sugarfragment tmp2.sugar.mol2
run
EOF
      cpptraj -i tmp.cpptraj.sugarmod.in
      if [ $? -ne 0 ] ; then
        echo "Error: Charge modification of sugar failed."
        exit 1
      fi
#exit 0
    else
      # No sugar modifications
      cp tmp.sugar.mol2 tmp2.sugar.mol2
    fi # END sugar modifications

    ### ADD correction factors to the base
    if [ $base_has_modifications -eq 1 ] ; then
      echo "Applying chi corrections to base"
      ### Charge for N9/N1 (head atom): 0.2
      cpptraj -p tmp.base.mol2 <<EOF
loadcrd tmp.base.mol2 name sugarfragment
### Joint correction terms ###
change crdset sugarfragment charge by 0.2000 of $HEAD01BASE
crdout sugarfragment tmp2.base.mol2
run
EOF
      if [ $? -ne 0 ] ; then
        echo "Error: Charge modification of base failed."
        exit 1
      fi
    else
      # No modification
      cp tmp.base.mol2 tmp2.base.mol2
    fi # END base modifications

    ## Strip fragments from capping groups
    ## Adjust the charge for each fragment
    cat > tmp.strip.cpptraj<<EOF
### BACKBONE
parm tmp.bb.mol2 name backbone
trajin tmp.bb.mol2 parm backbone
EOF

    if [ $IS_5CAP -eq 0 ]; then
      echo "strip $HEAD01BACKBONESTRIP" >> tmp.strip.cpptraj
    fi

    cat >> tmp.strip.cpptraj<<EOF
strip $TAIL01BACKBONESTRIP charge -0.8832
trajout tmp.backbone-striped.mol2 mol2
run
clear all
parm tmp2.sugar.mol2 name sugar
trajin tmp2.sugar.mol2 parm sugar
EOF

    if [ $IS_3CAP -eq 0 ]; then	
	echo "strip $TAIL01SUGARSTRIP" >> tmp.strip.cpptraj
    fi
    
    cat >> tmp.strip.cpptraj<<EOF
strip $ANCHOR03SUGARSTRIP
strip $HEAD01SUGARSTRIP charge -0.01191
trajout tmp.sugar-striped.mol2 mol2
run
clear all
parm tmp2.base.mol2 name base
trajin tmp2.base.mol2 parm base
strip $HEAD01BASESTRIP charge -0.10489
trajout tmp.base-striped.mol2 mol2
EOF
    
    ## Run CPPTRAJ, create stripped backbone and sugar
    cpptraj -i tmp.strip.cpptraj
    if [ $? -ne 0 ] ; then
      echo "Error: Creation of stripped backbone and sugar failed."
      exit 1
    fi
    
    ## Combine backbone and sugar fragments
    cat > tmp.combine.cpptraj<<EOF
parm tmp.sugar-striped.mol2
loadcrd tmp.sugar-striped.mol2 name Sugar parm tmp.sugar-striped.mol2
parm tmp.backbone-striped.mol2
loadcrd tmp.backbone-striped.mol2 name Backbone parm tmp.backbone-striped.mol2
dataset connect Sugar headmask $HEAD01SUGAR
dataset connect Backbone tailmask $TAIL01BACKBONE
sequence Backbone Sugar name BackboneSugar
change crdset BackboneSugar mergeres firstres 1 lastres 2
change crdset BackboneSugar resname from * to $RESNAME
change crdset BackboneSugar oresnums of :1 min 1 max 1
crdout BackboneSugar tmp.BackboneSugar.mol2
EOF
    ## Run CPPTRAJ, create sugar+base
    cpptraj -i tmp.combine.cpptraj
    if [ $? -ne 0 ] ; then
      echo "Error: Creation of stripped sugar+base failed."
      exit 1
    fi
    
    ## Combine BackboneSugar and base
    if [ -f 'tmp.Nucleotide.mol2' ] ; then
      rm tmp.Nucleotide.mol2
    fi
    cat > tmp.combine.cpptraj<<EOF
parm tmp.BackboneSugar.mol2
loadcrd tmp.BackboneSugar.mol2 name BackboneSugar parm tmp.BackboneSugar.mol2
parm tmp.base-striped.mol2
loadcrd tmp.base-striped.mol2 name Base parm tmp.base-striped.mol2
dataset connect BackboneSugar headmask $ANCHOR03SUGAR
dataset connect Base tailmask $HEAD01BASE
sequence Base BackboneSugar name Nucleotide
change crdset Nucleotide mergeres firstres 1 lastres 2
change crdset Nucleotide resname from * to $RESNAME
change crdset Nucleotide oresnums of :1 min 1 max 1
crdout Nucleotide tmp.Nucleotide.mol2
EOF
    cpptraj -i tmp.combine.cpptraj
    if [ $? -ne 0 ] ; then
      echo "Error: Creation of nucleotide failed."
      exit 1
    fi
    
    ###########################################
    ## OPTIMIZATION
    ############################################
    ## From the Nucleotide.mol2 file, generate
    ## topology/coordinates to run SANDER
    ## to optimize/clean the structure
    cat > tmp.opt.tleap<<EOF
###########################################
source leaprc.protein.ff19SB
source leaprc.DNA.OL15
#source leaprc.RNA.OL3
#source leaprc.gaff2
###########################################
addpath $DAT/dat
loadamberparams frcmod.modxna
###########################################
x = loadmol2 tmp.Nucleotide.mol2
charge x
check x
saveamberparm x tmp.opt.topo tmp.opt.coords
quit
EOF
    tleap -s -f tmp.opt.tleap

    ## Run 2000 frames of minimization
    cat > tmp.opt.in<<EOF
energy minimization
 &cntrl
  imin=1,maxcyc=2000,ntx=1,ntwr=100,ntpr=10,
  ioutfm=0,ntxo=1,cut=1000.0,ntb=0,igb=5,
 &end
EOF

    echo "   > Running GB optimization"
    sander -O -i tmp.opt.in -p tmp.opt.topo -c tmp.opt.coords -r tmp.opt.ncrst

    ## Generate check files from the optimization
    cpptraj -p tmp.opt.topo -y tmp.opt.ncrst -x tmp.opt.pdb
    cpptraj -p tmp.opt.topo -y tmp.opt.ncrst -x tmp.opt.mol2

    ## Create LIB file from optimized structure
    cat > tmp.lib.tleap<<EOF
###########################################
source leaprc.protein.ff19SB
source leaprc.DNA.OL15
###########################################
addpath $DAT/dat
loadamberparams frcmod.modxna
###########################################
$RESNAME = loadmol2 tmp.opt.mol2
set $RESNAME restype nucleic
set $RESNAME name $RESNAME
EOF

    if [ $IS_5CAP -eq 1 ]; then
      echo "set $RESNAME head 0" >> tmp.lib.tleap
      HEAD_ATOM="None"
    else
      echo "set $RESNAME head $RESNAME.1.P" >> tmp.lib.tleap
      HEAD_ATOM="P"
    fi

    if [ $IS_3CAP -eq 1 ]; then
      echo "set $RESNAME tail 0" >> tmp.lib.tleap
      TAIL_ATOM="None"
    else
      echo "set $RESNAME tail $RESNAME.1.O3'" >> tmp.lib.tleap
      TAIL_ATOM="O3'"
    fi

    echo "saveoff $RESNAME $RESNAME.lib" >> tmp.lib.tleap
    echo "quit" >> tmp.lib.tleap

    tleap -s -f tmp.lib.tleap
    
done<$INPUT

# Extract information for the summary
libfile="${RESNAME}.lib"
HEAD_ID="none"
TAIL_ID="none"

if [[ -f "$libfile" ]]; then
    mapfile -t connect_lines < <(awk "/!entry.${RESNAME}.unit.connect array int/{getline; head=\$1; getline; tail=\$1; print head; print tail}" "$libfile")
    HEAD_ID="${connect_lines[0]}"
    TAIL_ID="${connect_lines[1]}"
fi

echo
echo "====== modXNA Summary ======"
echo "Input file:     $INPUT"
echo "Residue name:   $RESNAME"
echo "5' capped:      $IS_5CAP"
echo "3' capped:      $IS_3CAP"
echo "HEAD atom:      $HEAD_ATOM (atom #: $HEAD_ID)"
echo "TAIL atom:      $TAIL_ATOM (atom #: $TAIL_ID)"
echo "----------------------------"
echo "Input content:"
echo "$SUMMARY_CONTENTS"
echo "============================"

exit 0
