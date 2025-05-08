export CALIB_WORKING_DIR=`pwd`
export DATA_PATH=$CALIB_WORKING_DIR/data/
export PLOT_PATH=$CALIB_WORKING_DIR/output/plots/
export OUTPUTROOT_PATH=$CALIB_WORKING_DIR/output/root/
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$CALIB_WORKING_DIR/include/
source $CALIB_WORKING_DIR/bin/BashColorSets.sh

#####################################################################################
## -- Host dependent settings

#### -- Setup root
MY_OS_REL=$(cat /etc/os-release | grep ^NAME | sed -e 's/NAME=//g' -e 's/"//g')
if [[ "$MY_OS_REL" == "AlmaLinux" && $(hostname) != *"dune-gpu01"* ]]; then
  source /cvmfs/larsoft.opensciencegrid.org/spack-v0.22.0-fermi/setup-env.sh
  spack load root@6.28.12 arch=linux-almalinux9-x86_64_v3
elif [ "$MY_OS_REL" = "Scientific Linux" ]; then

  #Check if PRODUCTS is undefined -- if so, set up relevant ups area
  if [[ -z $PRODUCTS ]]; then
    if [[ $HOSTNAME == "sbnd"* ]]; then
      echo "SBND"
      MY_EXPERIMENT="sbnd"
    elif [[ $HOSTNAME == "icarusgpvm"* ]]; then
      echo "ICARUS"
      MY_EXPERIMENT="sbnd"
    else
      echo "Warning: Unrecognized hostname $HOSTNAME"
    fi
  fi

  source /cvmfs/${MY_EXPERIMENT}.opensciencegrid.org/products/${MY_EXPERIMENT}/setup_${MY_EXPERIMENT}.sh
  setup root v6_28_12 -q e26:p3915:prof 
  setup xrootd v5_5_5a -q e26:p3915:prof
  setup cmake v3_27_4
else
  echo "WARNING: Seems you are using a private machine to run this repo"
  echo "I do not automatically set up ROOT. If ROOT is already setup, it shoould be okay"
fi

#### -- Calib ntuple list dir
export SAMPLE_PATH=$DATA_PATH/sample_list/sungbinosx/
if [[ `hostname` == *"sbnd"* ]]
then
    source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh
    export SAMPLE_PATH=$DATA_PATH/sample_list/sbndgpvm/
    export SBND_DATA_PATH=/cvmfs/sbnd.opensciencegrid.org/products/sbnd/sbnd_data/
    export SBNDDATA_VERSION=v01_28_00
fi

if [[ `hostname` == *"dune-gpu01"* ]]
then
    export SAMPLE_PATH=$DATA_PATH/sample_list/dune-gpu01/
    export FILELIST_LABEL=_dune_gpu01_run_
fi
#####################################################################################
