if [[ `hostname` == *"gpvm"* ]]
then
    source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh
    spack load root@6.28.12
fi

export CALIB_WORKING_DIR=`pwd`
export DATA_PATH=$CALIB_WORKING_DIR/data/
export PLOT_PATH=$CALIB_WORKING_DIR/output/plots/
export OUTPUTROOT_PATH=$CALIB_WORKING_DIR/output/root/
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$CALIB_WORKING_DIR/include/

source $CALIB_WORKING_DIR/bin/BashColorSets.sh
