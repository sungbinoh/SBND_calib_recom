source /cvmfs/larsoft.opensciencegrid.org/products/setup
source /cvmfs/sbn.opensciencegrid.org/products/sbn/setup
setup sbncode v09_75_03_02 -q e20:prof

export CALIB_WORKING_DIR=`pwd`
export DATA_PATH=$CALIB_WORKING_DIR/data/
export PLOT_PATH=$CALIB_WORKING_DIR/output/plots/
export OUTPUTROOT_PATH=$CALIB_WORKING_DIR/output/root/
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$CALIB_WORKING_DIR/include/
