import os,sys,time
import datetime
import argparse
import warnings

parser = argparse.ArgumentParser(
    description="SBND e-lifetime grid submission command: process input calilb ntuples and produce output root files with dQdx vs. t_drift.",
    epilog="""\
Examples:

  -- 
  $ python run_df_maker.py -ngrid 2 -c ./configs/cohpi_slcdf.py -o test_cohpi_slcdf -i input_0.root,input_1.root,...

  -- Note!!
  Output root files are sent to /pnfs/<exp>/scratch/users/<User>/sbnd_calib_out
""",
    formatter_class=argparse.RawTextHelpFormatter  # Ensures line breaks are preserved
)

parser.add_argument('-o', dest='output', default="", help="output root file name prefix")
parser.add_argument('-i', dest='inputfiles', default="", help="input root file path, you can submit multiple files using comma, i.e.) -i input_0.root,input_1.root")
parser.add_argument('-l', dest='inputfilelist', default="", help="a file of list for input root files")
parser.add_argument('-ngrid', dest='NGridJobs', default=0, type=int, help="Number of grid jobs. Default = 0, no grid submission.")
parser.add_argument('-nfile', dest='NFiles', default=0, type=int, help="Number of files to run. Default = 0, run all input files.")
args = parser.parse_args()


def run_grid(inputfiles):
    # 1) dir/file name style
    JobStartTime = datetime.datetime.now()
    timestamp =  JobStartTime.strftime('%Y_%m_%d_%H%M%S')

    # 2) Define MasterJobDir -- produce grid job submission scripts in $SBNDCALIB_GRID_OUT_DIR
    SBNDCALIB_GRID_OUT_DIR = os.environ['SBNDCALIB_GRID_OUT_DIR']
    MasterJobDir = SBNDCALIB_GRID_OUT_DIR + "/logs/" + args.output + "__" + timestamp + "_log"
    OutputDir = SBNDCALIB_GRID_OUT_DIR + "/etau/" + args.output + "__" + timestamp
    os.system('mkdir -p ' + MasterJobDir)

    # 3) grid job is based on number of files
    ngrid = args.NGridJobs
    if(len(inputfiles) <= ngrid):
        ngrid = len(inputfiles)

    NInputfiles = len(inputfiles)
    print("Number of Grid Jobs: %d, number of input calib ntuple files: %d" % (ngrid, NInputfiles))

    # 4) prepare bash scripts for each job and make tarball
    flistForEachJob = []
    for i in range(0,ngrid):
        flistForEachJob.append( [] )

    for i_line in range(0,len(inputfiles)):
        flistForEachJob[i_line%ngrid].append(inputfiles[i_line])

    for i_flist in range(0,len(flistForEachJob)):
        flist = flistForEachJob[i_flist]
        this_list = open(MasterJobDir + '/input_list_%s.txt'%(i_flist),'w')
        for i_f in range(0,len(flist)):
            this_list.write('%s\n'%(flist[i_f]))
        this_list.close()
        
    os.system('cp ./grid/bin/grid_executable_etau.sh %s' %MasterJobDir)

    # 5) copy scripts for running run_lifetime_loop.C
    CALIB_WORKING_DIR = os.environ['CALIB_WORKING_DIR']
    cp_include = "cp -r " + CALIB_WORKING_DIR + "/include " + MasterJobDir
    cp_setup = "cp " + CALIB_WORKING_DIR + "/setup.sh " + MasterJobDir
    cp_BashColorSets = "cp -r " + CALIB_WORKING_DIR + "/bin " + MasterJobDir
    cp_lifetime_loop = "cp " + CALIB_WORKING_DIR + "/src/run_lifetime_loop.C " + MasterJobDir
    os.system(cp_include)
    os.system(cp_setup)
    os.system(cp_BashColorSets)
    os.system(cp_lifetime_loop)

    yzunif_map_dir = os.environ['SBND_YZCORR_PATH']
    cp_yzunif_map = "cp " + yzunif_map_dir + "/*.root " + MasterJobDir
    os.system(cp_yzunif_map)
    
    os.chdir(MasterJobDir)
    tar_cmd = 'tar cf bin_dir.tar ./'
    os.system(tar_cmd)

    submitCMD = '''jobsub_submit \\
-G sbnd \\
--auth-methods="token" \\
-e LC_ALL=C \\
--role=Analysis \\
--resource-provides="usage_model=DEDICATED,OPPORTUNISTIC" \\
--lines '+FERMIHTC_AutoRelease=True' --lines '+FERMIHTC_GraceMemory=1000' --lines '+FERMIHTC_GraceLifetime=3600' \\
--append_condor_requirements='(TARGET.HAS_SINGULARITY=?=true)' \\
--tar_file_name "dropbox://$(pwd)/bin_dir.tar" \\
--email-to sungbin.oh555@gmail.com \\
-N %d \\
--disk 100GB \\
--expected-lifetime 10h \\
"file://$(pwd)/grid_executable_etau.sh" \\
"%s" \\
"%s"'''%(ngrid,OutputDir,args.output)

    print(submitCMD)
    os.system(submitCMD)
    
    # go back to working dir
    os.chdir(CALIB_WORKING_DIR)

if __name__ == "__main__":
    printhelp = ((args.inputfiles == "" and args.inputfilelist == "") or args.output == "")
    if printhelp:
        parser.print_help()
        print(parser.epilog)
        sys.exit(1)
        
    else:
        ### Organize input list 
        InputSamples = []
        StringForHash = ""
        if args.inputfilelist != "":
            lines = open(args.inputfilelist)
            for line in lines:
                if "#" in line:
                    continue
                line = line.strip('\n')
                InputSamples.append(line)
                StringForHash += line
        else:
            split_inputfiles = args.inputfiles.split(",")
            for split_inputfile in split_inputfiles:
                InputSamples.append(split_inputfile)
                StringForHash += args.inputfiles

        if(args.NFiles > 0 and len(InputSamples) > args.NFiles):
            InputSamples = InputSamples[:args.NFiles]
                
        ### check if it is grid mode for pool mode
        if args.NGridJobs == 0:
            print("0 grid jobs, not running");

        elif args.NGridJobs > 0:
            print("Runing Grid mode");
            run_grid(InputSamples)
            
        else:
            print("-ngrid must be greater than 0 to run.");
