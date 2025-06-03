# TIMBER {#mainpage}
[Full Documentation](https://lcorcodilos.github.io/TIMBER/)

TIMBER (Tree Interface for Making Binned Events with RDataFrame) is an easy-to-use and fast python analysis framework used to quickly process CMS data sets. 
Default arguments assume the use of the NanoAOD format but any ROOT TTree can be processed.

## Installation instructions for python3

These instructions use python3 and CMSSW. The instructions below are for the cmslpc-el9 cluster for the Run 3 BB -> (b tau tau)(b tau tau) search.

```
source /cvmfs/cms.cern.ch/cmsset_default.sh # this should go in your ~/.bashrc or ~/.bash_profile so it's done automatically when you log in
mkdir nobackup/BBto2b4tau
cd nobackup/BBto2b4tau/
cmsrel CMSSW_13_2_10
cd CMSSW_13_2_10
cmsenv
cd ..
python3 -m virtualenv timber-env
git clone --branch restframes_devel git@github.com:jmhogan/TIMBER.git  ## this requires an "SSH key" for cmslpc-el9. If you don't have it, use the https:// clone method
cd TIMBER/bin/
git clone git@github.com:fmtlib/fmt.git
cd ../../
```

Copy this entire code snippet below to the terminal and hit enter. It will add this information about the boost library paths to the activation script via the "cat" command.

```
cat <<EOT >> timber-env/bin/activate

export BOOSTPATH=/cvmfs/cms.cern.ch/el8_amd64_gcc10/external/boost/1.78.0-0d68c45b1e2660f9d21f29f6d0dbe0a0/lib
if grep -q '\${BOOSTPATH}' <<< '\${LD_LIBRARY_PATH}'
then
  echo 'BOOSTPATH already on LD_LIBRARY_PATH'
else
  export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:\${BOOSTPATH}
  echo 'BOOSTPATH added to PATH'
fi
EOT
```

Now you can activate the python3 environment, set a proper LD_LIBRARY_PATH for boost libraries and build the TIMBER binaries

```
source timber-env/bin/activate
cd TIMBER
source setup.sh
```

## Everyday run setup

You only need to do the installation and compile instructions once. To get back into your TIMBER environment for editing and running our code: 

```
source /cvmfs/cms.cern.ch/cmsset_default.sh # this should go in your ~/.bashrc or ~/.bash_profile so it's done automatically when you log in
voms-proxy-init --voms cms --valid 168:00 # will be valid for a week, only needed when it's needed
cd nobackup/BBto2b4tau/CMSSW_13_2_10
cmsenv
cd ../
source timber-env/bin/activate
cd TIMBER/
python3 BBTo2b4tau.py testfile_SIGNAL_2022.txt 0 0 2022 # args are: file list, first file number, last file number, year
```

Two test files have been made for the "2022" data segment: `testfile_DATA_2022.txt` and `testfile_SIGNAL_2022.txt`.

