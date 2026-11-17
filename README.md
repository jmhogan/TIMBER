# TIMBER {#mainpage}
TIMBER (Tree Interface for Making Binned Events with RDataFrame) is an easy-to-use and fast python analysis framework used to quickly process CMS data sets. 
Default arguments assume the use of the NanoAOD format but any ROOT TTree can be processed.

## Installation instructions for python3

These instructions use python3 and CMSSW. The instructions below have been tested on el8 and el9, both lxplus and lpc.

```
source /cvmfs/cms.cern.ch/cmsset_default.sh # this should go in your ~/.bashrc or ~/.bash_profile so it's done automatically when you log in
mkdir nobackup/BBto2b4tau
cd nobackup/BBto2b4tau/
cmsrel CMSSW_14_1_9
cd CMSSW_14_1_9
cmsenv
pip3 install --user --no-binary=correctionlib --upgrade correctionlib
cd ..
python3 -m virtualenv timber-env #If this step fails, you might need to do `python3 -m pip install --user virtualenv`
git clone git@github.com:jmhogan/TIMBER.git # this requires an "SSH key" for cmslpc-el9. If you don't have it, use the https:// clone method
cd TIMBER/
mkdir bin
cd bin
git clone git@github.com:fmtlib/fmt.git
cd ../..
```

Boost library path (the boost version as well!) may change depending on the CMSSW version so this may need to be modified by hand.

Copy the whole multi-line string to the environment activation script by copying this cell and executing it on the command line. NOTE: it will NOT WORK to copy and paste these contents into your activate script using a text editor!

```
cat <<EOT >> timber-env/bin/activate
export SCRAM_ARCH=${SCRAM_ARCH}
if [[ "\$SCRAM_ARCH" == "el8_amd64_gcc11" ]]; then
  BOOSTPATH=/cvmfs/cms.cern.ch/el8_amd64_gcc11/external/boost/1.78.0-dfb1dc972d1e1af822bb548909730506/lib
elif [[ "\$SCRAM_ARCH" == "el9_amd64_gcc11" ]]; then
  BOOSTPATH=/cvmfs/cms.cern.ch/el9_amd64_gcc11/external/boost/1.78.0-c49033d06e1a3bf1beac1c01e1ef27d6/lib
elif [[ "\$SCRAM_ARCH" == "el9_amd64_gcc12" ]]; then
  BOOSTPATH=/cvmfs/cms.cern.ch/el9_amd64_gcc12/external/boost/1.80.0-e4b2127ad40e46379cea2f9087beeda8/lib
else
  BOOSTPATH=/cvmfs/cms.cern.ch/el8_amd64_gcc10/external/boost/1.78.0-0d68c45b1e2660f9d21f29f6d0dbe0a0/lib
fi

if [[ ":\$LD_LIBRARY_PATH:" != *":\$BOOSTPATH:"* ]]; then
  export LD_LIBRARY_PATH="\${LD_LIBRARY_PATH:+\$LD_LIBRARY_PATH:}\$BOOSTPATH"
  echo "BOOSTPATH added to LD_LIBRARY_PATH"
else
  echo "BOOSTPATH already on LD_LIBRARY_PATH"
fi

EOT
```

This will activate the python3 environment, set a proper LD_LIBRARY_PATH for boost libraries and build the TIMBER binaries

```
source timber-env/bin/activate
cd TIMBER
source setup.sh
```

Tip: Add the lines below to the top of `timber-env/bin/activate` script. With this, one can skip doing `cmsenv` every time after opening a new shell and just activate the environment instead.
```
cd CMSSW_14_1_9
cmsenv
cd ..
```

## Everyday run setup

You only need to do the installation and compile instructions once. To get back into your TIMBER environment for editing and running our code: 

```
voms-proxy-init --voms cms --valid 168:00 # will be valid for a week, only need to run when the proxy's time is up
cd nobackup/BBto2b4tau/
source timber-env/bin/activate
cd TIMBER/
python3 BBTo2b4tau.py condor/NanoList/Bprime_M1000_2023NanoList.txt 0 0 2023 # args: file list, first file number, last file number, year
```

Many text files live in the `condor/NanoList` folder, for signal, data, and background simulations. Any will work to test the script interactively.

