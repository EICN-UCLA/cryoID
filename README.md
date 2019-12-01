## cryoID

cryoID is a python-based program that determines the unique identity of the protein(s) in unknown near-atomic resolution cryoEM density maps from a pool of candidate protein sequences. It takes cryoID two sequencial steps to get the job done, with two subprograms called get_queries and search_pool respectively. More technical details are given in our [Nature Methods paper](https://www.nature.com/articles/s41592-019-0637-y). If you find cryoID useful, please cite this paper.

cryoID was designed by Xiaorun Li (Lee), Chi-Min Ho and Mason Lai from Prof. Hong Zhou lab at UCLA, and developed by Xiaorun. We thank Tom Terwilliger for the development of new phenix tools used in cryoID subprogram get_queries.

cryoID is an open-source software under the MIT license.


## Installation

cryoID executables and source codes (python 2.7) are both provided here. Current version was developed for Linux-based systems (CentOS 6 & 7, etc.) with GNU C Library versions >= 2.17. To install the program, download cryoID from github:
```
git clone https://github.com/EICN-UCLA/cryoID.git
```
Make sure you have the privilege to run cryoID executables. If not, change the files' property.
```
cd /your/directory/to/cryoID/bin
chmod +x *
```

cryoID generates query sequences from cryoEM density maps by calling sequence_from_map, a new tool in PHENIX (version >= 1.14-3260. The version used in cryoID tutorial was 1.17.1-3660). If you don't have it, download [PHENIX](https://www.phenix-online.org/) and install one.
```
tar -zxvf phenix-installer-version.tar.gz
cd ./phenix-installer-version
./install --prefix=/your/directory/to/phenix
```

cryoID relies on COOT to open the query model for users' inspection. Most versions will do. If you do not have COOT installed, download it [HERE](https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/binaries/release/). 
```
tar -zxvf coot-version.tar.gz
```

Blastp and makeblastdb from BLAST are used for query searching in cryoID and have been incorperated into our program under cryoID/bin directory. In rare cases if they don't work on your system, you can install your own version [HERE](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).

The last step is to set up your environment for running cryoID. For example, we add the following lines to .bashr file for our bash shell.
```
# .bashrc
export PATH=/your/directory/to/cryoID/bin:/your/directory/to/coot/bin:$PATH
source /path/to/phenix/phenix_env.sh
```

Now you should be able to run cryoID.

If you want to work on the source codes (python 2.7), [Anaconda](https://www.anaconda.com/distribution/) is recommended.


## Usage

Please refer to the [tutorial](https://github.com/EICN-UCLA/cryoID/) for instructions on cryoID usage. The test data can be download [HERE](https://github.com/EICN-UCLA/cryoID_SM/).

