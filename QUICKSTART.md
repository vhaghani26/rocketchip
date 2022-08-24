Rocketchip Quick Start
======================

This is a small demo that shows you how to setup and run Rocketchip. This has
been tested in the following environments.

+ VirutalBox Linux-lite VM (2GB RAM, 1 CPU) - successful
+ VirtualBox Lubuntu 22.04.01
+ Others
+ Mac something or other - failed

## Overview ##

1. Install Conda (if not already installed)
2. Clone Rocketchip Repository
3. Create Rocketchip Conda Enviornment
4. Setup Rocketchip Test Environment
5. Run Rocketchip with Local Data
6. Run Rocketchip with Remote Data


## 1. Install Conda ##

Download a conda installer such as mini-conda and install it with a command
that looks something like the following:

```
sh Miniconda3-py39_4.10.3-Linux-x86_64.sh
```

Also install `mamba` which is better way to run conda installs

```
conda install mamba -n base -c conda-forge
```


## 2. Clone Rocketchip ##

Clone the Rocketchip repository to wherever you typically keep repositories.
For example this might be the `Code` directory in your home.

```
cd ~/Code
git clone https://github.com/vhaghani26/rocketchip
```


## 3. Create Rocketchip Conda Environment ##

Create the conda environment for Rocketchip. This will take a few minutes to
run.

```
cd rocketchip
mamba env create --file rocketchip.yaml
conda activate rocketchip
```


## 4. Setup Rocketchip Test Environment ##

Create a directory for Rocketchip testing. For example, `rocketchip_demo` in
your home directory. This directory will also need subdirectories for `data`,
`local`, and `remote`.

```
cd $HOME
mkdir rocketchip_demo
cd rocket_demo
mkdir data local_demo remote_demo

```

The `data` directory allows multiple projects to share the same files. For
example, all analyses using the human genome will share the same human genome
file in the `data` directory. Similarly, sequencing reads can be shared from
multiple projects. Genomes and reads are _copied_ to the `data` directory.
While this creates some duplication, it prevents errors associated with moving,
editing, or deleting the original files.

The `local_demo` and `remote_demo` directories are used later.

Set the `ROCKETCHIP` environment variable to point to the `rocketchip` GitHub
repository you downloaded previously and `ROCKETCHIP_DATA` to the`data`
directory you just created.

```
ROCKETCHIP=$HOME/Code/rockechip
ROCKETCHIP_DATA=$HOME/rocketchip_demo/data
```

Do these need to be exported?

## 5. Run Rocketchip on Local Data ##




## 6. Run Rocketchip on Remote Data ##








Make sure your Rocketchip conda environment is active. If not, activate it.

```
conda activate rocketchip
```




Stopped here, need to add some smaller genome and read files.



