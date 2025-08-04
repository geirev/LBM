# Lattice-Boltzmann-Model
This repository contains a 3D-implementation of a Lattice-Boltzmann-Model on a D3Q27 lattice for high Reynolds number flow.
The code is NVIDIA's CUDA Fortran and it runs on either single core, OPEN-MP multicore, or GPU.
The collision operator is a single relaxation time third order expansion with regularization similar to [Jacob et al (2018)](https://hal.science/hal-02114308)
and [Feng et al, (2018)](https://doi.org/10.1029/2020MS002107),
but with out the hybrid regularization, which I have not yet needed.
The turbulence closure scheme is the one of [Vreman (2004)](https://doi.org/10.1063/1.1785131).

The model boundary conditions are periodic or inflow-outflow in the i-direction, periodic or closed no-slip or free-slip two timestep bounceback in the j- and k- direction.

The code allows for inserting solid bodies within the model domain to simulated. e.g. flow around an airfoil or a cylinder.

Additionally, there is a complete implementation of an actuator line model for the NREL-5Mw wind turbine, and it is possible to include multiple turbines at any location of the model domain.

Inflow turbulence is mimiked or introduced at a section inside the inflow boundary at i=1, (typically at the slice i=10) by applying a smooth in space and time pseudo random force on the fluid.



<p align="center">
  <a href="#installation">Installation</a> *
  <a href="#code-standards">Code standards</a> *
  <a href="#setting-up-an-experiment">Experiment setup</a> *
  <a href="#plotting">Plotting</a> *
  <a href="#git-instructions">Git instructions</a> *
  <a href="https://github.com/geirev/LBM/blob/master/LICENSE">License</a>
</p>


<p align="center">
<img src="docs/cyl.png" width="300">
</p>

---

# Installation:

## 1. Building the Project

If you plan to collaborate or contribute anything to the project, use the <a href="#1b-advanced-installation">Advanced Installation</a> option.

### 1a. Basic installation

Create a directory to clone the three following repositories:

```bash
git clone git@github.com:geirev/LBM.git
```

### 1b. Advanced installation

Make a personal github account unless you already have one.
Fork the LBM repository.
Next clone the forked repositories and set upstream to the original repositories where
you need to replace <userid> with your github userid

```bash
git clone git@github.com:<userid>/LBM.git
pushd LBM
git remote add upstream https://github.com/geirev/LBM
```
or, if you have set up git-ssh
```bash
git clone git@github.com:<userid>/LBM.git
pushd LBM
git remote add upstream git://github.com:geirev/LBM
popd
```

If you are new to Git, read the section <a href="#git-instructions">Git instructions</a>

## 2. Required Packages

### Linux

```bash
sudo apt-get -y update
sudo apt-get -y install libblas-dev liblapack-dev libfftw3-dev gfortran
sudo apt-get -y install gnuplot  # Needed if you want to use the gnuplot plotting macro


curl https://developer.download.nvidia.com/hpc-sdk/ubuntu/DEB-GPG-KEY-NVIDIA-HPC-SDK | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg
echo 'deb [signed-by=/usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg] https://developer.download.nvidia.com/hpc-sdk/ubuntu/amd64 /' | sudo tee /etc/apt/sources.list.d/nvhpc.list
sudo apt-get update -y
sudo apt-get install -y nvhpc-25-7
sudo apt install nvidia-cuda-toolkit
```
gfortran allows for single core or OPEN-MP parallelization, but if you plan to run on GPU you should install the nvidia toolkit instead.
Currently only the libfftw3-dev library is linked.


## 3. Compile the `LBM` code

Navigate to the `src` folder:

```bash
cd LBM/src
```

then compile and install the executable in the target directory, defaulting to
`$HOME/bin`:

Single core compilation
```bash
make
```

OPEN-MP compilation
```bash
make MP=1
```

CUDA GPU compilation
```bash
make CUDA=1
```

For double precision code (-r8) add the flag

CUDA GPU compilation
```bash
make CUDA=1 DP=1
```

To compile from scratch add a -B flag. (Necessary if you change in between parallelization settings like CUDA or OPEN-MP or none of these).





## 6. Run the code

### Linux


Start by defining the required dimensions in the  mod_dimensions.F90 file, and compile.

Create a separate catalog preferably on a large scratch disk or work area, and copy the example/infile.in to this catalog.

The example/run.sh script may be required for large open-mp simulations, as it sets ulimit -s unlimited and defines the number of cores.

The example/uvel.orig file defines an atmospheric boundary layer if it is found in the run direcotry.

To execute the code run:

```bash
boltzmann
```


# Plotting

## Tecplot
If you have tecplot (tec360) there are `.lay` and `.mcr` files in the `run2` directory.


---

# Setting up an experiment

The following explains how to set up and configure your simulations.

Compile the code after setting the the number of countries you will include in (nc=1 in mod_dimensions.F90).

Make a directory where you will run the code.

Initially you only need the file

```bash
run2/infile.in
```

---

# Code standards

If you plan to change the code note the following:

I always define subroutines in new modules:

```Fortran90
module m_name_of_routine
! define global variables here
contains
subroutine name_of_sub
! define local variables here
...
end subroutine
end module
```

in the main program you write

```Fortran90
program seir
use m_name_of_routine
call  name_of_routine
end program
```

The main program then has access to all the global variables defined in the module, and
knows the header of the subroutine and the compiler checks the consistency between the call
and the subroutine definition.

make new -> updates the dependencies for the makefile
make tags -> runs ctags (useful if you use vim)

For this to work install the scripts in the ./bin in your path and install ctags

---

# Git instructions

When working with git repositories other than the ones you own, and when you expect to contribute to the code,
a good way got organize your git project is described in https://opensource.com/article/19/7/create-pull-request-github
This link is also a good read: <a href="https://dev.to/valeriavg/master-git-in-7-minutes-gai">Git tutorial</a>

This organization will allow you to make changes and suggest them to be taken into the original code through a pull request.

So, you need a github account.
Then you fork the repository to your account (make your personal copy of it) (fork button on github.com).
This you clone to your local system where you can compile and run.

```bash
git clone https://github.com/<YourUserName>/EnKF_seir
git remote add upstream https://github.com/geirev/EnKF_seir
git remote add origin git@github.com:<YourUserName>/EnKF_seir
git remote -v                   #   should list both your local and remote repository
```

To keep your local master branch up to date with the upstream code (my original repository)

```bash
git checkout master             #   unless you are not already there
git fetch upstream              #   get info about upstream repo
git merge upstream/master       #   merges upstream master with your local master
```

If you want to make changes to the code:

```bash
git checkout -b branchname      #   Makes a new branch and moves to it
```

Make your changes

```bash
git add .                       #   In the root of your repo, stage for commit
git status                      #   Tells you status
git commit                      #   Commits your changes to the local repo
```

Push to your remote origin repo

```bash
git push -u origin branchname   #   FIRST TIME to create the branch on the remote origin
git push                        #   Thereafter: push your local changes to your forked  origin repo
```

To make a pull request:

1. Commit your changes on the local branch

```bash
git add .                       #   In the root of your repo, stage for commit
git status                      #   Tells you status
git commit -m"Commit message"   #   Commits your changes to the local repo
git commit --amend              #   Add changes to prevous commit
git push --force                #   If using --amend and previous commit was pushed
```

2. Update the branch where you are working to be consistent with the upstream master

```bash
git checkout master             #   unless you are not already there
git fetch upstream              #   get info about upstream repo
git merge upstream/master       #   merges upstream master with your local master
git checkout brancname          #   back to your local branch
git rebase master               #   your branch is updated by adding your local changes to the updated master
```

3. squash commits into one (if you have many commits)

```bash
git log                      #   lists commits
git rebase -i indexofcommit  #   index of commit before your first commit
```

Change pick to squash for all commits except the latest one.
save and then make a new unified commit message.

```bash
git push --force             #   force push branch to origin
```

4. open github.com, chose your branch, make pullrequest, check that there are no conflicts

Then we are all synced.

If you manage all this you are a git guru.
Every time you need to know something just search for git "how to do something" and there are tons of examples out there.

For advanced users:
Set the sshkey so you don't have to write a passwd everytime you push to your remote repo: check settings / keys tab
Follow instructions in
https://help.github.com/en/github/using-git/changing-a-remotes-url

To make your Linux terminal show you the current branch in the prompt include the follwoing in your .bashrc

```bash
parse_git_branch() {
     git branch 2> /dev/null | sed -e '/^[^*]/d' -e 's/* \(.*\)/(\1)/'
}
export PS1="\[\033[01;32m\]\u@\h\[\033[00m\]:\[\033[1;31m\]\w\[\033[0;93m\]\$(parse_git_branch)\[\033[0;97m\]\$ "
```

---

# Definition of time

The model uses the following definitions for time:

"Startdate" e.g., 01/03-2020 denote that the model runs from 00:00am 01/03-2020, and this time is
set to t=0.0. (see the table below)

The system integreates the ensemble of models forward one day at the time, starting from 00:00 to
24:00 and thus outputs the solutions at midnight end of the day. Thus, the output files will
contain entries like

```
 model time                                                         Output date
 0.00000E+00  Corresponds to 00:00 in the morning of the start day  (29/02-2020)
 0.10000E+01  Corresponds to 24:00 in the night of the start day    (01/03-2020)
 0.20000E+01  Corresponds to 24:00 in the night of day 2            (02/03-2020)
 0.30000E+01  Corresponds to 24:00 in the night of day 3            (03/03-2020)
 0.40000E+01  Corresponds to 24:00 in the night of day 4            (04/03-2020)
 0.50000E+01  Corresponds to 24:00 in the night of day 5            (05/03-2020)
```

## Intervention times:

If an intevention time is defined as 15/03-2020 the intervention is avtive (through values of R)
from the morning on the 15/03 at 00:00am, which corresponds to model time t>14.0. This means that
R(t) switches from R1 to R2 from the morning of 15/03. Accordingly ir swiches from 1 to 2 for using
the correct Rmat(ir).

## Measurement times

A measurement for a particular day is located at midnight that day. Thus, given a measurement from
corona.in, e.g.,

```bash
 03/03-2020    7   154   2206
```

will be located at the same time as the output of day 3, i.e., the end of the day.
