# Lattice-Boltzmann-Model
This repository contains a 3D implementation of a Lattice-Boltzmann model on a D3Q19 or D3Q27 lattice for high Reynolds number flow.

## License

This project is dual-licensed:

- Open Source: [Apache License 2.0](./LICENSE-APACHE) for research and academic use.
- Commercial: A commercial license is available for proprietary use.

If you are interested in a commercial license, please contact:
**[Geir Evensen] – [geir.evensen@gmail.com]**

## Introduction

The code, NVIDIA's CUDA Fortran, runs on a single core, an OPEN-MP multicore, or a GPU, with optimization primarily for the GPU.

The code runs in single precision as default but a flag copiles a double precision version.

The collision operator is a single relaxation time third-order expansion with regularization similar to [Jacob et al (2018)](https://hal.science/hal-02114308) and
[Feng et al, (2018)](https://doi.org/10.1029/2020MS002107), but excluding the hybrid regularization, which I have not yet needed. It is also possible to use it for
viscous flow without regularization and turbulence closure, and using a second-order BGK expansion.

The turbulence closure scheme is the one described by [Vreman (2004)](https://doi.org/10.1063/1.1785131).

The model boundary conditions are periodic or inflow-outflow in the i-direction, periodic or closed no-slip or free-slip two-timestep bounceback in the j- and k-
directions.

The code allows for inserting solid bodies within the model domain to simulate, e.g., flow around an airfoil or a cylinder.

Additionally, there is a complete implementation of an actuator line model for the NREL-5Mw wind turbine, and it is possible to include multiple
turbines at any location of the model domain.

Inflow turbulence is mimicked or introduced at a section inside the inflow boundary at i=1 (typically at the slice i=10) by applying a smooth in space
and time, pseudo-random force on the fluid.

The forcing function for the inflow turbulence and the turbines is that of Kupershtokh (2009). For the turbine forcing, it is also possible to run with the forcing
formulations of Guo et al (2002). However, when using regularization, we project the non-equilibrium distribution onto the third-order Hermite polynomials. The
significant difference between the Guo scheme's equilibrium distribution computed on the forcing-updated velocities leads to a poor representation using Hermite
polynomials, and we partly lose the effect of the forcing. Thus, in the Guo scheme with regularization, it is necessary to compute the regularization first on
R(fneq)=R(f-feq(u)) and then recover f=feq(u)+R(fneq) before computing the updated forcing velocities u+du and the forcing distribution df. Next we compute
feq(u+du) and fneq(u+du)= f-feq(u+du), which goes into the collision and vreman calls. Thus, the cost of Guo is therefore much higher as it requires two calls to
fequil and an extra computation of f=feq +R(fneq) and extra computation of fneq= f-feq(u+du). It is possible to reduce the computational cost by updating only at
the turbine locations, but for now, it is not worth the effort.



<p align="center">
  <a href="#installation">Installation</a> *
  <a href="#code-standards">Code standards</a> *
  <a href="#setting-up-an-experiment">Experiment setup</a> *
  <a href="#plotting">Plotting</a> *
  <a href="#git-instructions">Git instructions</a> *
  <a href="https://github.com/geirev/LBM/blob/main/LICENSE">License</a>
</p>

<p align="center">
<img src="example/windturbine.png" width="300">
</p>

<p align="center">
<img src="example/city.png" width="300">
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
git remote add upstream git://github.com:geirev/LBM
popd
```
or, if you have not set up git-ssh
```bash
git clone git@github.com:<userid>/LBM.git
pushd LBM
git remote add upstream https://github.com/geirev/LBM
```

If you are new to Git, read the section <a href="#git-instructions">Git instructions</a>

## 2. Required Packages

### Linux

```bash
sudo apt-get -y update
sudo apt-get -y install libfftw3-dev               # fft library used when sampling pseudo-random fields

# nvidia Cuda fortran compiler and utilities installation (get the latest version, nvhpc-25-7 is just for illustration).
curl https://developer.download.nvidia.com/hpc-sdk/ubuntu/DEB-GPG-KEY-NVIDIA-HPC-SDK | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg
echo 'deb [signed-by=/usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg] https://developer.download.nvidia.com/hpc-sdk/ubuntu/amd64 /' | sudo tee /etc/apt/sources.list.d/nvhpc.list
sudo apt-get update -y
sudo apt-get install -y nvhpc-25-7
sudo apt install nvidia-cuda-toolkit
```

## 3. Compile the `LBM` code

For gpu compilation run 'nvidia-smi' or 'lshw -C display' to find you gpu-card.
The check the compute capability of your gpu in the table https://developer.nvidia.com/cuda-gpus.
Navigate to the `src` folder and open the makefile and specify the correct -gpu=ccXX flag.


```bash
cd LBM/src
```

then compile and install the executable in the target directory, defaulting to
`$HOME/bin`:

Default is compilation for single core in single precision on D3Q27 lattice.
```bash
make -B
```

Default is compilation for single core in single precision on D3Q27 lattice using gfortran.
```bash
make -B GFORTRAN=1
```

Default is compilation for OPEN-MP in single precision on D3Q27 lattice using gfortran.
```bash
make -B GFORTRAN=1 MP=1
```

Compilation for single core in single precision on D3Q19 lattice
```bash
make -B D3Q19=1
```

Compilation for OPEN-MP in single precision on D3Q27 lattice
```bash
make -B MP=1
```

Compilation for CUDA GPU in single precision and D3Q27 lattice
```bash
make -B CUDA=1
```

Compilation for CUDA GPU in double precision and D3Q27 lattice
```bash
make -B DP=1 CUDA=1
```

To recompile from scratch add a -B flag. (Necessary if you change in between parallelization settings like CUDA or OPEN-MP).

E.g., recomile for GPU in double precision using D3Q19 lattice
```bash
make -B CUDA=1 DP=1 D3Q19=1
```

Running in single precision is about twice as fast as using double precision.

Running on the D3Q19 lattice reduces the CPU time with around 40 % but seems to introduce more noise.

A test running 200 time steps on single CPU, with OPEN-MP, and GPU for a domain of 121x121x928 gave the following wall times:

```bash
single-core     (gfortran)  : 915.84 s   (make -B GFORTRAN=1)
single-core     (nvfortran) : 641.48 s   (make -B)

open-mp 6 cores (gfortran)  : 295.57 s   (make -B GFORTRAN=1 MP=1)

open-mp 6 cores (nvfortran) : 244.28 s   (make -B MP=1)
open-mp 13 cores (nvfortran): 204.78 s   (make -B MP=1)
open-mp 16 cores (nvfortran): 202.78 s   (make -B MP=1)
open-mp 24 cores (nvfortran): 253.74 s   (make -B MP=1)

GPU             (nvfortran) :  36.12 s   (make -B CUDA=1)
```
The simulations were run on a "Lenovo Legion 7 Pro" laptop with a "Core Ultra 9 275 HX" (having 24 indepenent cores) and the gpu card is "Nvidia RTX 5090."


<p align="center">
<img src="doc/gpu.png" width="400">
</p>
These plots clearly show that GPU is the optimal choice for heavy simulations, while CPU scaling beyond 10-16 cores is inefficient.

For the GPU simultion the timing of different routines were as follows after a first optimization pass doing the standards:
---
```bash
---------------------------------------------------------------------------------------------------
                           gfortran          nvfortran     nvfortran      nvidia            Speedup
                             1 core             1 core      16 cores         GPU     nvf 1 core/GPU
---------------------------------------------------------------------------------------------------
initialization     time =      0.46               0.69          0.70        0.70
turbine forcing    time =      9.37               5.35          2.42        0.60
turbulence forcing time =      0.84               0.49          0.17        0.06
equil              time =    344.73             186.11         59.24        8.35              41.28
regularization     time =    377.08             252.33         77.45       11.71              21.55
vreman             time =     87.00              97.88         20.06        1.19              82.25
collisions         time =     24.85              24.98         10.64        4.21               5.93
applyturbines      time =      0.49               0.39          0.60        0.06
applyturbulence    time =      0.07               0.06          0.08        0.01
solids             time =      0.00               0.00          0.00        0.00
boundarycond       time =      1.68               1.73          1.19        0.89
drift              time =     46.70              45.84         21.85        4.49              10.20
macrovars          time =     20.14              22.71          5.33        0.74              30.68
diag               time =      1.34               1.89          2.00        1.98
averaging and turb time =      0.19               0.22          0.21        0.20
final stuff        time =      0.82               0.75          0.76        0.86
Total wall time    time =    915.84             641.48        202.77       36.12              17.75
---------------------------------------------------------------------------------------------------
```





## 4. Run the code

Start by defining the required grid dimensions in the src/mod_dimensions.F90 file, and compile.

Create a separate catalog preferably on a large scratch disk or work area, and copy the example/infile.in to this catalog.
The example infile.in corresponds to a  2D city case with flow through the city, and the program should run without any other input files. Just choose
the city2 case in mod_dimensions.F90. You can also use this grid for the cylinder case by editing the infile.in. For more realistic 3D runs increase the
vertical resolution.

The example/run.sh script may be required for large grids, as it sets ulimit -s unlimited and it also defines the number of cores used in OPEN-MP simuilations.

The example/uvel.orig file defines an atmospheric boundary layer if it is found in the run direcotry.

To execute the code run:

```bash
boltzmann
```

## 5. Code profiling
To profile the code run
```bash
nvprof boltzmann
```
which gives a detailed listing of the CPU time used by each kernel.
It is clear that the third order expansion and the regularization are by far the most expensive computations in the code.


## 6. Plotting

The current code version outputs Tecplot plt files read by tec360.

The plotting routine is m_diag.F90 that dumpis a file tecGRID.plt containing the grid layout, i.e., the i, j, and k indices and the blanking variable.
For each solution time the routine saves the density and three velocity components in each grid point.
When using Tecplot one must load the tecGRID.plt file as the first file and add any number of solution files.
Thus, the diagnostic saved is minimal, and we compute absolute velocity and vorticity, and Q-criterion diagnostics within Tecplot using the loaded velocity
fields.

Additionally, it is possible to compute the averages over any number of time-steps and these are then saved to tecAVERAGE.plt. This file also contains the
turbulent kinetic energy.

If you are not using Tecplot you can chose another output format by replacing the call to the tecout routine in m_diag.F90 to your choice.

---

## 7. Code standards

If you plan to change the code note the following:

I always define subroutines in new modules:

```Fortran90
module m_name_of_subroutine
! define global variables here
contains
subroutine name_of_subroutine
! define local variables here
...
end subroutine
end module
```

in the main program you write

```Fortran90
program name
use m_name_of_subroutine
call  name_of_subroutine
end program
```

The main program then has access to all the global variables defined in the module, and
knows the header of the subroutine and the compiler checks the consistency between the call
and the subroutine definition.

make new -> updates the dependencies for the makefile
make tags -> runs ctags (useful if you use vim)

The current makefile updates the dependencies at every compilation, so if you add a file with a new subroutine you can just type make and it will be included
in the compilation.

For this to work install the scripts in the ./bin in your path and install ctags

---

## 7. Git instructions

When working with git repositories other than the ones you own, and when you expect to contribute to the code,
a good way got organize your git project is described in https://opensource.com/article/19/7/create-pull-request-github
This link is also a good read: <a href="https://dev.to/valeriavg/main-git-in-7-minutes-gai">Git tutorial</a>

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

To keep your local main branch up to date with the upstream code (my original repository)

```bash
git switch main             #   unless you are not already there
git fetch upstream              #   get info about upstream repo
git merge upstream/main       #   merges upstream main with your local main
```

If you want to make changes to the code:

```bash
git switch -c branchname      #   Makes a new branch and moves to it
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

2. Update the branch where you are working to be consistent with the upstream main

```bash
git switch main             #   unless you are not already there
git fetch upstream              #   get info about upstream repo
git merge upstream/main       #   merges upstream main with your local main
git switch brancname          #   back to your local branch
git rebase main               #   your branch is updated by adding your local changes to the updated main
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
