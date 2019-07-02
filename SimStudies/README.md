# How to run the scripts

All simulation studies are to be run from their respective folders, 
i.e., study A is to run from SimA_IM as its working dir.

All figures from the paper will be deposited in a figures folder.

Before you run a simulation study, make sure you have installed the 
packages used within it. This is most easily done by installing the 
[`tsars` package](https://github.com/nielsrhansen/SLODE) first. 
This auxilary package is also used for the simulations.

Some of the simulations are computationally heavy, in particular Study B
and Study E. For the paper, they were run on a 128 core server using multiple
cores. 

The following gives the commands for running the scripts as R batch processes
on a Linux server supposing that you start from the SimStudies directory.

## Sim A:

```
cd SimA_IM
R CMD BATCH test_IM.R &
cd ..
```

## Sim B:

```
cd SimB_MAK
R CMD BATCH main.R &
R CMD BATCH R/read_results.R &
cd ..
```

## Sim D

```
cd SimD_GRADE
R CMD BATCH main.R &
cd ..
```

## Sim E

```
cd SimE_Glycolysis
R CMD BATCH main.R &
R CMD BATCH read_results.R &
cd ..
```

## Sim F

```
cd SimF_EnvZOmpR
R CMD BATCH main.R &
cd ..
```









