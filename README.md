# Knitogami

To build, see the original [wiki](https://github.com/wimvanrees/growth_SM2018/wiki) for instructions on how to compile and run.

# Workflow

The build process outputs an executable, `shell`, in the `/bin` directory.

Create a folder within `/bin` to store all files related to a particular run, e.g.:
```
cd bin
mkdir run_squiggles
cd run_squiggles
```

Copy the executable and any relevant knit pattern files to the run directory:
```
cp ../shell .
cp ../patterns/knit_squiggles.txt .
```

Run the executable with any specified parameters:
```
./shell -sim knit -filename knit_squiggles.txt -xCurv -1.0 -yCurv 1.0 \
    -res 0.1 -dx 1.4 -dy 1.2 -h 0.5 -E 1.0 -nu 0.4 \
    -boundarySpringConstant 0.0 -nSteps 10
```

The result is a set of VTP files containing the resulting meshes at each quasistatic simulation step, as values of $\kappa_x$ are increased from $0$ to $xCurv$ and $\kappa_y$ is increased from $0$ to $yCurv$. VTP files can be viewed with [Paraview](https://www.paraview.org).
