<<<<<<< HEAD
# phasta01
Modified PHASTA solver to handle mixture of non-reacting thermally perfect gases in thermal non-equilibrium.

Cloned from PHASTA and modified various files to expand capability. 
=======
#build and test

    wget www.scorec.rpi.edu/~cwsmith/phastaChefTests.tar.gz .
    tar xzf phastaChefTests.tar.gz # use for CASES path below
    
Note, the following disables the SVLS and PETSC solvers and relies on LESLIB for the incompressible solver and the native compressible solver.

    cmake \
    -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_CXX_COMPILER=g++ \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DCMAKE_BUILD_TYPE=Debug \
    -DPHASTA_INCOMPRESSIBLE=ON \
    -DPHASTA_COMPRESSIBLE=ON \
    -DPHASTA_USE_LESLIB=ON \
    -DLESLIB=/path/to/libles.a \
    -DPHASTA_USE_SVLS=OFF \
    -DPHASTA_USE_PETSC=OFF \    
    -DPHASTA_TESTING=ON \
    -DCASES=/path/to/phastaCases/ \
    ..

    make

    ctest
>>>>>>> 61040f0701708be686128d6eaf805539d6d7981e
