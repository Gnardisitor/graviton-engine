# Graviton Engine

## About

An n-body simulation engine made to learn and apply multiple different algorithms in C++.

## Usage

Currently, all of the different settings are inside the `main.hpp` file as macro definitions.

## Compilation

### Linux

Make sure that GCC or LLVM is installed on your machine. It also requires OpenMP to be installed. To build, use:

```sh
cmake --preset linux-gcc                # Can replace with linux-clang or debug
cmake --build --preset linux-gcc
```

### Windows

Make sure LLVM is installed on your machine. There is no need to install OpenMP since it is already included with the LLVM installation.

```powershell
winget install -i LLVM.LLVM             # Install LLVM for Windows

cmake --preset windows-clang            # Can replace with debug
cmake --build --preset windows-clang
```

MSVC is not supported due to it having a much older OpenMP implementation that struggles with SIMD instructions and vectorization, therefore making it way too slow. On a test version of the code, the MSVC compiled binary was nearly three times slower.

## Performance

On a Ryzen 7 9800X3D on Linux using 16 threads, a 5000 body cluster simulation using position verlet takes around 78 seconds to complete. Performance can greatly vary across systems depending on processors and systems. On a Ryzen 7 7840HS laptop, it performed better with 8 threads instead of 16, most likely due to it being unable to keep up with two threads performing operations per core. This does not seem to be an issue with the 9800X3D most likely due to its higher power envelope since it is a desktop processor.

## Future

The remaining empty classes for different initial conditions, physics engine, and ODE solvers will slowly be implemented. The goal in the end is for this project to become a library that can be used for different projects, but currently is still used as an exectuable for basic testing.
