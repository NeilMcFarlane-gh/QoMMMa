
![logo](logo.png)

QoMMMa is an additive QM/MM package maintained within the Harvey group at KU Leuven. It uses Tinker for any MM-side calculations, and QM-side calculations can be performed with a variety of quantum chemistry codes including: Gaussian, ORCA, xTB, Jaguar, and Molpro. The source code is written in Fortran90, and Python3 is used as a scripting language to wrap the source code and generate job files for the QM and MM packages.

QoMMMa can perform QM/MM optimisations, reaction pathway optimisations (with nudged elastic band, adiabatic mapping, and growing string method), and frequency analyses. As a code, it is continuously a work-in-progress, but the techniques described above *should* have functionality.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Examples](#examples)
- [Citations](#citations)

### Prerequisites

To install QoMMMa, you need a Fortran90 compiler such as gfortran or ifort (please note that the makefiles in the installation procedure are in gfortran format, so if you wish to use another compiler, these should be changed appropriately).
To run QoMMMa, you need any version of Python3.
To run QM calculations, you will need *at least* one of the options for the QM calculation, which are:
* xTB
* Gaussian
* ORCA
* Molpro
* Jaguar

## Installation

To install QoMMMa, starting in the QoMMMa directory, use the following bash commands:

```bash
cd src
make
cd freq
make
cd ../tinker
./qommma_compile.make
./qommma_library.make
./qommma_link.make
cp analyze_qommma.x ../../bin/analyze_grad
cp minimize_qommma.x ../../bin/minimize
cd ../../
```

You will then need to change the QOMMMA variable in the file qommma, to reflect the installed directory directory.

## Usage

Instructions on how to use the software, including code examples and explanations.

```python
# Example code
```

## Examples

Provide examples of how to use the software, including input data and expected output.

```python
# Example code with input and output
```

## Citations

If you use this software in your research, please cite it as follows:

```
@article{YourCitation,
  author    = {Author Name and Author Name},
  title     = {Project Title},
  journal   = {Journal Name},
  year      = {202X},
  volume    = {X},
  number    = {X},
  pages     = {X--X},
  doi       = {XX.XXXX/XXXXX},
}
```

```
BibTeX or other citation formats as applicable
```
