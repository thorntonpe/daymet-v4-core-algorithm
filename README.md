# daymet-v4-core-algorithm
The purpose of this repository is to provide a definitive reference for the core algorithms used to produce the Daymet V4 dataset, as described in a manuscript currently undergoing review for Nature Scientific Data.
The core algorithms are defined by the code in interpolate.d, predict.d, and xval.d subdirectories.
Include files needed for the interpolation, prediction, and cross-validation routines are in the include.d subdirectory.
Some basic utilities that provide important context for the core algorithms are defined by code in libfio.d, libgeo.d, and tools.d.
The bin and library.d subdirectories are empty and will be populated if the provided make system is executed.
NOTE: This repository is not intended to support an operational or production environment for generating the Daymet dataset or for applying the Daymet core algorithms in new locations or with new station inputs. A full operational environment requires many additional utilities for input data processing and metadata support. The repository's sole purpose is to provide the definitive expression of the Daymet V4 algorithms. The Daymet V4 dataset was generated using this exact collection of code, with the addition of the utilities needed to operationalize the processing of the data workflow.
