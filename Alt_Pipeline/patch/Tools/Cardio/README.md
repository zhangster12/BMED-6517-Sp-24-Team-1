# Cardio
MATLAB package for processing cardiovascular signals. Admin: [Jonathan Zia](mailto:zia@gatech.edu).

## Table of Contents

1. [Summary](#summary)
2. [How to Use](#how-to-use)
3. [Policies](#policies)

## Summary

This is a MATLAB package for processing cardiovascular signals of the following modalities: EEG, ICG, SCG, and PPG. This README file contains general information for using the package; see each sub-folder for a sub-package-specific README file.

## How to Use

Note that library and data files should be added to the MATLAB path before using his library. This script may be modified and either run before use or added to your `startup.m` file. Using `genpath` will also add custom classes located in the `Classes` folder to your path.
```matlab
addpath(genpath('/path/to/Cardio/'))
addpath(genpath('/path/to/Data/'))
```

Functions in MATLAB packages may be called using the synatx `<package>.<sub_package>.<function>`. For example, the function located in `Cardio/+cardio/+general/ema.m` may be called via `cardio.general.ema()`. For any function in this package, a description of the function and its input arguments may be obtained with the command `help cardio.<sub_package>.<function>`. If the function is a member of a class, help may be obtained with the command `help <class>.<function>`. 

Member functions are organized into sub-libraries based on the primary modality that they process. Dataset-specific processing scripts are contained in the `+data` directory. The `Classes` folder contains custom classes in sub-folders (denoted by the syntax `@ClassName`) and enumerations. Enumerations are data structures which may take values from a specified list. This is useful for avoiding syntax errors when setting flags throughout the package. 

For example, consider a function which may compute the distance between two signals in one of two ways. Having the user type the flag `euclidian` or `manhattan` may result in syntax errors that go unnoticed. Instead, we define the enumeration `Distance`, which can take the values `Distance.euclidian` or `Distance.manhattan`. This auto-fills for the user in the MATLAB IDE, preventing syntax errors. Enumerations are used throughout the package, and may be called via `<name>.value`. A description for enumerations used in this package are given below.

## Policies

Modifications to the Cardio package should not be pushed directly to the master branch. To make modifications, create a branch from the most recent commit in the master branch. When the modifications are complete, a merge request should be made, which will be reviewed and approved by the admin.

Code in the master branch should be **clean**, **modular**, **portable**, and **thoroughly-commented**. 

**Clean** code (1) is free of extraneous scripts, (2) is clearly-divided into sections, (3) uses descriptive variable names, and (4) is concise; e.g. avoids multi-level nests.

**Modular** code divides repeated tasks into sub-functions -- copy/paste should not be used, when possible.

**Portable** code should function seamlessly with the rest of the package without throwing errors.

**Thoroughly-commented** code should be clearly-understood with a quick pass over the comments. The proper heuristic for this package is one line of comments per one line of code.

Additionally:
1. Functions should make use of `varargin` and `varargout` when possible and should not have long lists of input or output arguments. 
2. Functions must include a block comment at the top of the function which may be queried with the `help` command. These blocks should include a summary of the function complete description of all input arguments.
3. When possible, enumerations should be used rather than text-based flags, as aforementioned.  
4. For groups of scripts that are intended for a specific use -- such as processing a specific dataset -- object-oriented programming is preferred and should be employed. 
5. When possible, functions should be **generalizable** so that they may be used for a multitude of projects or purposes.

In general, please leave this package better than you found it.
