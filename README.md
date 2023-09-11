# Article Pipeline Data Validation Scripts (ap_validation)
This repository contains the scripts needed to ensure the quality of the zipped NMR data deposited into the Article Pipeline.

#### WARNING: The main branch of this repo is a python package with a setup.py file to allow for pip installation. Changes to this branch could effect production versions of applications that depend on it.

# Installation
To install the apvalidation package into your own Python project, just use pip along with the git link to this project.
```
pip install @git+https://github.com/liningtonlab/apvalidation.git
```

## (REQUIRED FOR DEVELOPMENT) Submodule Installation
If you are planning to clone this repo to your computer it is important to note that it relies on external files in order to perform standardization steps. This package was added to this repo like such...

```
git submodule add https://github.com/np-mrd/npmrd_data_exchange.git npmrd_data_exchange
```

In order to update this submodule whenever these files change remember to run...

```
git submodule update --remote
```

# Functionality
The code in this repository has four main functionalities.
1. Traverse the deposited NMR zip folder to determine the type of machine used in the experiments (Bruker, Varian, or JEOL) and the filetype of the data.
2. Read the parameter files to extract information about the experiment.
3. Standardize and display the structures of the compounds using the smiles string submitted by the user.
4. Validate the format of the submitted peak lists.

The first two sections of this documentation will focus on the extraction of parameters from the data files.

# Extracting Parameters from the Data
This function uses functionality one and two from the four functionalities above.

The high-level steps taken to extract this data is outlined in the diagram below.


<img src="https://user-images.githubusercontent.com/55040326/183517545-1a7ce3ea-137b-4238-8488-650d1dfc5d67.png" />

Each step in the above diagram will be covered in more detail below.

## 1. Determine the NMR Vendor & Filetype.

simple_file_finder.py 

MetaFinder class has three main components; finding a vendor for each experiment, finding a parameter file path, and searching for errors if exist. 
find_meta function returns a dictionary with vendor and parameter file paths. 

#### Example Return
```
MetaFinder.find_meta("/Downloads/NMR/Aspochalasin I.zip")
{ 'vendor_name': ['bruker','bruker','bruker'], 'meta_file': [['/nrm/1/acqu',/nrm/1/acqu2'],['/nrm/2'],['/nrm/3']]}
```

validator function append found error to error_message. If a new error is found, update private functions to handle the error

## 2. Handling Different Filetypes.
This section exists to deal with the need to support the JDX filetype.

### What is a JDX file?
A JDX file is a nested xml-like filetype which uses tags to nest it's data. If you want to learn more about the JDX filetype you can read about it [here](http://www.jcamp-dx.org/).

### Where do JDX files come from?
For the purpose of the deposition site we will be receiving JDX files from depositors who have either .mnova, or .jdf files. These are filetypes which CANNOT be taken by the deposition site but which CAN be exported as a JDX files to be processed that way instead.

These files are to be exported to JDX using the [MestreNova Software](https://mestrelab.com/download/mnova/).

### How are JDX files handled?
JDX files are not immediately compatible with the python package we use to read parameter files [nmrglue](https://www.nmrglue.com/). 
The two reasons these files do not work with nmrglue are:

1. JDX files have extra padding tags on the top and the bottom of the file.
2. JDX files can have information about more than one experiment in the same file. These experiments are separated via tags in the file.

Using the two reasons above the [mnova_jdx_reader.py](https://github.com/liningtonlab/apvalidation/blob/main/apvalidation/mnova_jdx_reader.py) file strips these extra tags and separates any muli-experiment JDX files into their own individual files.

## 3. Reading the Parameter Files.
Once the extraction is at this step, all files should be in an acceptable format for [nmrglue](https://www.nmrglue.com/) to process.

The final step is to now extract the parameters from the data file to determine information about the experiment. This is done through the Varian, Bruker, and Jcampdx classes in [extract.py](https://github.com/liningtonlab/apvalidation/blob/main/apvalidation/extract.py). The parameters to be extracted are:

- experiment type
- F2 nucleus
- F1 nucleus
- frequency
- temperature (kelvin)
- solvent


### The Extraction Code
The following code can all be found in [extract.py](https://github.com/liningtonlab/apvalidation/blob/main/apvalidation/extract.py).

There are three main classes defined in extract.py. These are Varian, Bruker, and Jcampdx. Each of these classes are specialized to handle different types of data files. 

- procpar --> Varian
- acqu --> Bruker
- jdx --> Jcampdx

Each of these classes contain methods which are able to read and parse the corresponding file. Most of the methods in the classes are helper functions which are not meant for user use. As a user there are two main methods that you should be concerned with; find_params() and read(). 
- read(): This method does exactly what is sounds like, it reads the raw parameter file. The input to this method should be a string representation of the filepath to the desired parameter file. The output is a python dictionary object which stores the contents of the file in key, value pairs.
- find_params(): This method does the bulk of the heavy lifting. It takes the parameter dictionary returned by the read() method and extracts all the needed parameters. A lot of the work is done by calling additional helper functions internally, but what matters is that the output is a python dictionary containing just the parameters we want in a standardized format.

#### Example Usage (Varian Machine)
```
  import extract
  
  filepath = "./procpar"
  parameter_dict = extract.Varian.read(filepath)
  desired_parameters = extract.Varian.find_params(parameter_dict)
```
# Validating SMILES Strings
When a SMILES string is submitted to the deposition site, it must be checked for validity. 
Validating SMILES and producing images of those structures is the third functionality performed by apvalidation.

The code to perform this function can be found in the [smiles_validation.py](https://github.com/liningtonlab/apvalidation/blob/main/apvalidation/smiles_validation.py) file. 

The function validate_struct(), contained in this file does both the validation, and the structure image creation at the same time. The inputs to this function are:

validate_struct(smiles, img_path, asInchiKey=False):

- str smiles: The SMILES submitted by the depositor
- str img_path: The path to the where the structure image should be saved after it is created.
- bool asInchiKey: When True, the function saves the image with the InchiKey as the file name. The default uses the SMILES as the filename.

Example Usage:

```
validate_struct(smiles="CC=C", img_path="./image_folder")

# This ex. checks that CC=C is a valid SMILES and then saves the structure as
# an image in the ./image_folder using SMILES as the image filename.
```

# Validation Peak Lists

The final functionality implemented by apvalidation is the peak list validation. When a peak list is submitted to the deposition site, that list must be checked for the proper format.

The [peak_validator.py](https://github.com/liningtonlab/apvalidation/blob/main/apvalidation/peak_validator.py) contains all the logic to validate these peak lists.

Below you can see a diagram which leads through the different checks performed on a peak list before it is confirmed to be valid.


<img src="https://user-images.githubusercontent.com/55040326/183711596-f9ff52e4-794e-4766-832e-f12db03967ad.png"/>



## Unzip files.
extract_core.py

extract_core_file function finds a necessary file in the submitted zip file and extracts them
 
```
necessary file list
varian = ["fid", "procpar", "log", "text", "log 2", "procpar 2", "text 2"]
bruker = ["fid", "ser", "acqu", "acqu2", "acqus", "acqu2s"]
jcamp = ["jdx"]
```

Directory name convention : {f1_nucleus}_{f2_nucleus}_{experiment_type}
