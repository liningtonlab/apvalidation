# Article Pipeline Data Validation Scripts (ap_validation)
This repository contains the scripts needed to ensure the quality of the zipped NMR data deposited into the Article Pipeline.

#### WARNING: The main branch of this repo is a python package with a setup.py file to allow for pip installation. Changes to this branch could effect production versions of applications that depend on it.


## Functionality
The code in this repository has four main functionalities.
1. Traverse the deposited NMR zip folder to determine the type of machine used in the experiments (Bruker, Varian, or JEOL) and the filetype of the data.
2. Read the parameter files to extract information about the experiment.
3. Standardize and display the structures of the compounds using the smiles string submitted by the user.
4. Validate the format of the submitted peak lists.

The first two sections of this documentation will focus on the extraction of parameters from the data files.

# Extracting Parameters from the Data (Sections 1&2)
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
For the purpose of the deposition site we will be recieving JDX files from depositors who have either .mnova, or .jdf files. These are filetypes which CANNOT be taken by the depostion site but which CAN be exported as a JDX files to be processed that way instead.

These files are to be exported to JDX using the [MestreNova Software](https://mestrelab.com/download/mnova/).

### How are JDX files handled?
JDX files are not immediately compatible with the python package we use to read parameter files [nmrglue](https://www.nmrglue.com/). 
The two reasons these files do not work with nmrglue are:

1. JDX files have extra padding tags on the top and the bottom of the file.
2. JDX files can have information about more than one experiment in the same file. These experiments are separated via tags in the file.

Using the two reasons above the [mnova_jdx_reader.py](https://github.com/liningtonlab/apvalidation/blob/main/apvalidation/mnova_jdx_reader.py) file strips these extra tags and separates any muli-experiment JDX files into their own individual files.

## 3. Reading the Parameter Files.
Once the extraction is at this step, all files should be in an acceptable format for [nmrglue](https://www.nmrglue.com/) to process.

The final step is to now extract the parameters from the data file to determine information about the experiment. This is done through the Varian, Bruker, and Jcampdx_Handler classes in [extract.py](https://github.com/liningtonlab/apvalidation/blob/main/apvalidation/extract.py). The parameters to be extracted are:

- experiment type
- F2 nucleus
- F1 nucleus
- frequency
- temperature (kelvin)
- solvent


### The Extraction Code
The following code can all be found in [extract.py](https://github.com/liningtonlab/apvalidation/blob/main/apvalidation/extract.py).

There are three main classes defined in extract.py. These are Varian, Bruker, and Jcampdx_Handler. Each of these classes are specialized to handle different types of data files. 

- procpar --> Varian
- acqu --> Bruker
- jdx --> Jcampdx_Handler

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

#### extract.py Execution Flowchart
The following chart shows the execution of the parameter extraction scripts from calling the read function to how all the helper functions contribute to find_params().
To get a full view of the flowchart please download the image and view it locally.
![text](https://github.com/liningtonlab/ap_validation/blob/main/NP%20Validation%20extract.py.png)

## 3. Unzip files.
extract_core.py

extract_core_file function finds a necessary file in the submitted zip file and extracts them
 
```
necessary file list
varian = ["fid", "procpar", "log", "text", "log 2", "procpar 2", "text 2"]
bruker = ["fid", "ser", "acqu", "acqu2", "acqus", "acqu2s"]
jcamp = ["jdx"]
```

Directory name convention : {f1_nucleus}_{f2_nucleus}_{experiment_type}
