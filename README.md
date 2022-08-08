# Article Pipeline Data Validation Scripts (ap_validation)
This repository contains the scripts needed to ensure the quality of the zipped NMR data deposited into the Article Pipeline.

#### WARNING: The main branch of this repo is a python package with a setup.py file to allow for pip installation. Changes to this branch could effect production versions of applications that depend on it.

## Dependencies
First and most importantly, these scripts will require an installation of the python package rdkit. This is important to note since rdkit's installation process is particularily complicated if you are not using an Anaconda environment.
## Functionality
The code in this repository has three main functionalities.
1. Traverse the deposited NMR zip folder to determine the type of machine used in the experiments (Bruker, Varian, or JEOL) and retrieve the parameter files.
2. Read the retrieved parameter files to extract information about the experiment.
3. Standardize and display the structures of the compounds using the smiles string submitted by the user.

## 1. Traversing the ZIP Files.

simple_file_finder.py 

MetaFinder class has three main components; finding a vendor for each experiment, finding a parameter file path, and searching for errors if exist. 
find_meta function returns a dictionary with vendor and parameter file paths. 

#### Example Return
```
MetaFinder.find_meta("/Downloads/NMR/Aspochalasin I.zip")
{ 'vendor_name': ['bruker','bruker','bruker'], 'meta_file': [['/nrm/1/acqu',/nrm/1/acqu2'],['/nrm/2'],['/nrm/3']]}
```

validator function append found error to error_message. If a new error is found, update private functions to handle the error

## 2. Reading the Parameter Files.

Once the parameter file is retrieved from the code in section (1.), some information from this file must be extracted. The experiment parameters that are needed from this file are the following.
- experiment type
- F2 nucleus
- F1 nucleus
- frequency
- temperature (kelvin)
- solvent

### User Methods
The functions for retrieving these parameters can be found in the "extract.py" file under the "paramExtract/packages" folder. In this file you will notice three different classes; Varian, Bruker, and Jcampdx. These three classes reference to the three different types of machines that produce the NMR data. Since the machines produce different parameter files  we need different functions to deal with each type.

Each of these classes contain methods which are able to read and parse the parameter file produced by the correpsonding NMR machine. Most of the methods in the classes are helper functions which are not meant for user use. As a user there are two main methods that you should be concerned with; find_params() and read(). 
- read(): This method does exactly what is sounds like, it reads the raw parameter file produced by the corresponding NMR machine. The input to this method should be a string representation of the filepath to the desired parameter file. The output is a python dictionary object which stores the different parameters in key, value pairs.
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
