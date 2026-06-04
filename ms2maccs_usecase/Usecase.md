# MS2X for MS2MACCS

[MSMACCS repository](https://github.com/j-a-dietrich/MS2MACCS)

MS2MACCS predicts MACCS fingerprints for positive and negative mode from MS2 data. It relies on "substructure MACCS" (MACCS fingerprints for a specific fragment/mz). These "substructure MACCS" were approximated with MS2X.

## Data Mining for MS2MACCS

The file data_mining_ms2maccs.py is a python command line program that was used to generate a dictionary object (key is mz rounded to the second digits after the comma; value is summed MACCS fingerprints). The resulting dictionary is stored as a pickle object.

### ⚡Multi-core application

The script utilizes all cores available to be able to generate fragments for more than 20 000 spectra with on average over 100 fragments in 8 hours.

### 🗂️ Variables to adapt
*path_name*: mgf file that contains peaks and SMILES for the query spectra <br>
*save_name*: name for pkl dictionary file that will be generated <br>
*charge*: 1 for positive mode; -1 for negative mode <br>

### ✨ Adaptability
The current version also includes a subformula vector attached to the maccs fingerprint vector (10 + 167) resulting in a dimension of 177 for the vector in the stored dictionary.

### 🔍 Data used
Soon we will upload the input data used in this study

## Combining dictionaries

The file combine_fps.py import all stored dictionaries used in the study (training, validation and test data set) and combined them to one dictionary for positive mode and one dictionary for negative mode. Additionally it combines the positive with negative mode .mgf files for the training data sets.