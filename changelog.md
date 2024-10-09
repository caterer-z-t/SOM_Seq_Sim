# Changelog

## [Version 0.0.0] - 2024-09-16 @caterer-z-t
- Added github/workflows to automatically test the code as it is being pushed and commited to the repo.
- Added code to `utils.py` to create genetic sequencing data and plot the new sequencing data compared to the original data.
- Updated `dev_sequencing_code.py` to take positional arguments and allow for help argumengts for usage. Modified `README.md` to show usage and implementation of this.
- Converted `utils.py` to be a class rather than multiple functions, subsequently changed `dev_sequencing_code.py` to implement this change
- added `__init__.py` files for easier referencing files and functions
- added comments and files and directories to .gitignore to be left out from uploads to github
- updated `README.md` file to contain proper structure of package and software thus far
- modified `utils.py` to take the arguments of random_row which will highlight the row in the distributions and proportions figures, if not specified (FALSE) the traditional distributions are still generates
- updated `dev_sequencing_code.py` to take the argument from the command line for the file path of the data used, and updated the use of the functions in `utils.py` to correctly call each function with updated parameters and usage
- created `tests/seq_tests.py` to test the base functions of `utils.py` and initiated the testing scripts for this project.

## [Version 0.0.0] - 2024-09-15 @caterer-z-t
- Added dependencies to the `env.yml` file
- Added directory of data to `.gitignore` file
- Included how to obtain example data in `README.md` file
- Wrote a few examples of generated unimodal, bimodal, and trimodal distributions in the files of `unimodal_example.py`, `bimodal_example.py` and `trimiodal_example.py`. 
- Wrote a few functions in `utils.py` to calculate optimal numbers of rows and columns to be used in a figure. Additionally wrote functions to plot and show the distributions of numerical and categorical data from a file. (use case as of now is metadata from experiment). 
- Running the functions from `utils.py` in the python script of `dev_sequencing_code.py`, using the file path from the first argument in the bash script. To run this script please use the code
  ``` bash
  python3 {filepath to `dev_sequencing_code.py`} {filepath to metadata.tsv}
  ```

## [Version 0.0.0] - 2024-09-12 @caterer-z-t
- Created a `Chaneglog.md` file and moved the changelog informaiton from the `README.md` file to this file. 
- Added readthedocs, used base template from [wrmXpress-gui](https://github.com/wheelerlab-uwec/wrmXpress-gui) repo
- Initializing a `README.md ` template to be updated and solidified later.
  - Updated `README.md` file to have proper md formatting. 
- Created a template `env.yml` file for the virtual environment for this project.