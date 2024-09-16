# Changelog

## [Version 0.0.1] - 2024-09-15 @caterer-z-t
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