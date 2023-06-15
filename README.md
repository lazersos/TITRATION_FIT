# TITRATION_FIT
Code to fit a H2S04 titration curve using archaic methods.

# Requirements
The following python packages are required to run this program: `numpy`, `scipy`, `matplotlib`
The setup script (if it works on your machine) should install these automatically

`python3 setup.py install`

or

`python3 setup.py install --user`

# Running the code
The code takes a two column CSV file with the first column being your ordinate in mL and the second column being the abscissa in mV.

`./titration_fit.py --file myfile.csv`

or 

`./titration_fit.py --help`

