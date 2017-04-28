# alignTrack

### Welcome to this collaborative repo for UCL g-2 Tracker team! ###

### DIR Structure ###
1. mpIIDESY/ - main directory for C++ Tracker MC code, and pede. 
2. C_mpIIDESY/ - for experimenting with toy C++ model MC and pede.
3. C_mpIIDESY/ - for experimenting with toy Fortran model MC and pede.
4. PEDEv04-03-04/ - clean pede (most recent) checkout from DESY
5. python_toy_tracker/ - John's python code for toy tracker (with correct geometry)

### INSTALLATION ###
Here are the instructions to get the code (requires c++0x compiler suport for Logger) on working on gm2ucl at Fermilab: 

1.  `git clone https://github.com/glukicov/alignTrack.git`
to get the latest code from our repository 
2. `cd alignTrack/mpIIDESY`
3. `make`
to build the pede executable 
4. test that it works by `pede -t`
(should give a terminal output [last 2 lines]):
 Millepede II-P ending   ... Mon Dec 12 12:31:15 2016 
 Peak dynamic memory allocation:    0.100512 GB
 
To run PEDE algorithm in general case:
` ./pede str.txt  ` [where e.g. str.txt is a steering file, which specifies both - a data.bin file and a constraint file (e.g. con.txt)]

#### Random (Integer) Number Generation ####
   * `python randomIntGenerator.py -u True -o uniform_ran.txt -s <random_seed> -p <randoms_decimal_places> -n <number_of_randoms_to_generate>`
   * `python randomIntGenerator.py -g True -o gaussian_ran.txt -s <random_seed> -p <randoms_decimal_places> -n <number_of_randoms_to_generate>`
Please note that the output filenames shown are the default random number files for the modified Fortran version of *Test 1/2*, and C++ MC. Random numbers are generated using Python's built-in *Marsenne Twister* generator. Generally, ~5,000,000 random numbers in each file seems to be more than enough for the default *Test 1/2*, and C++ MC. 

 1.* `python randomIntGenerator.py -u True -o uniform_ran.txt -s 123456789 -n 5000000`
 
 2. * `python randomIntGenerator.py -g True -o gaussian_ran.txt -s 987654321 -n 5000000`
 
#### Running C++ version of mptest2.f90: ####
1. Compile code with `make -f MakeMp2test.mk`
2. Generate data by running `./Mptest2`. This generates:
   * `C_Mp2tst.bin`, `C_Mp2con.txt`, `C_Mp2str.txt`
3. Fit data by running `./pede C_Mp2str.txt`.

### Producing PEDE Histograms ### 
1. ` root -l `
2. root [0]  ` .L readPedeHists.C+`
3. root [1] ` gStyle->SetOptStat(1111111)` [to see Under/Overflows and Integrals]
4. root [2] ` readPedeHists()` [possible options inisde () "write" "nodraw" "print"] 


### To run PEDE algorithm for Fortran version of Mptest2 ###
1. ` ./pede -t=track-model`
where track-model = SL0, SLE, BP, BRLF, BRLC [see p. 154 of refman.pdf] 

e.g. ./pede -t=SL0 [check the correct parameters, aslo option for flag -ip] 


#### Running C++ Port of Test 1 ####
1. Compile code with `make -f MakeMpTest1.mk`
2. Generate random plaintext files of uniform and gaussian distributed random numbers, according to instructions above.
3. Generate data by running `./MpTest1 <uniform_randoms_file> <gaussian_randoms_file>`. This generates:
   * `mp2test1_true_params_c.txt`: A plaintext file containing true values of parameters, with their labels, for comparison with fitted values.
   * `mp2tst1_c.bin`: A binary file containing fitting data.
   * `mp2test1con_c.txt`: A plaintext file containing parameter constraints.
   * `mp2test1str_c.txt`: A plaintext file containing steering information for `pede`. (Note - it is important to use the steering file generated here, rather than that generated by `./pede -t`, as this steering file flags the data binary as being generated by C, rather than Fortran.)
4. Fit data by running `./pede mp2test1str_c.txt`. 

#### Python Analysis Scripts ####
Various Python scripts are provided for analysis of *Test 1* fits using C++ and Fortran. 

##### Comparing True, Fitted Parameter Values #####
1. Ensure Python use is enabled.
2. Generate test data, using either Fortran or C++ programme. This will generate a file containing true parameter values, for example `mp2test1_true_params_c.txt`.
3. Fit data using `./pede str.txt`. This will generate a results file `millepede.res`.
4. Run python script to compare parameter values, using `python compareParams.py -f <pede results file> -t <true parameter values file>`

##### Constraint Comparison #####
`compareConstraints.py` is designed to carry out fits of *Test 1* C++ data, varying the constraints used by `pede` for fitting. Just ensure *Test 1* is properly built, with the binary `MpTest1* available, then run the script. A plot of the differences between fitted and true plane displacements, against the parameter label for each plane displacement, will be shown, with one series for each set of constraints applied in the fit.

##### Reading Parameters Into Root #####
Two python scripts are supplied to read parameter values into a Root `TTree` - `readCParamsToRoot.py`, and `readFortranParamsToRoot.py`. These both generate `mille` data, then carry out `pede` fits using different fitting methods. Parameter values, and associated uncertainties (if available), are output in a `TTree`, with a number indexing the fitting routine used (true parameter values are denoted with `0`). These Root files are used for subsequent analysis.

##### Comparing Labelled Parameter Values #####
`compareLabelsRoot.py` creates various comparison plots between true and fitted parameter values. By default, parameters found by the inversion `pede` method are examined. Parameters must be read into a Root file using the scripts above, then the script is used as: `python compareLabelsRoot.py -i <input_file>`.

##### Comparing C++, Fortran Parameter Values #####
`compareCppFortranFits.py` creates similar plots to above, with two series - one for values found with the C++ version of *Test 1*, and one for those found with the Fortran version. Parameter values are also output in the console as a table. Two Root files must be created using the scripts above, which are then used with the script using: `python compareCppFortranFits.py -c <C++_Root_file> -f <Fortran_Root_file>`.

##### Examining Alignment Matrix Eigenvectors #####
Weak modes of alignment (overall detector shifts, shears *etc.*) may be examined by diagonalising the alignment matrix, and examining its eigenvectors. This can be done using `pede` and the `eigenSpectrumAnalysis.py`, as follows:

1. Generate `mille` data using `MpTest1` binary.
2.  Modify output steering file so that `pede` will use diagonalisation.
   * Set so that the line `method diagonalization` is the last method shown in the steering file.
3. Run `pede` using the modified steering file. A file of eigenvectors (`millepede.eve`) will be produced.
4. Run script using `python eigenSpectrumAnalysis.py -i <eigenvector_file> -o <output_directory>`

This script will plot the eigenspectrum for the matrix, showing the eigenvalue associated to each eigenvector. It will also plot values for each alignment parameter in each eigenvector, showing plane displacements and drift velocity deviations separately. All plots are output as `.pdf` files in the specified directory. Eigenvectors with low eigenvalues correspond to weak modes of alignment. Applying proper constraints should suppress these weak modes.
