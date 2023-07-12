# DNA Strand Displacement System

## Data Process

The data process was used to process the data generate from the "Q" qpcr machine. The data from the machine is stored in a csv file formate with the general naming convension of "4WJ_HEX\_{test-type}\_{plate number}\_{trigger}\_{trigger}". There are slight variation for the naming convention depending on the test type.

### General processing

The data processing script will detemine the test gorups and condtions based on the file name and search the input files for any csv files that matches the search string. 

## Test type present

There are currently three test types: concentration (Conc), ratio (Ratio), and screening (screen)

### Concentration study

The concentration study is used to examine how the trigger concentration affects the release profile.

### Screening study

The screening study is used to examine how the trigger design affects the release profile. The script in screen.py imports some custom scripts used to deteremine the test groups and conditions. It first search the input folder for any csv files that matches the prescribed search string ("4WJ_HEX_Screen*.csv") and combined them into a master list of data file names. The data file name master list is then processed to generate a dictionary of unique test conditions and presented to the user for selection. If there is no file that matches the search string, the script will end.

After the user selected the desired condition for processing, a list of different paths for downstream file storage was generated. Lists of group (numbers and names) and a list of experimental time point were also generated. 

Data combination was carried out by searching the directory for file matching the testing condition and convert the desired files into a dataframe for further processing. The system will determine the location of both postsive and negative control for the trial and subtract the average negative control from all conditions. The data is then normalized against the average positive control condition. The samples that does not start or is negative after the background subtraction was omitted from the dataframe. The dataframes from each trial are combined into one dataframe and export into a pickle file format for a fast retrival and preserved the format within Python. Prior the the export, the script will check to see if there is an existing pickle file with the same condition and if those two dataframe is exactly the same. If there is an existing pickle file and they are not the same dataframe, it will concatenate the two dataframes together. This check is to ensure that the file is not override if there is a new trial added after the initial trials.

Since the outliers were dropped in the data combination process, new group lists were generated to update the list of samples with the sample that are still presented in the dataframe. The dataframe is further normalized against the maximum values of the averaged 100% release condition, which is the perfectly complementory trigger. The dataframe from this operation is also exported into a seperate pickle file format for a fast retrival and preserved the format within Python. 

The dataframe from the normalization step will be further analysis by calculating the average and standard deviation (std) for each condition. The column name in this dataframe correspondes to the condition and what the value is (mean or std). The dataframe is again exported into a seperate pickle file format for a fast retrival and preserved the format within Python. 

The data files that are processed under this will be moved to a seperate folder after the whole script is completed.

### Ratio study

The ratio study is used to examine how different ratio of triggers affects the release profile.

## Curve fitting

The curve fitting was carried out using the lmfit package using two models: one phase association and a simplified one step, irreversible second order rate law kinetic. Model class from lmfit was used for the weighted curve fitting. The model takes in the time points (t), mean release quantities (y), and the standard deviations for the weighted fitting. Those values for each condition were calculated prior in the data processing step.

The one phase association fitting was a classic curve fitting procedure by adjusting the two parameteres, platuea and rate, to minimize the error between the predicted values from the observed values. The kinetic model was similar with an extra step of solving the series of ordinary differential equations as a initial value problem. This is acheived by using the solve_ivp funciton in Scipy.integrate with LSODA method. The initial condition was set to be 500 nM of initial species and  0 nM of products to be consistant with the experimental procedures. The 3-$\sigma$ uncertainty bands for both models were calculated using the eval_uncertainty function from lmfit. The parameters for both models were extract along with the goodness-of-fit value ($r^2$) into a seperate dataframe and exported into a pickle file for a fast retrival and preserved the format within Python. 

## Nupack simulation
Concurrent to the experiemental data analysis, Nupack (v4.0.0.27) was downloaded and used as a package with Python3. Nupack was used to simulate the theoratical concetnration for each component in the system at the specific ion concentrations (monovalent ion conc. is 161.5 mM  and divalent ion conc. is 1.2 mM). The monovalent ion includes sodium and potassium; divalent ion includes magnesium. The ion concentration calculation was based on the buffer condition used in the experiment. The temperature was set to be 40 $^\circ C$ for the Nupack model to be consistant with the experimental condition. The stacking option of the ensemble was selected for the Nupack model. The tube class was defined with three DNA strands, namely the base, incumbant (incumb), and trigger (trig), at the starting concentration of 500 nM. The max complexes size was set to be 2 for the tube class. Both complex analysis and complex concentration function from Nupack with the tube class and the model specified were used to calculate the structure free energy and simulate the theoratical complex concentration of each complex. The free energy and the theoratical concentration were combined into a dataframe and exported into a pickle file for a fast retrival and preserved the format within Python. 