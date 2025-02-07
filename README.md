# Ringberg DuCKLinG session

Below you find the step by step instructions on how to run this session.  

During this project you will fit a mock observation that was created by DuCKLinG.

The goal is to answer the following questions:
- Which molecules are included in the fitting?
- How would you add another molecule?
- How would you change the prior ranges?
- What dust composition is retrieved?
- What molecular conditions are retrieved?
- How much noise was added to the mock observation?

## Steps

These steps guide you through the tutorial. A more detailed explanation is given below. 

- Downloading the data
- Installation
- Running a retrieval
- Plotting the retrieval results
- Creating new input files and exploring the effects (optional)

## Detailed explanation of the individual steps

### Downloading the data
Clone this repository `git clone https://github.com/tillkaeufer/Ringberg_duckling_session`.

Alternitavely, you can click on the green 'code' botton on the top right and then on 'download zip'.
If you do this you need to unpack your zip file.

At the end of this step you should have somewhere a folder called 'Ringberg_duckling_session'.

### Installation
We will be working in a Python 3 enviroment and need a few python packages.
First, make sure that you have python 3 running.
Second, go in the terminal to the 'Ringberg_duckling_session' folder and run the following command:
`python requirements.txt` (ish)

This should hopefully install a bunch of packages.
If you have any problems with individual packages...

### Running a retrieval

A single command will start the retrieval (hopefully).
Before you execute it have a look at the input file that is used (found in ...).
Try to understand what the individual lines in the input files are doing and answering the following questions:

- What distance to the object is used for the fitting?
- Which parameters are fixed during the fitting?
- Which dust species are included in the fitting?
- Which molecules are included in the fitting?
- How would you add another molecule?
- How would you change the output folder that the results are saved in?
  
If you feel comfortable with the file let's run a retrieval.  
For doing so, make your that you are in the 'Ringberg_duckling_session' folder and run the following command:

`python retrieval-input.py ./Input_files/example1.txt`

Now you to be a bit patient depending on your machine.
In the meanwhile you can add your email address to ... if you want to be included in the DuCKLinG repository to use it in your future research.


### Plotting the retrieval results

When the retrieval finished you can plot the results using the following command:

`python retrieval-input.py ./Input_files/example1 all`

This can again take a little while.  

Afterwards there should be a folder called './Output/example1/figures'.  
Open this folder and explore the figures that where created.  
Use them to answer the following questions:

- How large are the deviations between model posterior and observation (estimate based on figure)?
- What is the dominant dust species that has been retrieved?
- What temperatures have been retrieved for H2O?
- What column densities have been retrieved for H2O?
- What molecule is causing the features at ...?
- Which molecule is emitting at smaller radii (CO2 or H2O)?

If you made it this far, **congratulations**!!!  
Feel free to continue playing around by changing things in the input file and exploring the effect on the fit.  

The optional things I planed are listed below.

### Optional next steps

#### How to model dust absorption features?

#### How many water components where used here?

##### How to model gas absorption features?

##### Fitting the continuum subtracted spectrum


