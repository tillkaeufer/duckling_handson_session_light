# DuCKLinG hands-on session

Below you find the step by step instructions on how to run this session.  

During this project you will fit a mock observation that was created by DuCKLinG.

The goal is to answer a few questions during the process.

The quiz can be accessed using the following link: https://www.canva.com/design/DAGe7S1YAoM/UEEmq63JZn2cyPYtBioAJQ/edit?utm_content=DAGe7S1YAoM&utm_campaign=designshare&utm_medium=link2&utm_source=sharebutton

At different steps during this tutorial, I'll ask you to go to the quiz and answer a few questions. We can discusss common difficulties at the end of the session.

## Steps

These steps guide you through the tutorial. A more detailed explanation is given below. 

- Downloading the data
- Installation
- Exploring the observations
- Running a retrieval
- Plotting the retrieval results
- Creating new input files and exploring the effects (optional)

## Detailed explanation of the individual steps

### Downloading the data
Clone this repository `git clone https://github.com/tillkaeufer/duckling_hands-on_session`.

Alternitavely, you can click on the green 'code' botton on the top right and then on 'download zip'.
If you do this you need to unpack your zip file.

Since the slab grids are part of this repository, it has a size of a few hundred megabits. 
Make sure you have a good internet connection when downloading.

At the end of this step you should have somewhere a folder called 'duckling_hands-on_session'.

### Installation
We will be working in a Python 3 environment and need a few python packages.
First, make sure that you have python 3 running. If you are using conda and want to create a new environment you can do this with the follwoing command:
`conda create --name duckling_handson python=3.9`  
`conda activate duckling_handson`

Second, go in the terminal to the 'duckling_hands-on_session' folder and run the following command:
`pip install -r requirements.txt` 

This should hopefully install a bunch of packages.

You can test if everything is installed correctly by running:
`python test_installation.py` 

### Exploring the observations

At this point you can have a look at the individual observations that are available for fitting.  
You find them in the 'Observations' folder. 
Look at the images and try to answer the first four questions of the [quiz](https://www.canva.com/design/DAGe7S1YAoM/UEEmq63JZn2cyPYtBioAJQ/edit?utm_content=DAGe7S1YAoM&utm_campaign=designshare&utm_medium=link2&utm_source=sharebutton).

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


