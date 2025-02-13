### LineData

In this folder you find the slab grids that are used to interpolate and integrate.  
The grid present here are very small (to decrease the used storage and make the fitting faster).  
All grids here span a temperature range from 200K to 800K and a column density range from $10^15 \rm cm^(-2)$ to $10^18 \rm cm^(-2)$.  

The binned_data folder will contrain all the binned molecular data.  
This means that everytime you are fitting a new observation, the programm will search this folder and check if the observation has been fitted previously.  
In that case the slab grids are not rebinned to the observation, but simply the previously binned data is used.  
This can make things much faster in case you are fitting the same observations over and over again. 

