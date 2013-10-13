#Orbital parameters file

object='Theta Ori E, comp 2'    # Name of the object
file='ResG5C2D.txt'             # Name of the input file
psname='orie_1172'              # Name of output plot file 
output='pdf'                    # Can choose between: pdf, eps or png
data=False                      # If True then exit data will put in file *.txt
plot=True                       # Plot in Python window
plotlim=1.3                     # Plot limits. 1 = close fit.
porb=9.895                      # Period in days
hjd0=53284.9766                 # HJD Zero point for fit
gama=35.0                       # Systemic Velocity in km/s
k1=84.0                         # Semi-amplitude velocity in km/s
phaseoff=0.0001                 # Phase Delay between photo and Spectra (Sepctroscopic CVs)
sigma=3.2                       # Default sigma if errors are unknown. For estimating chi2
fix_porb=True                   # Change to False to fix parameter
fix_hjd0=True                   # Change to False to fix parameter
fix_gama=True                   # Change to False to fix parameter
fix_k1=True                     # Change to False to fix parameter
fix_phaseoff=False              # Change to False to fix parameter
#For the looping program
variable='porb'                 # Variable to loop: porb, hjd0, gama, k1
initial=9.600                    # Initial value
delta=0.008                      # Step value for loop
num=100                           # Number of loops