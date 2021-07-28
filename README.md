# Pulseq sequences for field-monitoring

Includes sequences:

 - Off-resonance and position calibration
 - Local eddy current calibration

Please visit https://github.com/pulseq/pulseq for further information about Pulseq.

## Getting started

### Clone this repositiory 

    git clone https://github.com/SkopeMagneticResonanceTechnologies/Pulseq-LocalEddyCurrentCalibration.git

Add Pulseq as a submodule

    cd .\Pulseq-LocalEddyCurrentCalibration\
    git submodule add https://github.com/pulseq/pulseq
    
### Creating the Pulseq sequence files in MATLAB  

 - Open MATLAB 
 - Run offres_pos_calib.m to create the sequence file for the off-resonance and position calibration
 - Run local_ec_calib.m o create the sequence file for local eddy current calibration

Note that the polarity of the x-axis has been reversed to be played out correctly on Siemens scanners.
