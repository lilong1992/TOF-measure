This program is used to take a quick look at neutron time of flight (TOF), and if distance (dist) and gamma time (tgam) are provided, a quick neutron energy plot (h12) will be viable.

Function triggertimefastcal() can be used separately to find PMT pulse time. The method implements constant fraction discrimination (CFD) to find 50% of peak as the onset time of a pulse. Check: https://www.sciencedirect.com/science/article/pii/S0168900207006213 for more mathmatical info. If the pulse baseline is not well defined (e.g. not enough pretrace), the default 'subtract' should set to false, and an estimate of half baseline should be provided by the user.

The TOF program will need MGM's TWaveForm to work. To install Twaveform, see https://github.com/lilong1992/TWaveform

The program takes the rootified data from SIS3316 digitizer. Assuming channel 0 is the PMT signal and channel 4 is the BPM signal. To run the program, simply type beamenergytof(filename) in ROOT. Sample data is also provided.

