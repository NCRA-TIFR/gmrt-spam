# gmrt-spam
An installation of SPAM for reducing GMRT observations using Direction Dependent Ionospheric Calibration (customized from SPAM implementation by Huib Intema)

## Downloading and installing SPAM

Here are some instructions on how to install SPAM on your Linux 64-bit system. 

Download this repository to your local device.

Make sure that the following dependencies have been met. If not install, and then proceed with the installation procedure.

* Python version 2.7 (plus C include files)
* Python packages: pylab / matplotlib, numpy, scipy
* gcc
* make
* swig
* optional: mencoder / mplayer (for making phase screen movies)
* optional: ImageMagick convert (for making phase screen movies)

**For the rest of the installation procedure we will assume ROOTPATH to be the path to the gmrt-spam folder**

## Edit /spam/setup.sh file

1. Open the setup.sh file using your preferred text editor.
2. Modify the SPAM_PATH variable to point to SPAM folder (i.e. ROOTPATH/spam
3. Modify the SPAM_HOST variable to contain the hostname of the device (can be found out by running the command 'hostname')
4. Modify the variable PYTHON to point to the location of the Python2.7 executable (can be found out by running the command 'whereis python')

## Install AIPS

1. Change the working directory to 'ROOTPATH/spam/'
2. Extract the AIPS package (run the command "tar -xvf 'name of the AIPS package'")
3. Run the command 'perl install.pl -n'

To quickly get through the AIPS installation:
- screen 0: <enter>
- screen 4: <enter>
- screen 4b: <enter>
- screen 5: <e.g., your institute name in capitals (no spaces), and enter>
- screen 5a: <enter>
- screen 5b: <enter>
- screen 6: <enter>
- screen 7: <copy/paste suggested path and enter>
- screen 8: <2x enter>
- screen 9: <enter>
- screen 9B: <enter>
- screen 11: <2x enter>
- during installation: <3x enter>

4. We will now increase the number of interactive AIPS sessions to 16. Type the following command 'source LOGIN.SH'
5. Run the following command 'RUN SETPAR' 
6. Within SETPAR

- 2 <enter>
- 10 <enter>
- 16 <enter>
- -1 <enter>
When prompted to enter passsword:
- AMANAGER <enter>
- 4 <enter>

7. Test AIPS by running aips tv=local.
8. Press enter in the terminal, and then enter a few test commands i.e.

- 11 <enter>
- print 2+2 <enter>
- kleenex <enter>

## Install ParselTongue

1. Change working directory to 'ROOTPATH/spam/parseltongue-2.3a'
2. Run the command 'sh ./configure --prefix=${SPAM_PATH}/ParselTongue --with-obit=${SPAM_PATH}/Obit PYTHON=${PYTHON}'
*This should return without errors. If it does return with errors, there is most likely to be a problem with the Obit installation. (See below)*
3. Run the following commands
- make
- make install

## Install SPAM

1. Change working directory to 'ROOTPATH/spam/python/spam/'
2. Open the 'makefile' using your preferred text editor
3. Edit the makefile so that the SWIGFLAGS and CCFLAGS point to the correct Python C header file directory
4. Run make
