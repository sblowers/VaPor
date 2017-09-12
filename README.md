# VaPor
Vascular Porous Model for Simulating Biological Temperatures

Authors: Stephen Blowers, Ian Marshall, Michael Thrippleton, Bridget Harris, Peter Andrews, Iain Bethune, and Prashant Valluri.

![vaporimage1](https://cloud.githubusercontent.com/assets/25664298/22834082/eba3f48e-efab-11e6-8cfc-6f906a6092a4.png)
![vaporimage2](https://cloud.githubusercontent.com/assets/25664298/22835187/cc857236-efaf-11e6-8a37-0c2dfe847d4b.png)

The Vascular Porous (VaPor) model is an open source software designed to simulate biological mass flows and heat transfer. Blood vessels are modelled as 1D line segments embedded within a 3D domain representing tissue. All mass and energy balances between these domains are conserved. The current version is set up to focus on cerebral temperatures. It includes domain reading, vessel generation, as well as flow and temperature solvers. Additionally, some display features are included to help read the results.

Note: this software requires the Imaging Processing Toolbox in Matlab to work.

Note: The linear solver in this program can use a considerable amount of RAM to solve. It is not recommended to use the DomainCoarsening = 'none' option with the example data, nor generate more than 200,000 vessels. Additionally, the program sometimes cannot be aborted using Ctrl+c due to the large memory usage whilst the linear solver is in progress and the only way to exit during the run is to kill Matlab itself. 

A brief overview of the the initialisation functions is given in the .pdf ProgramDocumentation. 


Data for the example arterial vessel tree provided in Vessels\BG001.swc is taken from the Brain Vascular (BraVa) database located at http://cng.gmu.edu/brava [1]. The specific reconstruction is BG001. 

Data for the example venous vessel tree has been extracted from a MRI scan provided by the Brain Research Imaging Centre, Edinburgh Imaging, University of Edinburgh.

The example image sequences used are taken from the tissue probability maps provided by SPM12 software located at http://www.fil.ion.ucl.ac.uk/spm/. To incorporate these files into the solver, they first must be extracted from the repository within the SPM12 files found on the above website. The raw image files are stored in spm12/tpm/TPM.nii and must be extracted through an appropriate image processeing software (e.g. ImageJ found free at https://imagej.nih.gov/ij/download.html). The corresponding directories for the image files should then be specified within "domainReadingInitialisation". Alternatively, the image files can be found hosted at https://github.com/sblowers/ImageFiles under a GNU Public License v2.0.

These are supplied for the purposes of exhibiting the application of the software to cerebral temperatures. They do not represent work developed by the authors of this software. Alternative data inputs can be used within the model. 

[1] Susan N. Wright, Peter Kochunov, Fernando Mut Maurizio Bergamino, Kerry M. Brown, John C. Mazziotta, Arthur W. Toga, Juan R. Cebral, Giorgio A. Ascoli. Digital reconstruction and morphometric analysis of human brain arterial vasculature from magnetic resonance angiography. NeuroImage, 82, 170-181, (2013). http://dx.doi.org/10.1016/j.neuroimage.2013.05.089
