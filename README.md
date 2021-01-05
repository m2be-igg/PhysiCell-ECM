# PhysiCell-ECM: A PhysiCell extension to account for the extracellular matrix

> This extension was developed for PhysiCell 1.7.1 (most recent at the time of publication)

## Overview:
This extension aims to implement the effect of the mechanical properties of the extracellular matrix (ECM) into the PhysiCell framework [1]. To do so, we have extended the PhysiCell code to take into account the viscosity of the ECM on individual cell migration.

We have used previously published experimental data [2], obtained for collagen matrices of different collagen densities, to characterize the effect of matrix density on both single cell motility and its subsequent effect The mechanical properties of the extracellular matrix for the collagen matrices used experimentally have also been characterized previously [3].

Moreover, we have defined a new function to randomly generate cell-generated locomotive force values, to fit the migration patterns observed in [2].

### References
[1] Ghaffarizadeh, Ahmadreza, et al. "PhysiCell: an open source physics-based cell simulator for 3-D multicellular systems." PLoS computational biology 14.2 (2018): e1005991.

[2] Plou, J., et al. "From individual to collective 3D cancer dissemination: roles of collagen concentration and TGF-β." Scientific reports 8.1 (2018): 1-14

[3] Valero, Clara, et al. "Combined experimental and computational characterization of crosslinked collagen-based hydrogels." PLoS One 13.4 (2018): e0195820.

## New Makefile rules:

**make ECM-single**: populates a project for single cell motility.

**make ECM-multi**: populates a project for cluster growth.

## Documentation
For general information on PhysiCell, look at QuickStart.pdf and UserGuide.pdf in the documentation/PhysiCell folder.
For more information on the data analysis procedures used to study the effect of collagen density on cell motility and cluster growth, see the Jupyter Notebooks in the documentation/data_analysis folder.
