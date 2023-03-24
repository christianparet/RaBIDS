# raBIDS
rapid, automated conversion to Brain Imaging Data Structure (BIDS) 

## Welcome!
Thank you for looking at this code repository! I hope you find the content helpful.

## Description
raBIDS is Matlab-based software supporting the import of human neuroimaging data to the BIDS format (visit https://bids.neuroimaging.io/ for information on BIDS). The 'Import' part of raBIDS is a wrapper built around the dicm2nii toolbox (https://github.com/xiangruili/dicm2nii). The 'Create-SOTs' part provides an automated routine to read stimulus onset times from an experiment logfile. The vision of this project is to provide a tool to researchers without programming background who are used to work in a Windows environment and with MS Office applications. In its current version (v0.X), RaBIDS uses Excel-2010-based tables to load experiment information that is used to assign scans to tasks.
There are several limitations and software requirements that may prevent users from using the current beta version v0.2.X. Please consult the manual for more information.
This project is currently under development. Before using raBIDS please consult the manual. A tutorial is provided. Data for this tutorial can be requested from the author.
Make sure to use the latest release.

## On a personal note
I am fascinated by the idea that all neuroimaging data is archived following a generic structure such as BIDS. A generic data structure can massively improve data re-use, data sharing, and reproducibility of results. Conversely, it is a powerful instrument of quality assurance for researchers, as it aids to keep overview over your data and it helps you and your colleagues to cooperate on a dataset. Besides, the plenty tools that are available to work with BIDS data are a valuable resource for neuroimaging researchers and a great chance for the field to improve research quality.

All the taken advantages aside, it can be a painstaking experience to implement BIDS in the lab. I had to learn this myself as a lay programmer at a medical faculty with limited open-science infrastructure, by that time. No question, our research IT is great and very helpful. But when I started with BIDS, the tools that were available for BIDS conversion (i.e., the import of fMRI data acquired at the scanner into the BIDS format) required advanced programming skills, or where not so easy to run within our institutional IT infrastructure. This is where I started off to write the code in this repository.

Eventually, the code works well for my lab, but it may not work for you if you are using a different operating system and different kind of data than the one it was developed on. Also, you will miss features that are important for you but that were not (yet) required for data analysis in my lab.

If you find the code helpful and you are excited about the idea behind this endeavour, you are very warmly invited to contribute to the project. The current version is miles away from the vision of a easy to use, generic toolbox for BIDS conversion of fMRI data. But it is a first step. Maybe it can attract other researchers who are more skilled in programming, who share this vision, and who are willing to collaborate.

## Issue reporting
Please report issues via Github.

## License
The program is provided under the CC BY 4.0 license. The developer does not guarantee functioning of the code. There is no guarantee that the results obtained with this code are valid.

## Author information
Christian Paret

RG Psychobiologie of Selfregulation

Department of Psychosomatic and Psychotherapeutic Medicine

Central Institute of Mental Health

christian . paret [at] zi-mannheim . de
