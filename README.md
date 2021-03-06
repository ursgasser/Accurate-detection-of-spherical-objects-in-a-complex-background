# Accurate-detection-of-spherical-objects-in-a-complex-background

IDL code for the detection of spherical particles in 3d microscopy images.
This code was published in (see citation.bib) U. Gasser and B. Zhou, OpticsExpress, DOI: https://doi.org/10.1364/OE.434652 (2021).

We plan to also give a python version.

The particle detection code is given in the files find_spheres_div.pro and find_spheres_y2.pro.

find_spheres_div.pro: B&gradient method explained in the paper mentioned above.

find_spheres_y2.pro: B&Y method explained in the paper mentioned above.

See the IDL code in the files 'img_4769_example.txt' and 'img_4400_example.txt' for particle-detection examples. Calculated images are used in these two examples and the determined positions can be compared with the known positions of the spheres. Both the B&$\nabla$ and B&Y methods are used with and without fracshift and suppression of bridges to compare methods. 
An example of the B&Y method with a measured confocal image is given in file 'img_asm533_example.txt'.
