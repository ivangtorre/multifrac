# MULTIFRAC

MULTIFRAC is a plugin for ImageJ for multifractal analysis of 2D and 3D images.  
ImageJ can be download from its official site: https://imagej.nih.gov/ij/index.html

This plugin is released uncer license Creative Commons 4.0: Attribution   
You must give appropriate credit citing

### CITE 

* Torre, I.G., Heck R.J. & Tarquis, A.M. (2020). MULTIFRAC: An ImageJ plugin for multiscale characterization of 2D and 3D stack images.

Classes have developed during several projects and more information could be found in:

* Torre, I. G., Martín-Sotoca, JJ;  Losada, J. C., Lopez, P., & Tarquis, A. M. (2018). Multiscaling properties of soil images. biosystems engineering, 168, 133-141.
* Torre, I. G., Martín-Sotoca, J. J., Losada, J. C., López, P., & Tarquis, A. M. (2020). Scaling properties of binary and greyscale images in the context of X-ray soil tomography. Geoderma, 365, 114205.
* Torre, I. G., Losada, J. C., Heck, R. J., & Tarquis, A. M. (2018). Multifractal analysis of 3D images of tillage soil. Geoderma, 311, 167-174.


### INSTALLATION and USE
* (i) Open ImageJ2 or FIJI.
* (ii) RunHelp|Update...and choose Manage  update  sites.
* (ii)  Activate Multifrac checkbox  in  the alphabetically-sorted list of update sites.  Press OK, Apply changes 
* (iv)  Restart  ImageJ  or  FIJI.  MULTI-FRAC  will  be  integrated  under  the  plugins  menu.


### ANALYSIS
More info: https://imagej.net/Multifrac
* 2D fractal dimension analysis: Implements fractal dimension computation through box counting algorithm. It admits non squared gray images, but in that case the plugin will ask the user to resize or to cut the image with included functionalities and to binarize it. Given that fractal dimension depends on the underlaying interpretation of the image, it is possible to apply it counting black or white pixels. Then, the plugin will ask the user to restrict the calculations to some particular scales or not. Finally, it will automatically save the results and the analyzed image in sub-folders where the image was loaded.
* 2D multifractal analysis box counting: It adresses multifractal analysis with BC sliding boxes. Similar to fractal dimension analysis, the method resize or cut the image if necessary and ask if the user want to invert the image. Aditionally it includes $q$ range selection and select minimum scale size for chi(q, l). chi(q, l), tau(q), D_q and f(alpha)$representations and table results will be saved in sub-folders located in image directory.
* 2D multifractal analysis gliding box: It includes similar funcionalities than in 2D multifractal analysis box counting method, but here there resize or cutting the image is not necessary. Otherwise it is similar.
* 2D Lacunarity analysis.
* 2D configuration entropy and characteristic length.
* 3D stack fractal dimension analysis: It includes similar functionalities as in 2D fractal dimension analysis but for 3D images.
* 3D stack multifractal analysis using box counting algorithm.
* 3D stack multifractal analysis using gliding box algorithm.

### SUPPORT
There is not official support but feel free to contact me at:  
ivan.gonzalez.torre(at)upm.es

### DEVELOPMENT
This plugins has been developed with Eclipse and Mavens. 

