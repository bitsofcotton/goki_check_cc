# Goki Check
These program aims to get context(s) from a still image.
Now, only pseudo 2D->3D and tilt program is published. 

Please refer older information at https://sourceforge.net/p/gokicheck/wiki/Home/ .

# How to use
Please touch Makefile for Eigen and stdc++ enabled.
This program needs ascii raw ppm files to input/output.
To convert image files to raw ppm, it is powerful tool that imagemagick with 'convert from.image -compress none to.ppm'. 

# Parameters
* fisheye.hh
** z_max  : z-index resolution.
** stp    : number of points to be used in detecting edges.
** nlevel : ratio to auto level.
* tilt.hh
** z_atio : [0,1] to [0,z_atio].
** psi    : tilt as &pi;/2 * psi.
** samples: number of rotation samples.

# Context
This program is inspired from re-focus photo softwares.

# Status
Deducting 3D to 3D match with many different origin rotations. 

# Usage
    make tools
    
    # enlarge
    ./tools enlarge input.ppm output.ppm
    
    # enlarge for left differential
    ./tools enlargeds input.ppm output.ppm
    
    # detect edges with absolute differential with dft value.
    ./tools collect input.ppm output.ppm
    
    # make 2d to 3d pseudo bumpmap
    ./tools bumpscale input.ppm output.ppm
    
    # make tilts from original and bumpmap images.
    ./tools tilt input.ppm output-base input-bump.ppm

# How to use as library (sample code).
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> input;
    
    #include "enlarge.hh"
    enlarger2ex<float, complex<float> > enlarger;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> enlarged(enlarger.enlarge2(input, enlarger2ex<float, complex<float> >::ENLARGE_BOTH));
    
    enlarger2exds<float, complex<float> > enlargerds;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> enlargedds(enlargerds.enlarge2ds(input, enlarger2exds<float, complex<float> >::ENLARGE_BOTH));
    
    #include "edgedetect.hh"
    edgedetect<float, complex<float> > detect;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> edgecollect(detect.detect(input, edgedetect<float, complex<float> >::COLLECT_BOTH));
    
    #include "fisheye.hh"
    PseudoBump<float> bump;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> bumpped(bump.getPseudoBump(input));
    
    #include "tilt.hh"
    tilter<float> tilt;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> tilted(tilt.tilt(input, bumpped, 0));
    
    // If you need, please scope with namespace block.
    // but include guard definition may harms.

# Demos
http://services.limpid-intensity.info/ have a sample interface to working demos.
Please bookmark output directory page after sending images, then please wait some minutes.  
N.B. a file per 5 minutes, sha 256 order.  
N.B. up to 20 Mo total upload size.

# Tips
These program's enlarge and enlarge-ds is based on pseudo DFT half space plausible one.  
These program's collect is based on DFT differential.  
These program's bump assumes F=∞ and y-axis polarized graphics.   

# Another downloads
* https://ja.osdn.net/projects/goki-check/
* https://www.sourceforge.net/projects/gokicheck/
* https://files.limpid-intensity.info/goki_check_cc-1.00-rc4.tar.gz (preparing...)
