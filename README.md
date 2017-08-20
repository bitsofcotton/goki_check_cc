# Goki Check
These program aims to get context(s) from a still image.
Now, testing solid whole to rotate partials match, and, deducting z-axis variable scaled one.

Please refer older information at https://sourceforge.net/p/gokicheck/wiki/Home/ .

# How to use
Please touch Makefile for Eigen and stdc++ enabled.  
This program needs ascii raw ppm files to input/output.  
To convert image files to raw ppm, it is powerful tool that imagemagick with 'convert from.image -compress none to.ppm'. 

# Parameters
* fisheye.hh
* * z_max  : z-index resolution.
* * stp    : number of points to be used in detecting edges.
* * nlevel : ratio to auto level.
* tilt.hh
* * z_atio : [0,1] to [0,z_atio].
* scancontext.hh
* * matchPartialPartial::thresh  : threshold to detect parallel vectors.
* * matchPartialPartial::threshp : threshold for matched points.
* * matchPartialPartial::threshr : threshold for matched size ratios.

# Context
This program is inspired from re-focus photo softwares.

# Status
Writing whole to rotated partials match.

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
    
    # list matches.
    ./tools match input-matchbase.ppm output-base input-to-bematched.ppm ref-to-be-matched.ppm ref-matchbase.ppm

# How to use as library (sample code).
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> input;
    
    #include "enlarge.hh"
    enlarger2ex<float> enlarger;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> enlarged(enlarger.enlarge2(input, enlarger2ex<float>::ENLARGE_BOTH));
    
    enlarger2exds<float> enlargerds;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> enlargedds(enlargerds.enlarge2ds(input, enlarger2exds<float>::ENLARGE_BOTH));
    
    #include "edgedetect.hh"
    edgedetect<float> detect;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> edgecollect(detect.detect(input, edgedetect<float>::COLLECT_BOTH));
    
    #include "fisheye.hh"
    PseudoBump<float> bump;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> bumpped(bump.getPseudoBump(input));
    
    #include "tilt.hh"
    tilter<float> tilt;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> tilted(tilt.tilt(input, bumpped, 0, 8, .999));
    
    #include "scancontext.hh"
    lowFreq<float> lf;
    matchPartialPartial<float> statmatch;
    std::vector<Eigen::Matrix<float, 3, 1> > shape(lf.getLowFreq(input, 300));
    std::vector<Eigen::Matrix<float, 3, 1> > shape2(lf.getLowFreq(input2, 300));
    statmatch.init(shape0, .85, .25, .125);
    std::vector<match_t<float> > matches(statmatch.match(shape1, 20));
    // match operations.
    
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
These program's match assumes input file as bump map. And, that matches including z-depth.
These program's match3d assumes input file as bump map and .obj 3d file. And for now, it considers only vertexes (not planes). So redig result gains strange one.

# Another downloads
* https://ja.osdn.net/projects/goki-check/
* https://www.sourceforge.net/projects/gokicheck/
* https://konbu.sakura.ne.jp/files/goki_check_cc-1.00-stable.tar.gz
* http://files.limpid-intensity.info/goki_check_cc-1.00-stable.tar.gz
