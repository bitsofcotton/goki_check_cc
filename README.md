# Goki Check
These program aims to get context(s) from a still image.  
And then, aiming to get scene from roughly prepared dictionaries, and, reproduce the scene with another prepared dictionaries, and, get a copy of image produced 3D model into a dictionary.  
Now, implementing bone enabled match.

Please refer older information at https://sourceforge.net/p/gokicheck/wiki/Home/ .

# How to use
Please touch Makefile for Eigen and stdc++ enabled.  
This program needs ascii raw ppm files to input/output.  
To convert image files to raw ppm, it is powerful tool that imagemagick with 'convert from.image -compress none to.ppm'. 

# Parameters
* fisheye.hh
* * z_max   : z-index resolution.
* * stp     : number of points to be used in detecting edges.
* * crowd   : box size for reducing crowded point return.
* * vmax    : max points to be returned in pseudoBumpVec.
* * nloop   : number of loops with tilting and bumps.
* * ndiv    : tilt angle base.
* * rdist   : ratio of distance for camera and plane.
* tilt.hh
* * z_ratio : [0,1] to [0,z_atio].
* scancontext.hh
* * matchPartialPartial::ndiv    : number of divides that match angles.
* * matchPartialPartial::thresh  : threshold to detect parallel vectors.
* * matchPartialPartial::threshp : threshold for matched points.
* * matchPartialPartial::threshr : threshold for matched size ratios.
* * matchPartialPartial::threshs : threshold for operator ==.

# Context
This program is inspired from re-focus photo softwares.  
And around this, there's many preceders that many approach to get bump maps with certain conditions
(such as multiple camera conditions, or, with layered objects, or, spherical, or, from movie, etcetc).
(And I have no clue which is a correct official page, no hyperlink from here.).   
This get bump maps in a different way.

# Status
Searching bone-enabled 3d model simple format. Writing whole to rotated partials match.

# Usage
    make tools
    
    # enlarge
    ./tools enlarge input.ppm output.ppm
    
    # enlarge for left differential
    ./tools enlargeds input.ppm output.ppm
    
    # detect edges with absolute differential with dft value.
    ./tools collect input.ppm output.ppm
    
    # make 2d to 3d pseudo bumpmap
    ./tools bump input.ppm output.ppm output-delaunay.ppm output.obj
    
    # make tilts from original and bumpmap images.
    ./tools tilt input.ppm output-base input-bump.ppm
    
    # make (pseudo) lowpoly and get match.ppm and .obj file.
    ./tools lpoly input.ppm output-match.ppm output.obj
    
    # list matches.
    ./tools match input-matchbase.ppm output-base input-tobematched.ppm
    
    # list matches 2d - 3d.
    ./tools match3d input-matchbase.ppm output-base input-tobematched.obj

# How to use as library (sample code).
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> input;
    
    #include "enlarge.hh"
    enlarger2ex<float> enlarger;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> enlarged(enlarger.enlarge2(input, enlarger2ex<float>::ENLARGE_QUAD));
    
    enlarger2exds<float> enlargerds;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> enlargedds(enlargerds.enlarge2ds(input, enlarger2exds<float>::ENLARGE_QUAD));
    
    #include "edgedetect.hh"
    edgedetect<float> detect;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> edgecollect(detect.detect(input, edgedetect<float>::COLLECT_BOTH));
    
    #include "fisheye.hh"
    PseudoBump<float> bump;
    std::vector<Eigen::Matrix<float, 3, 1> > points;
    std::vector<Eigen::Matrix<int,   3, 1> > delaunay;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> bumps;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> bumpd(bump.getPseudoBumpVec(input, points, delaunay, bumps));
    
    #include "tilt.hh"
    tilter<float> tilt;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> tilted(tilt.tilt(input, bumpped, 0, 8, .95));
    
    #include "scancontext.hh"
    lowFreq<float> lf;
    matchPartialPartial<float> statmatch;
    std::vector<Eigen::Matrix<float, 3, 1> > shape0(lf.getLowFreq(input, 300));
    std::vector<Eigen::Matrix<float, 3, 1> > shape1(lf.getLowFreq(input2, 300));
    std::vector<match_t<float> > matches(statmatch.match(shape0, shape1));
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
These program's match3d assumes input file as bump map and .obj 3d file.  
These program's match assumes one of vertices is full and another is lowPoly, but lowFreq implementation now, it worse generate lowPoly.

# Known bugs
PseudoBump makes a pseudo plausible things.

# Another downloads
* https://ja.osdn.net/projects/goki-check/
* https://www.sourceforge.net/projects/gokicheck/
* https://konbu.sakura.ne.jp/files/goki_check_cc-1.01-lack-rotate-release3.tar.gz
* http://files.limpid-intensity.info/goki_check_cc-1.01-lack-rotate-release3.tar.gz
