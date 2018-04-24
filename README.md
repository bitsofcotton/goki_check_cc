# Goki Check
These program aims to get prepared model geometry in A still image in a deterministic way.  
Please refer older information at https://sourceforge.net/p/gokicheck/wiki/Home/ .

# How to use
Please touch Makefile for Eigen and stdc++ enabled.  
This program needs ascii raw ppm files to input/output.  
To convert image files to raw ppm, it is powerful tool that imagemagick with 'convert from.image -compress none to.ppm'. 

# Parameters
* fisheye.hh
* * z_max  : z-index resolution.
* * stp    : number of points to be used in detecting edges.
* tilt.hh
* * z_ratio : [0,1] to [0,z_atio].
* scancontext.hh
* * matchPartialPartial::ndiv    : number of divides that match angles, effects the matching errors.
* * matchPartialPartial::threshp : ratio of threshold for matched points.
* * matchPartialPartial::threshs : ratio of threshold for operator ==.

# Context
This program is inspired from re-focus photo softwares.  
And around this, there's many preceders that many approach to get bump maps with certain conditions
(such as multiple camera conditions, or, with layered objects, or, spherical, or, from movie, etcetc).
There's a defocus photo algorithms http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.100.2308&rep=rep1&type=pdf some words with googled. So it's accurate for most cameras, goki_check_cc is standing on another hypothesis that is not widely  used (in another words, some transform or special camera is needed for photos...).  
There's preceders to match 3D to 2D with many approaches. (s.t. detecting topology of junction point, or, machine learning, and so on.)  
Searching the Internet more...

# Status
Searching bone-enabled 3d model simple format -> may be .fbx format is widely used.
Writing whole to rotated partials match.
Rewriting delaunay triangulation with divide and conquerer method.

# Usage
    make tools
    
    # enlarge
    ./tools enlarge input.ppm output.ppm
    
    # detect edges with absolute differential with dft value.
    ./tools collect input.ppm output.ppm
    
    # make 2d to 3d pseudo bumpmap
    ./tools bump input.ppm output.ppm
    
    # bumpmap to .obj file.
    ./tools obj input-bump.ppm output.obj
    
    # mask bumpmapped obj files and add support.
    ./tools maskobj input-mask.ppm input.obj output.obj size-ratio thin-or-thickness-size
    
    # make tilts from original and bumpmap images.
    ./tools tilt input.ppm output-base input-bump.ppm
    
    # list matches.
    ./tools match input-matchbase.ppm output-base input-tobematched.ppm matchbase-bump.ppm tobematched-bump.ppm
    
    # list matches 2d - 3d.
    ./tools match3d input-matchbase.ppm output-base input-tobematched.obj matchbase-bump.ppm
    
    # list matches 2d - 2d with hidden 3d.
    ./tools match2dh3d input-matchbase.ppm output-base input-tobematched.ppm hidden-3dmodel.obj

# How to use as library (sample code).
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> input;
    
    #include "enlarge.hh"
    enlarger2ex<float> enlarger;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> enlarged(enlarger.compute(input, enlarger.ENLARGE_BOTH));
    
    enlarger2ex<float> detect;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> edgecollect(detect.compute(input, detect.COLLECT_BOTH));
    
    #include "fisheye.hh"
    PseudoBump<float> bump;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> bumpd(bump.getPseudoBump(input));
    
    #include "tilt.hh"
    tilter<float> tilt;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> tilted(tilt.tilt(input, bumpped, 0, 8, .95));
    
    #include "scancontext.hh"
    lowFreq<float> lf;
    matchPartialPartial<float> statmatch;
    std::vector<Eigen::Matrix<float, 3, 1> > shape0;
    std::vector<Eigen::Matrix<float, 3, 1> > shape1;
    // init shape0, shape1.
    std::vector<match_t<float> > matches(statmatch.match(shape0, shape1));
    // match operations.
    
    // If you need, please scope with namespace block.
    // but include guard definition may harms.

# Demos
https://services.limpid-intensity.info/ have a sample interface working demos.
Please bookmark output directory page after sending images, then please wait some minutes.  
N.B. a file per 5 minutes, up to 20 Mo total upload size.

# Tips
These program's enlarge is based on pseudo DFT half space plausible one.  
These program's collect is based on DFT differential.  
These program's bump assumes F=∞ graphics.   
These program's match matches with calculated pseudo z-depth.  
These program's match3d assumes input file as a bump map and .obj 3d file.  
These program's match and match3d assumes one of vertices is full and another is lowPoly.

# Specification
PseudoBump generates the bumpmap that is pseudo plausible one because of one image condition.  
This is avoidable on very wide case with multiple camera conditions,
if it fails, the image resolution or color depth resolution lacks.

# Another downloads
* https://ja.osdn.net/projects/goki-check/
* https://www.sourceforge.net/projects/gokicheck/
* https://files.limpid-intensity.info/

# Pre calculated samples
![photosample pseudo bumpmap](https://files.limpid-intensity.info/photosample-bump.jpeg)
