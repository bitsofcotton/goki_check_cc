# goki_check_cc
Goki Check on the implementation of C plus plus.

# How to use.
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
