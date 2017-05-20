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
    ./tools bump input.ppm output.ppm
    
    # make 2d to 3d pseudo bumpmap on dft differential.
    ./tools bumpds input.ppm output.ppm
    
    # make tilts from original and bumpmap images.
    ./tools tilt input.ppm output-base input-bump.ppm
