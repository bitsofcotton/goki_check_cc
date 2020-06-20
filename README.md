# Goki Check
These program aims to implement one of a complement to ongoing other utilities.

# How to use
Please touch Makefile for libc++ enabled.  
This program needs ascii raw ppm files to input/output.  
For speed, we need http://eigen.tuxfamily.org/ library.
We need imagemagick for normal use, and ffmpeg for making .mp4.  
And, this command needs CONTINUOUS input. So without them, this program returns bugly result, the same reason causes pextend or enlarge recursive returns bugly.

# Parameters
* redig.hh
* * vbox : size of vector gathering rectangle.
* * rz   : z-axis output ratio.
* match.hh
* * ndiv    : number of divides that match angles, effects the matching errors.
* * threshr : error tolerance rate for matching, smaller is tight.
* * threshp : ratio of threshold for number of matched points.
* * threshs : ratio of threshold for operator ==.

# Context
This program is inspired from re-focus photo softwares.  
And around this, there's many preceders that many approach to get bump maps with certain conditions
(such as multiple camera conditions, or, with layered objects, or, spherical, or, from movie, etcetc).
There's a defocus photo algorithms http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.100.2308&rep=rep1&type=pdf some words with googled. So it's accurate for most cameras, goki_check_cc is standing on another hypothesis that is not widely used (in another words, some transform or special camera is needed for photos...). And, there exists preceders that make it from multiple pint images, this makes very accurate results.  
There's preceders to match 3D to 2D with many approaches. (s.t. detecting topology of junction point, or, machine learning, and so on.). And it is fater than this that PnP problem and specific point based matching.  
By searching with some word that is not common, there exists the article https://ryo620.org/2018/02/to-gltf-from-fbx-by-blender/ that I firstly know the gltf format by this. There's a https://github.com/jessey-git/fx-gltf/ library, but compatibility for this is abandoned.  
From some news, there exists the machine learning method that is making 3D model from single picture condition.
Searching the Internet more...

# Status
Archived.

# Usage
    make gokicheck
    
    gokicheck (collect|sharpen|bump|enlarge|pextend) <input.ppm> <output.ppm>
    gokicheck pcopy <vbox> <thresh> <zratio> <num_of_emph> <outbase> <inputdst.ppm> <inputdst-bump.ppm> <inputsrc.ppm> <inputsrc-bump.ppm>
    gokicheck ppred <vbox> <thresh> <zratio> <num_of_emph> <outbase> <input0.ppm> <input0-bump.ppm> ...
    gokicheck pred  <output.ppm> <input0.ppm> ...
    gokicheck obj   <gather_pixels> <ratio> <zratio> <thin> <input.ppm> <mask.ppm>? <output.obj>
    gokicheck (tilt|sbox)    <index> <max_index> <psi> <input.ppm> <input-bump.(ppm|obj)> <output.ppm>
    gokicheck (match0|match) <num_of_res_shown> <num_of_hidden_match> <vbox_dst> <vbox_src> <zratio> <dst.ppm> <src.ppm> <dst-bump.(ppm|obj)> <src-bump.(ppm|obj)> (<dst-mask.ppm> <src-mask.ppm>)? <output-basename>
    gokicheck matcho  <match> <nemph> <vbox_dst> <vbox_src> <zratio> <dst.ppm> <src.ppm> <dst-bump.(ppm|obj)> <src-bump.(ppm|obj)> (<dst-mask.ppm> <src-mask.ppm>)? <output-basename>
    gokicheck habit   <in0.obj> <in1.obj> (<index> <max_index> <psi>)? <out.obj>
    gokicheck reshape <num_shape_per_color> <input_color.ppm> <input_shape.ppm> <output.ppm>
    python2 test.py ./gokicheck match     input0.png input1.(png|obj)
    python2 test.py ./gokicheck match0    input0.png input1.(png|obj)
    python2 test.py ./gokicheck matcho    input0.png input1.(png|obj) match.txt
    python2 test.py ./gokicheck pred      input0.png input1.png ...
    python2 test.py ./gokicheck pcopy     input0.png input1.png
    python2 test.py ./gokicheck ppred     input0.png input1.png ...
    python2 test.py ./gokicheck pextend   input0.png
    python2 test.py ./gokicheck penetrate input0.png
    python2 test.py ./gokicheck enlarge   input.png
    python2 test.py ./gokicheck obj       input.png
    python2 test.py ./gokicheck jps       input.png
    python2 test.py ./gokicheck tilt      input.png
    python2 test.py ./gokicheck btilt     input.png
    python2 test.py ./gokicheck btilt2    input.png
    python2 test.py ./gokicheck flicker   input.png
    python2 test.py ./gokicheck sbox      input.png
    python2 test.py ./gokicheck demosaic  input.png
    python2 test.py ./gokicheck extend    input.png
    python2 test.py ./gokicheck prep      input.png
    python2 test.py ./gokicheck prepsq    input.png
    python2 test.py ./gokicheck mask      input.png
    python2 test.py ./gokicheck mask0

# How to use as library (sample code).
Please refer tools.cc, and please include with namespace directive (but include guard definition should harms). They are NOT thread-safe.

# Demos
https://konbu.azurewebsites.net/ have a sample interface working demos.

# Tips
These program's sharpen is based on pseudo DFT half space plausible one, and, in the ideal case, fixed image from recursive sharpen command returns the ideal result.  
These program's collect is based on DFT differential.  
These program's bump assumes F=∞ graphics.
These program's match matches with calculated pseudo z-depth. There's a hidden parameters.

# Specification
filter generates the bumpmap that is pseudo plausible one because of one image condition and hypothesis, but this is correct if the hypothesis, if it's in the focal point, edge is better clear than other places, is correct.  
Pseudo condition is avoidable on very wide cases with multiple camera conditions or multiple pint conditions,
if it fails, the image resolution or color depth resolution lacks, or, something like different colours with each angle like mirrors, or, because of scattering or fog things.  

If we use this program for 3D model input like blender, please make mask with mask0 or another softwares, and convert them
with mask command, then, please use input-mask.obj . If we are lucky, in the blender example, large scaled input will inputted,
then G, R, S command adjust, please use modifier as mirror -> solidify -> skin -> lowpoly , single sided experimental result
will shown. And if we use another input instead of mirror, double sided experimental result will shown.
And if we're using with rig and so on, the existance of z-axis cover harms.

matchPartial default threshr is severe, so please configure before to match.

# Another downloads
* https://ja.osdn.net/projects/goki-check/
* https://www.sourceforge.net/projects/gokicheck/
* https://sites.google.com/view/bitsofcotton

