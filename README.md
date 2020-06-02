# Goki Check
These program aims to implement one of a complement to ongoing other utilities.

# How to use
Please touch Makefile for libc++ enabled.  
This program needs ascii raw ppm files to input/output.  
For speed, we need http://eigen.tuxfamily.org/ library.
We need imagemagick for normaluse, and ffmpeg for making .mp4.

# Parameters
* enlarge.hh
* * dratio  : z-axis step number, larger combines sqrt(dratio) number of z-axis steps.
* * dist    : z-axis scan max distance.
* * plen    : extend pixels.
* * lrecur  : recursive number for sharpen.
* * bumpd   : z-axis differential pixel number.
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
Searching the Internet more...

# Status
Freezed.

# Usage
    make tools
    
    gokicheck bump    <psi> <input.ppm> <output.ppm>
    gokicheck pextend <pixels> <input.ppm> <output.ppm>
    gokicheck collect <input.ppm> <output.ppm>
    gokicheck light   <n_recur> <input.ppm> <output.ppm>
    gokicheck reshape <num_of_color_depth> <input_color.ppm> <input_shape.ppm> <output.ppm>
    gokicheck obj     <gather_pixels> <ratio> <zratio> <thin> <input.ppm> <mask.ppm>? <output.obj>
    gokicheck tilt    <index> <max_index> <psi> <input.ppm> <input-bump.(ppm|obj)> <output.ppm>
    gokicheck sbox    <index> <max_index> <input.ppm> <input-bump.(ppm|obj)> <output.ppm>
    gokicheck match   <num_of_res_shown> <num_of_hidden_match> <vbox_dst> <vbox_src> <dst.ppm> <src.ppm> <dst-bump.(ppm|obj)> <src-bump.(ppm|obj)> (<dst-mask.ppm> <src-mask.ppm>)? <output-basename>
    gokicheck match0  <num_of_res_shown> <num_of_hidden_match> <vbox_dst> <vbox_src> <dst.ppm> <src.ppm> <dst-bump.(ppm|obj)> <src-bump.(ppm|obj)> (<dst-mask.ppm> <src-mask.ppm>)? <output-basename>
    gokicheck matcho  <match> <num_emph> <vbox_dst> <vbox_src> <dst.ppm> <src.ppm> <dst-bump.(ppm|obj)> <src-bump.(ppm|obj)> (<dst-mask.ppm> <src-mask.ppm>)? <output-basename>
    gokicheck pose    <vbox> <thresh> <pose.txt> <input.ppm> <input-bump.ppm>
    gokicheck poso    <vbox> <thresh> <pose.txt> <input.ppm> <input-bump.ppm> <num_of_res_shown> <output-base>
    gokicheck habit   <in0.obj> <in1.obj> (<index> <max_index> <psi>)? <out.obj>
    python test.py ./gokicheck col  input.png
    python test.py ./gokicheck penetrate input.png
    python test.py ./gokicheck bump input.png
    python test.py ./gokicheck enlp input.png
    python test.py ./gokicheck light input.png
    python test.py ./gokicheck pextend input.png
    python test.py ./gokicheck mask0 input.png
    python test.py ./gokicheck mask  input.png
    python test.py ./gokicheck obj  input.png
    python test.py ./gokicheck mtl  input.png
    python test.py ./gokicheck tilt input.png
    python test.py ./gokicheck btilt input.png
    python test.py ./gokicheck btilt2 input.png
    python test.py ./gokicheck flicker input.png
    python test.py ./gokicheck jps input.png
    python test.py ./gokicheck match input0.png input1.(png|obj)
    python test.py ./gokicheck match0 input0.png input1.(png|obj)
    python test.py ./gokicheck matcho input0.png input1.(png|obj) match
    python test.py ./gokicheck pose input.png
    python test.py ./gokicheck poso input.png
    python test.py ./gokicheck pcomp input0.png input1.png ...
    python test.py ./gokicheck pcomp0 input0.png input1.png ...

# How to use as library (sample code).
Please refer tools.cc, and please include with namespace directive (but include guard definition should harms).

# Demos
https://konbu.azurewebsites.net/ have a sample interface working demos.

# Tips
These program's light is based on pseudo DFT half space plausible one.  
These program's collect is based on DFT differential.  
These program's bump assumes F=∞ graphics.
These program's match matches with calculated pseudo z-depth, please configure in reDig class initializer. And, that assumes one of vertices is full and another is lowPoly but now, it isn't.

# Specification
filter generates the bumpmap that is pseudo plausible one because of one image condition and hypothesis, but this is correct if the hypothesis, if it's in the focal point, edge is better clear than other places, is correct. And therefore, z-axis sign from some base plane has not to be determinated as unique one.  
Pseudo condition is avoidable on very wide cases with multiple camera conditions or multiple pint conditions,
if it fails, the image resolution or color depth resolution lacks, or, something like different colours with each angle like mirrors, or, because of scattering or fog things.  

If we use this program for 3D model input like blender, please make mask with mask0 or another softwares, and convert them
with mask command, then, please use input-mask.obj . If we are lucky, in the blender example, large scaled input will inputted,
then G, R, S command adjust, please use modifier as mirror -> solidify -> skin -> lowpoly , single sided experimental result
will shown. And if we use another input instead of mirror, double sided experimental result will shown.
And if we're using with rig and so on, the existance of z-axis cover harms.

Match matches including z-axis. So around this, if we match with .obj file, we don't have accurate z-axis ratio, bugly result returns.
If we should have correct matches, please configure z-axis ratio in tools.cc via redig.hh .
matchPartial default threshr is severe, so please configure before to match.

filter's enlarge generates a little blurred result in many cases (and this is reduced in larger images).
This is because we generate the complement light image with DFT half space plausible ones by shifting frequency space intensities, then, getting recursive enlarged image with DFT.

# Another downloads
* https://ja.osdn.net/projects/goki-check/
* https://www.sourceforge.net/projects/gokicheck/
* https://sites.google.com/view/bitsofcotton

