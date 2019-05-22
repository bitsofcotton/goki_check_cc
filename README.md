# Goki Check
These program aims to get prepared model geometry in A still image in a deterministic way.  
Please refer older information at https://sourceforge.net/p/gokicheck/wiki/Home/ .

# How to use
Please touch Makefile for stdc++ enabled.  
This program needs ascii raw ppm files to input/output.  
To convert image files to raw ppm, it is powerful tool that https://www.imagemagick.org/ with 'convert from.image -compress none to.ppm'.   
For speed, we need http://eigen.tuxfamily.org/ library.

# Parameters
* enlarge.hh
* * dratio : z-axis step ratio relate to 1., must be &lt; 1.
* * offset : color value ratio around convergences in making bump map.
* redig.hh
* * vbox : size of vector gathering rectangle.
* * rz   : z-axis output ratio.
* match.hh
* * matchPartial::ndiv    : number of divides that match angles, effects the matching errors.
* * matchPartial::threshr : error tolerance rate for matching, smaller is tight.
* * matchPartial::threshp : ratio of threshold for number of matched points.
* * matchPartial::threshs : ratio of threshold for operator ==.

# Context
This program is inspired from re-focus photo softwares.  
And around this, there's many preceders that many approach to get bump maps with certain conditions
(such as multiple camera conditions, or, with layered objects, or, spherical, or, from movie, etcetc).
There's a defocus photo algorithms http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.100.2308&rep=rep1&type=pdf some words with googled. So it's accurate for most cameras, goki_check_cc is standing on another hypothesis that is not widely used (in another words, some transform or special camera is needed for photos...). And, there exists preceders that make it from multiple pint images, this makes very accurate results.  
There's preceders to match 3D to 2D with many approaches. (s.t. detecting topology of junction point, or, machine learning, and so on.). And it is fater than this that PnP problem and specific point based matching.  
By searching with some word that is not common, there exists the article https://ryo620.org/2018/02/to-gltf-from-fbx-by-blender/ that I firstly know the gltf format by this. There's a https://github.com/jessey-git/fx-gltf/ library, but compatibility for this is abandoned.  
Searching the Internet more...

# Status
Checking details of implementation before to freeze the library.  

# Usage
    make tools
    
    gokicheck enlarge <ratio>  <input.ppm> <output.ppm>
    gokicheck pextend <pixels> <input.ppm> <output.ppm>
    gokicheck collect <input.ppm> <output.ppm>
    gokicheck idetect <input.ppm> <output.ppm>
    gokicheck bump    <input.ppm> <output.ppm>
    gokicheck obj     <shift_x_pixels> <gather_pixels> <zratio> <input.ppm> <mask.ppm>? <output.obj>
    gokicheck obj     stand <gather_pixels> <thin> <ratio> <zratio> <input.ppm> <mask.ppm>? <output.obj>
    gokicheck tilt    <index> <max_index> <psi> <shift_x_pixels> <input.ppm> <input-bump.(ppm|obj)> <output.ppm>
    gokicheck draw    <input-mask.ppm> <input-obj.(obj|gltf)> <output.ppm>
    gokicheck match   <num_of_res_shown> <num_of_hidden_match> <vbox_dst> <vbox_src> <dst.ppm> <src.ppm> <dst-bump.(ppm|obj)> <src-bump.(ppm|obj|gltf)> (<dst-mask.ppm> <src-mask.ppm>)? <output-basename>
    gokicheck matcho  <match> <num_emph> <vbox_dst> <vbox_src> <dst.ppm> <src.ppm> <dst-bump.(ppm|obj)> <src-bump.(ppm|obj|gltf)> (<dst-mask.ppm> <src-mask.ppm>)? <output-basename>
    gokicheck habit   <in0.obj> <in1.obj> (<index> <max_index> <psi>)? <out.obj>
    python test.py ./gokicheck col  input.png
    python test.py ./gokicheck bump input.png
    python test.py ./gokicheck emph input.png
    python test.py ./gokicheck emphe input.png
    python test.py ./gokicheck obj  input.png
    python test.py ./gokicheck mtl  input.png
    python test.py ./gokicheck scn  input.png
    python test.py ./gokicheck tilt input.png
    python test.py ./gokicheck btilt input.png
    python test.py ./gokicheck flicker input.png
    python test.py ./gokicheck pnga input.png
    python test.py ./gokicheck jps input.png
    python test.py ./gokicheck match input0.png input1.(png|obj)
    python test.py ./gokicheck matcho input0.png input1.(png|obj) match

# How to use as library (sample code).
Please refer tools.cc, and please include with namespace directive (but include guard definition should harms).

# Demos
https://services.limpid-intensity.info/ have a sample interface working demos.
Please bookmark output directory page after sending images, then please wait some minutes.  
N.B. a file per 5 minutes, up to 20 Mo total upload size.

# Tips
These program's enlarge is based on pseudo DFT half space plausible one.  
These program's collect is based on DFT differential.  
These program's bump assumes F=∞ graphics.  
These program's match matches with calculated pseudo z-depth, please configure in reDig class initializer.  
These program's match assumes one of vertices is full and another is lowPoly but now, it isn't.

# Specification
enlarger2ex generates the bumpmap that is pseudo plausible one because of one image condition and hypothesis, but this is correct if the hypothesis, if it's in the focal point, edge is better clear than other places, is correct.  
Pseudo condition is avoidable on very wide cases with multiple camera conditions or multiple pint conditions,
if it fails, the image resolution or color depth resolution lacks, or, something like different colours with each angle like mirrors, or, because of scattering or fog things.
And, generated bump map is NOT correct in global, so please change local-global ratio in enlarge.hh:enlarger2ex::compute::BUMP_Y.  
And, there's surely small chances that enhances bump command on emphe command processed image.

If we use this program for 3D model input like blender, please make mask with mask0 or another softwares, and convert them
with mask command, then, please use input-mask.obj . If we are lucky, in the blender example, large scaled input will inputted,
then G, R, S command adjust, please use modifier as mirror -> solidify -> skin -> remesh , single sided experimental result
will shown. And if we use another input instead of mirror, double sided experimental result will shown.
And if we're using with rig and so on, the existance of z-axis cover harms.

And if we use this program enlarge for vectorize software input like inkscape, region will gain better result, but,
colour itself worse results. So first calculate region, then, please use average of the region for colors.

Match matches including z-axis. So around this, if we match with .obj file, we don't have accurate z-axis ratio,
bugly result returns. If we know angle with projected plane, please tilt with objtilt command and make one-sided result,
match matches. If we should have correct matches, please configure z-axis ratio in tools.cc via redig.hh .
matchPartial default threshr is severe, so please configure before to match.

enlarger2ex 's enlarge generates a little blurred result in many cases (and this is reduced in larger images).
This is because we generate the one with DFT half space plausible ones by shifting frequency space intensities.

# Another downloads
* https://ja.osdn.net/projects/goki-check/
* https://www.sourceforge.net/projects/gokicheck/
* https://files.limpid-intensity.info/

# Pre calculated samples
![photosample pseudo bumpmap](https://files.limpid-intensity.info/photosample-bump.jpeg)

# Being examined
It is searching that can we convert 3d model into the forests that is the center lines of metaball and depth simply.
