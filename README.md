# Goki Check
These program aims to implement one of a complement to ongoing other utilities.  
And this library is written in a deterministic way.

# How to use
Please touch Makefile for libc++ enabled.  
This program needs ascii raw ppm files to input/output.  
We can use /dev/stdin, /dev/stdout to in/output.  
We need imagemagick for normal use, and ffmpeg to make .mp4.  

# Context
This program is inspired from re-focus photo softwares.  
And around this, there's many preceders that many approach to get bump maps with certain conditions
(such as multiple camera conditions, or, with layered objects, or, spherical, or, from movie, etcetc).
There's a defocus photo algorithms http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.100.2308&rep=rep1&type=pdf some words with googled. So it's accurate for most cameras, goki_check_cc is standing on another hypothesis that is not widely used (in another words, some transform or special camera is needed for photos...). And, there exists preceders that make it from multiple pint images, this makes very accurate results.  
There's preceders to match 3D to 2D with many approaches. (s.t. detecting topology of junction point, or, machine learning, and so on.). And it is fater than this that PnP problem and specific point based matching.  
By searching with some word that is not common, there exists the article https://ryo620.org/2018/02/to-gltf-from-fbx-by-blender/ that I firstly know the gltf format by this. There's a https://github.com/jessey-git/fx-gltf/ library, but compatibility for this is abandoned.  
From some news, there exists the machine learning method that is making 3D model from single picture condition.
Searching the Internet more...

# Usage
    make gokicheck
    
    gokicheck (collect|sharpen|bump|enlarge|flarge|pextend|blink|represent) <input.ppm> <output.ppm> <recursive_num> <rot_num>
    gokicheck (pred|composite|lenl)  <output.ppm> <input0.ppm> ...
    gokicheck (cat|catr) <input0.ppm> ...
    gokicheck obj   <gather_pixels> <ratio> <input.ppm> <output.obj>
    gokicheck (tilt|sbox) <index> <max_index> <psi> <input.ppm> <input-bump.ppm> <output.ppm>
    gokicheck match <num_of_match> <num_of_emph> <vbox_dst> <vbox_src> <dst.ppm> <src.ppm> <dst-bump.ppm> <src-bump.ppm> <output-basename>
    gokicheck recolor  <dimension> <input.ppm> <input-copy.ppm> <output.ppm> <intensity>
    gokicheck recolor2 <dimension> <input.ppm> <output.ppm> <intensity>
    gokicheck recolor3 <dimension> <input.ppm> <input-shape> <output.ppm>
    gokicheck retrace  <dimension> <input-mask.ppm> <input-mask-copy.ppm> <intensity>
    gokicheck retrace2 <dimension> <input-mask.ppm> <outbase> <intensity>
    gokicheck newtrace <dimension> <size> <output.ppm>
    gokicheck reimage  <dimension> <input.ppm> <input-src.ppm> <output.ppm> <intensity>
    gokicheck reimage2 <dimension> <input.ppm> <output.ppm> <intensity>
    python2 test.py ./gokicheck match input0.png input1.png
    python2 test.py ./gokicheck (tilecat|tilecatb|tilecatr|tilecatbr) <tile count> < cat.txt
    python2 test.py ./gokicheck (cat|catb|catr|catbr|pred|composite|lenl) input0.png input1.png ...
    python2 test.py ./gokicheck newtrace dimension size
    python2 test.py ./gokicheck retrace  dimension input.ppm intensity
    python2 test.py ./gokicheck retrace2 dimension dst.ppm src.ppm intensity
    python2 test.py ./gokicheck (pextend|sharpen|bump|penlarge|penl|flarge|represent|jps|tilt|btilt|obj|sbox|prep|presq) input.png

# How to use as library (sample code).
Please refer tools.cc, and please include with namespace directive
(but include guard definition should harms).  

# Known Bug not to be fixed.
reDig::edge returns relaxed edges, this can be fixed with multiple on counting
+/- on same scan line.

reDig::mesh2 returns something buggy when we input non lattice points.

The bump command returns tilt ignored result.
This is because we integrate local focuses and integration itself ignores differential offset part.

# Another downloads
* https://konbu.azurewebsites.net/ (Sample Site)
* https://drive.google.com/drive/folders/1B71X1BMttL6yyi76REeOTNRrpopO8EAR?usp=sharing
* https://1drv.ms/u/s!AnqkwcwMjB_PaDIfXya_M3-aLXw?e=qzfKcU
* https://ja.osdn.net/users/bitsofcotton/
* https://www.sourceforge.net/projects/gokicheck/ (abandoned)

