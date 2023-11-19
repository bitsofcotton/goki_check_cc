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
From some news, there exists the normal vector of the light based method that is making 3D model from single picture condition with some machine learning.  
Searching the Internet more...

# Usage
    make gokibin
    
    gokibin (collect|sharpen|bump|enlarge|flarge|blink|represent|nop|limit|bit) <input.ppm> <output.ppm> <recursive_num> <rot_num>
    gokibin (cat|catr) <input0.ppm> ...
    gokibin (tilt|sbox) <index> <max_index> <psi> <input.ppm> <input-bump.ppm> <output.ppm>
    gokibin obj <input.ppm> <output.obj>
    gokibin match <num_of_match> <num_of_emph> <vbox_dst> <vbox_src> <dst.ppm> <src.ppm> <dst-bump.ppm> <src-bump.ppm> <output-basename>
    gokibin recolor  <dimension> <input.ppm> <input-copy.ppm> <output.ppm> <intensity>
    gokibin recolor2 <dimension> <input.ppm> <output.ppm> <intensity>
    gokibin recolor3 <dimension> <input.ppm> <input-shape> <output.ppm>
    gokibin habit <in0.obj> <in1.obj> <out.obj>
    python3 test.py ./gokibin (sharpen|bump|enlarge|flarge|represent|jps|tilt|obj|sbox|prep|presq|nop|limit|bit|cleansq) <param> input0.png ...
    python3 test.py ./gokibin (cat|catb|catr|catbr) input0.png input1.png ...
    python3 test.py ./gokibin (tilecat|tilecatb|tilecatr|tilecatbr) <tile count> < cat.txt
    python3 test.py ./gokibin match input0.png input1.png <vboxdst> <vboxsrc> <number_of_submatches> <number_of_emphasis>
    python3 test.py ./gokibin i2i <param> img0.png ...
    python3 test.py qredg pextend in.png
    python3 test.py predg pred in0.png ...

# Another downloads
* https://konbu.azurewebsites.net/ (Sample Site)
* https://drive.google.com/drive/folders/1B71X1BMttL6yyi76REeOTNRrpopO8EAR?usp=sharing
* https://1drv.ms/u/s!AnqkwcwMjB_PaDIfXya_M3-aLXw?e=qzfKcU
* https://ja.osdn.net/users/bitsofcotton/
* https://www.sourceforge.net/projects/gokicheck/ (abandoned)

# Real close
2023/03/03
2023/03/13 integrate some files into lieonn.hh after close #1.
2023/03/19 add per depth predbit command, after close #2.
2023/03/20 elim predbit, add bit cmd, after close #3.
2023/03/24 code clean, after close #4.
2023/04/02 merge catg fix.
2023/04/03 merge.
2023/04/05 improve accuracy stability.
2023/04/20 obj+ command.
2023/04/29 theoretical fix and their close on bump, afterbump, gettilevec chain.
2023/05/07 add functions around cvstereo.py.
2023/06/24 fix around z-axis scale on afterbump, get...vec, obj, tilt, sbox commands.
2023/06/29 fix same place z-axis with reasonable.
2023/07/08 integrate gokicvs to test.py.
2023/09/04 large fix around bump/obj/jps, better reasonable on logic.
2023/09/06 bump local/global algorithm strategy change.
2023/10/05 update lieonn to latest one.
2023/10/06 curvature fix, also don't need to local to global trans in bump_both.
2023/10/07 test.py bump command change, add clean command. update readme.
2023/10/18 update around bump, we had should to cleansq after bump. really close around bump with this, this is the reasonable one for one picture condition. code clean.
2023/10/20 enlarge cmd improve and fix.
2023/10/24 pextend command retry.
2023/10/27 pextend command close. prepare to all close.
2023/10/30 copy structure reliably with randtools meaning.
2023/11/19 add pred command.

