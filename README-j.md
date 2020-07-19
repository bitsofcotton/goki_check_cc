# Goki Check
現行のツールセットで提供される画像処理を補うことが目的のライブラリ群です。

# 使い方
Makefile を libc++ を使えるように変更してください。  
このプログラムは ascii 形式の ppm ファイルを入出力に使用します。  
また、速度を担保するには http://eigen.tuxfamily.org/ ライブラリが必要です。  
通常の使用で imagemagick が、動画の作成に ffmpeg が必要です。  
また、このプログラムは入力にある程度連続な画像を仮定します。つまり、細かい成分がたくさんある画像に関してはあまり良い結果を返しません。また、pextend, enlarge コマンドを繰り返し使用することも同様の理由で良い結果を返しません。

# 調整可能なパラメタ
* enlarge.hh
* * recur : filter 適用の擬似繰り返し回数です。
* redig.hh
* * vbox : ベクタ生成の際にまとめるピクセルの数です。
* * rz   : 奥行きの乗数です。
* match.hh
* * ndiv    : 合致する角度の分割数です。合致の誤差にも影響します。
* * threshr : 合致に許容する誤差の比率です。
* * threshp : 検出される最小の合致する点の数の比率です。
* * threshs : 合致が同じかどうか判定する閾値の比率です。

# 文脈
写真の後でのピント調整プログラムに刺激されました。また、この分野に関して、様々な(これと異なる)付帯条件での先行がたくさんありました。
(例えば、たくさんのカメラによる画像を使用するものや、あらかじめレイヤ毎に用意しておくもの、球状に膨らませるもの、動画から生成するものなどです。)
検索結果に defocus photo アルゴリズムがありました。 http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.100.2308&rep=rep1&type=pdf 。これはほとんどの場合のカメラ撮影について、goki_check_cc よりも正確です。goki_check_cc は異なる仮定を用いていますが、一般的に使用されるカメラ撮影の場合には、事前のいくつかの変換が必要です。また、フォーカスをずらして 2 枚以上の写真を撮影できる場合には、より正確な先行結果があります。  
また、合致の分野に対して、様々な(これと異なる)付帯条件での先行がたくさんありました。(例えば、トポロジを検出して端を合致するものや、機械学習を使うもの、および座標変換の変数を点の数で固定するものなどです。)。また、特徴点を利用した PnP 問題の合致の方がこちらよりも高速です。  
検索によって限られた組み合わせの語句では見つからなかった https://ryo620.org/2018/02/to-gltf-from-fbx-by-blender/ という記事さんがありました。gltf2 のライブラリは初めて知りました。https://github.com/jessey-git/fx-gltf/ さんへの対応を検討していましたが、見送られました。  
ニュースサイトさんによって機械学習を使用した1枚の画像からの 3D モデルの算出メソッドがあることを知りました。側面と背面の情報を学習したモデルによって補完するようです。
さらに検索中です。

# 状態
バグ情報を受け付けています。

# 使い方
    make gokicheck
    
    gokicheck (collect|sharpen|bump|enlarge|pextend) <input.ppm> <output.ppm> <recursive_num> <rotate_num>
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
    python2 test.py ./gokicheck ppred     input0.png input1.png ...
    python2 test.py ./gokicheck pextend   input.png
    python2 test.py ./gokicheck penetrate input.png
    python2 test.py ./gokicheck sharpen   input.png
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
    python2 test.py ./gokicheck mask0     input.png

# ライブラリとしての使い方
tools.cc を参照してください。また、必要であれば namespace ブロックでスコープしてください。
ただし、高確率でインクルードガードの定義が有害です。  
また、これらのライブラリはスレッドセーフではありません。

# デモ
https://konbu.azurewebsites.net/ にサンプルがあります。

# Tips
デフォルトでは threshr の値がシビアです。

# その他のダウンロードサイト
* https://ja.osdn.net/projects/goki-check/
* https://www.sourceforge.net/projects/gokicheck/
* https://sites.google.com/view/bitsofcotton

