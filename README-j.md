# Goki Check
現行のツールセットで提供される画像処理を補うことが目的のライブラリ群です。

# 使い方
Makefile を libc++ を使えるように変更してください。  
このプログラムは ascii 形式の ppm ファイルを入出力に使用します。  
また、速度を担保するには http://eigen.tuxfamily.org/ ライブラリが必要です。
通常の使用で imagemagick が、動画の作成に ffmpeg が必要です。

# 調整可能なパラメタ
* enlarge.hh
* * dratio  : z 軸を走査する際のステップ数です。sqrt(dratio) ステップがジッタ除去のために使われます。
* * dist    : z 軸方向をスキャンする際の最大の距離です。
* * offset  : bump マップを作成する際の発散を防ぐための最小値の比率です。0 より大きく、また、小さい値の方が良い値を返します。
* * plen    : EXTEND の際に補完するピクセル数です。
* * lrecur  : SHARPEN の際に再帰的に繰り返す回数です。
* * bumpd   : z 軸を計算する際に傾きを取得するピクセルの幅数です。
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
さらに検索中です。

# 状態
Freeze 中です。

# 使い方
    make gokicheck
    
    gokicheck bump    <psi> <input.ppm> <output.ppm>
    gokicheck pextend <pixels> <input.ppm> <output.ppm>
    gokicheck collect <input.ppm> <output.ppm>
    gokicheck light   <n_recur> <input.ppm> <output.ppm>
    gokicheck reshape <num_of_color_depth> <input_color.ppm> <input_shape.ppm> <output.ppm>
    gokicheck obj     <gather_pixels> <ratio> <zratio> <thin> <input.ppm> <mask.ppm>? <output.obj>
    gokicheck tilt    <index> <max_index> <psi> <input.ppm> <input-bump.(ppm|obj)> <output.ppm>
    gokicheck sbox    <index> <max_index> <input.ppm> <input-bump.(ppm|obj)> <output.ppm>
    gokicheck match   <num_of_res_shown> <num_of_hidden_match> <vbox_dst> <vbox_src> <dst.ppm> <src.ppm> <dst-bump.(ppm|obj)> <src-bump.(ppm|obj)> (<dst-mask.ppm> <src-mask.ppm>)? <output-basename>
    gokicheck match0   <num_of_res_shown> <num_of_hidden_match> <vbox_dst> <vbox_src> <dst.ppm> <src.ppm> <dst-bump.(ppm|obj)> <src-bump.(ppm|obj)> (<dst-mask.ppm> <src-mask.ppm>)? <output-basename>
    gokicheck matcho  <match> <num_emph> <vbox_dst> <vbox_src> <dst.ppm> <src.ppm> <dst-bump.(ppm|obj)> <src-bump.(ppm|obj)> (<dst-mask.ppm> <src-mask.ppm>)? <output-basename>
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

# ライブラリとしての使い方
tools.cc を参照してください。また、必要であれば namespace ブロックでスコープしてください。
ただし、高確率でインクルードガードの定義が有害です。

# デモ
https://konbu.azurewebsites.net/ にサンプルの準備中です。  

# Tips
light は DFT 時の半整数空間から擬似的にとってきています。  
collect は単純に DFT 微分の後、abs をとってきています。  
bump は F=&infin; を仮定しています。また、疑似的に傾けた画像に対して擬似的に 2 つのカメラで撮影したと仮定しています。  
match は z 軸方向まで含めて合致する部分を探します。reDig クラス内部で深度の比率を調節してください。また、片方が稠密な頂点、もう片方が lowPoly された頂点になっている入力を仮定していますが、現状そうなってはいません。

# 仕様
filter2ex はもっともらしいバンプマップを返しますが、正しくない場合があります。
これは、1 枚の画像のみを使用する事と仮定している構造である、もし焦点が合っていればより角が立ってみえる、という構造によるもので、これによると、ある基準の面に対しての z 軸の符号は一意には定まりません。  
しかし、ほとんどの場合では複数台のカメラあるいは複数のピントを使えば正しいバンプマップを返すことができます。
(それができない場合には画像や色の解像度が足りないか、鏡など角度によって異なる色を返すものや、散乱やフォグなどの光学的事象がある際です)

また、blender など 3D 編集ソフトへの入力として使用する際には mask0 コマンドでマスクを作ったあともしくは、他でマスクを作り、
mask コマンドで変換後、filename-mask.obj ファイルを入力に使用してください。
うまくいくと、blender の場合には scale がとても大きい結果として入力され、G, R, S コマンドで調整したあとに、
モディファイアで mirror -> solidify -> skin -> lowpoly などの経路できちんとした片面の結果が得られる可能性があります。
mirror を裏表にすることにより、もしかするときちんとした両面の結果が得られます。
また、Rig などとともに使用する際には、z 軸方向に重なりのない入力を使用しないとおかしな結果になります。

match コマンドは深度情報を含めて合致します。その際に、.obj データなどは若干おかしな結果が返ることがあります。
より正確に合致するには深度情報の比率を適宜変えて合致してみてください。
また、デフォルトでは threshr の値がシビアなので、適宜調整してください。

filter2ex による拡大は大きな画像に対しては比較的安定な画像を返しますが、現行使用されている拡大のアルゴリズムよりは劣ります。
これは半整数空間の周波数成分を参照しながら、色深度を補正しつつ、DFT 拡大を再起的に行うことによって拡大を得ているためです。

# その他のダウンロードサイト
* https://ja.osdn.net/projects/goki-check/
* https://www.sourceforge.net/projects/gokicheck/
* https://sites.google.com/view/bitsofcotton

# 検討中のもの
メタボールの中心となるような線分と重み付けの森を 3D 模型から計算することが簡単にできるか散策しています。

