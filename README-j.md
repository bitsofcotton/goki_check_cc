# Goki Check
静止画像からあらかじめ想定される模型を座標と共に得ることが目的のプログラムです。このプログラムは決定論的な手法を用いています。  
古い情報に関しては、https://sourceforge.net/p/gokicheck/wiki/Home/ を参照してください。

# 使い方
Makefile を stdc++ を使えるように変更してください。  
このプログラムは ascii 形式の ppm ファイルを入出力に使用します。  
変換の際は、例えばパワルフなツールである https://www.imagemagick.org/ の 'convert from.image -compress none to.ppm' が使用できます。  
また、速度を担保するには http://eigen.tuxfamily.org/ ライブラリが、bone 情報を担保した合致をする際には https://github.com/jessey-git/fx-gltf/ ライブラリが必要です。

# 調整可能なパラメタ
* enlarge.hh
* * dratio : z 軸を走査する際のステップ幅です。1 に対して記述し、1 より小さくなくてはいけません。
* * offset : bump マップを作成する際の発散を防ぐための最小値の比率です。0 より大きく、また、小さい値の方が良い値を返します。
* redig.hh
* * vbox : ベクタ生成の際にまとめるピクセルの数です。
* * rz   : 奥行きの乗数です。
* match.hh
* * matchPartialPartial::ndiv    : 合致する角度の分割数です。合致の誤差にも影響します。
* * matchPartialPartial::threshr : 合致に許容する誤差の比率です。
* * matchPartialPartial::threshp : 検出される最小の合致する点の数の比率です。
* * matchPartialPartial::threshs : 合致が同じかどうか判定する閾値の比率です。

# 文脈
写真の後でのピント調整プログラムに刺激されました。また、この分野に関して、様々な(これと異なる)付帯条件での先行がたくさんありました。
(例えば、たくさんのカメラによる画像を使用するものや、あらかじめレイヤ毎に用意しておくもの、球状に膨らませるもの、動画から生成するものなどです。)
検索結果に defocus photo アルゴリズムがありました。 http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.100.2308&rep=rep1&type=pdf 。これはほとんどの場合のカメラ撮影について、goki_check_cc よりも正確です。goki_check_cc は異なる仮定を用いていますが、一般的に使用されるカメラ撮影の場合には、事前のいくつかの変換が必要です。また、フォーカスをずらして 2 枚以上の写真を撮影できる場合には、より正確な先行結果があります。  
また、合致の分野に対して、様々な(これと異なる)付帯条件での先行がたくさんありました。(例えば、トポロジを検出して端を合致するものや、機械学習を使うもの、および座標変換の変数を点の数で固定するものなどです。)。また、特徴点を利用した PnP 問題の合致の方がこちらよりも高速です。  
検索によって限られた組み合わせの語句では見つからなかった https://ryo620.org/2018/02/to-gltf-from-fbx-by-blender/ という記事さんがありました。gltf2 のライブラリは初めて知りました。  
さらに検索中です。

# 状態
Freeze 前の細かな実装のチェックをしています。特に、gltf2 ライブラリへの対応はまだバグが潜在しています。また、合致周りは beta です。

# 使い方
    make gokicheck
    
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
    python test.py ./gokicheck bump input.png
    python test.py ./gokicheck obj  input.png
    python test.py ./gokicheck mtl  input.png
    python test.py ./gokicheck scn  input.png
    python test.py ./gokicheck tilt input.png
    python test.py ./gokicheck btilt input.png
    python test.py ./gokicheck flicker input.png
    python test.py ./gokicheck pnga input.png
    python test.py ./gokicheck jps input.png
    python test.py ./gokicheck match input0.png input1.(png|obj|gltf)
    python test.py ./gokicheck matcho input0.png input1.(png|obj|gltf) match

# ライブラリとしての使い方
tools.cc を参照してください。また、必要であれば namespace ブロックでスコープしてください。
ただし、高確率でインクルードガードの定義が有害です。

# デモ
https://services.limpid-intensity.info/ にサンプルがあります。  
画像をアップロードした後、リンク先のページをブックマークしてください。その後、計算終了までに空いていて数分かかります。  
N.B. 5分に一回バッチが回ります。アップロードのサイズを 20Mo までに制限してあります。

# Tips
enlarge は DFT 時の半整数空間から擬似的にとってきています。  
collect は単純に DFT 微分の後、abs をとってきています。  
bump は F=&infin; を仮定しています。  
match は z 軸方向まで含めて合致する部分を探します。reDig クラス内部で深度の比率を調節してください。  
match は片方が稠密な頂点、もう片方が lowPoly された頂点で有る入力を仮定していますが、現在そうなってはいません。

# 仕様
enlarger2ex はもっともらしいバンプマップを返しますが、正しくない場合があります。
これは、1 枚の画像のみを使用する事と仮定している構造である、もし焦点が合っていればより角が立ってみえる、という構造によるもので、
ほとんどの場合では複数台のカメラあるいは複数のピントを使えば正しいバンプマップを返すことができます。
(それができない場合には画像や色の解像度が足りないか、鏡など角度によって異なる色を返すものや、散乱やフォグなどの光学的事象がある際です)
また、このプログラムは大域的には正しくないバンプマップを返します。補正するには enlarger2ex 内部で局所大域の比率を調整してください。

また、blender など 3D 編集ソフトへの入力として使用する際には mask0 コマンドでマスクを作ったあともしくは、他でマスクを作り、
mask コマンドで変換後、filename-mask.obj ファイルを入力に使用してください。
うまくいくと、blender の場合には scale がとても大きい結果として入力され、G, R, S コマンドで調整したあとに、
モディファイアで mirror -> solidify -> skin -> remesh などの経路できちんとした片面の結果が得られる可能性があります。
mirror を裏表にすることにより、もしかするときちんとした両面の結果が得られます。
また、Rig などとともに使用する際には、z 軸方向に重なりのない入力を使用しないとおかしな結果になります。

match コマンドは深度情報を含めて合致します。その際に、.obj データなどは若干おかしな結果が買えることがあります。
その場合、objtilt コマンドで射影を生成してから使うとうまくいく場合がありますが、この場合角度がある程度見当がついている
場合に限ります。より正確に合致するには深度情報の比率を適宜変えて合致してみてください。
また、デフォルトでは threshr の値がシビアなので、適宜調整してください。

enlarger2ex による拡大は大きな画像に対しては比較的安定な画像を返しますが、現行使用されている拡大のアルゴリズムよりは劣ります。
これは半整数空間の周波数成分を参照しながら、周波数空間をずらすことによって拡大を得ているためです。

# その他のダウンロードサイト
* https://ja.osdn.net/projects/goki-check/
* https://www.sourceforge.net/projects/gokicheck/
* https://files.limpid-intensity.info/

# あらかじめ計算されたサンプル
![photosample pseudo bumpmap](https://files.limpid-intensity.info/photosample-bump.jpeg)
