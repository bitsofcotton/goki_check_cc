# Goki Check
静止画像からあらかじめ想定される模型を座標と共に得ることが目的のプログラムです。このプログラムは決定論的な手法を用いています。  
古い情報に関しては、https://sourceforge.net/p/gokicheck/wiki/Home/ を参照してください。

# 使い方
Makefile を stdc++ を使えるように変更してください。  
このプログラムは ascii 形式の ppm ファイルを入出力に使用します。  
変換の際は、例えばパワルフなツールである https://www.imagemagick.org/ の 'convert from.image -compress none to.ppm' が使用できます。  
また、速度を担保するには http://eigen.tuxfamily.org/ ライブラリが、bone 情報を担保した合致をする際には https://github.com/jessey-git/fx-gltf/ ライブラリが必要です。

# 調整可能なパラメタ
* redig.hh
* * vbox : ベクタ生成の際にまとめるピクセルの数です。
* * rz   : 奥行きの乗数です。
* scancontext.hh
* * matchPartialPartial::ndiv    : 合致する角度の分割数です。合致の誤差にも影響します。
* * matchPartialPartial::threshp : 検出される最小の合致する点の数の比率です。
* * matchPartialPartial::threshs : 合致が同じかどうか判定する閾値の比率です。

# 文脈
写真の後でのピント調整プログラムに刺激されました。また、この分野に関して、様々な(これと異なる)付帯条件での先行がたくさんありました。
(例えば、たくさんのカメラによる画像を使用するものや、あらかじめレイヤ毎に用意しておくもの、球状に膨らませるもの、動画から生成するものなどです。)
検索結果に defocus photo アルゴリズムがありました。 http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.100.2308&rep=rep1&type=pdf 。これはほとんどの場合のカメラ撮影について、goki_check_cc よりも正確です。goki_check_cc は異なる仮定を用いていますが、一般的に使用されるカメラ撮影の場合には、事前のいくつかの変換が必要です。  
また、合致の分野に対して、様々な(これと異なる)付帯条件での先行がたくさんありました。(例えば、トポロジを検出して端を合致するものや、機械学習を使うもの、および座標変換の変数を点の数で固定するものなどです。)。また、特徴点を利用した PnP 問題の合致の方がこちらよりも高速です。  
検索によって限られた組み合わせの語句では見つからなかった https://ryo620.org/2018/02/to-gltf-from-fbx-by-blender/ という記事さんがありました。gltf2 のライブラリは初めて知りました。  
さらに検索中です。

# 状態
Freeze 前の細かな実装のチェックをしています。

# 使い方
    make tools
    
    # enlarge
    ./tools enlarge input.ppm output.ppm
    
    # detect edges with absolute differential with dft value.
    ./tools collect input.ppm output.ppm
    
    # make 2d to 3d pseudo bumpmap
    ./tools bump input.ppm output.ppm
    
    # make 2d to 3d persistent pseudo bumpmap
    ./tools pbump input.ppm output.ppm
    
    # bumpmap to .obj file
    ./tools obj input-bump.ppm output.obj
    
    # mask .obj files and add support.
    ./tools maskobj input-mask.ppm input.obj output.obj size-ratio thin-or-thickness-size
    
    # make tilts from original and bumpmap images.
    ./tools tilt input.ppm output-base input-bump.ppm
    
    # list matches.
    ./tools match input-matchbase.ppm output-base input-tobematched.ppm matchbase-bump.ppm tobematched-bump.ppm input-mask.ppm tobematched-mask.ppm
    
    # list matches 2d - 3d.
    ./tools match3d input-matchbase.ppm output-base input-tobematched.obj matchbase-bump.ppm matchbase-mask.ppm
    
    # list matches 2d - 2d with hidden 3d.
    ./tools match2dh3d input-matchbase.ppm output-base input-tobematched.ppm bump-matchbase.ppm bump-tobematched.ppm mask-matchbase.ppm mask-tobematched.ppm hidden-object.obj
    
    # list matches 2d - 3d.
    ./tools match3dbone input-matchbase.ppm output-base input-tobematched.gltf matchbase-bump.ppm matchbase-mask.ppm
    
    # list matches 2d - 2d with hidden 3d.
    ./tools match2dh3dbone input-matchbase.ppm output-base input-tobematched.ppm bump-matchbase.ppm bump-tobematched.ppm mask-matchbase.ppm mask-tobematched.ppm hidden-object.gltf

    # habit
    ./tools habit mask.ppm output.obj input0.obj input1.obj
    
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
ほとんどの場合では複数台のカメラを使えば正しいバンプマップを返すことができます。
(それができない場合には画像や色の解像度が足りないか、鏡など角度によって異なる色を返すものや、散乱やフォグなどの光学的事象がある際です)

また、このプログラムは大域的には平坦なバンプマップを返します。(これはスケール別に適用し直すと恐らく補正できます pbump コマンドを試してみてください)。
これにより、合致の際に使用される模型が薄くなりすぎて合致されない可能性があります。

# その他のダウンロードサイト
* https://ja.osdn.net/projects/goki-check/
* https://www.sourceforge.net/projects/gokicheck/
* https://files.limpid-intensity.info/

# あらかじめ計算されたサンプル
![photosample pseudo bumpmap](https://files.limpid-intensity.info/photosample-bump.jpeg)
