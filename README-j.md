# Goki Check
静止画像からあらかじめ想定される文脈を得ることが目的のプログラムです。  
その後、得た文脈から他のモデルを使って再描画、あるいは、得た文脈のモデルを辞書に追加することが計画されています。  
現在はボーン情報のある固いオブジェクト同士の合致を実装しています。

# 使い方
Makefile を Eigen と stdc++ を使えるように変更してください。  
このプログラムは ascii 形式の ppm ファイルを入出力に使用します。  
変換の際は、例えばパワルフなツールである imagemagick の 'convert from.image -compress none to.ppm' が使用できます。

# 調整可能なパラメタ
* fisheye.hh
* * z_max   : 出力する z 軸の解像度です。
* * stp     : ぼやけ具合を検出する際に使用される点の数です。
* * crowd   : 混みすぎる点を緩和する正方形のサイズです。
* * vmax    : pseudoBumpVec で返される点の最大数です。(若干上下します)
* * rdist   : カメラと距離を置いて計算される平面との距離の比です
* * nloop   : 傾けながら z 軸を補正していく際のループ回数です
* * ndiv    : 傾ける際の最小の角度の分割数です
* tilt.hh
* * z_atio : [0,1] から [0,z_atio] への線形写像。
* scancontext.hh
* * matchPartialPartial::thresh  : 平行なベクトルのための誤差です。1 - &epsilon;
* * matchPartialPartial::threshl : 平行なベクトルのノルムのための誤差です。 1 - &epsilon;
* * matchPartialPartial::threshN : 平行なベクトルを検出する際の角度のノルムの閾値です。
* * matchPartialPartial::threshp : 検出される最小の合致する点の数です。
* * matchPartialPartial::threshr : 検出される最小の画像倍率です。

# 文脈
写真の後でのピント調整プログラムに刺激されました。
また、この分野に関して、様々な(これと異なる)付帯条件での先行がたくさんありました。
(例えば、たくさんのカメラによる画像を使用するものや、あらかじめレイヤ毎に用意しておくもの、球状に膨らませるものなど。
 ただし、正しい公式ページがどれなのかわからないために、リンクは差し控えさせていただきます。)

# 状態
回転のある場合の全体の一致の実装をしています。

# 使い方
    make tools
    
    # enlarge
    ./tools enlarge input.ppm output.ppm
    
    # enlarge for left differential
    ./tools enlargeds input.ppm output.ppm
    
    # detect edges with absolute differential with dft value.
    ./tools collect input.ppm output.ppm
    
    # make 2d to 3d pseudo bumpmap
    ./tools bump input.ppm output.ppm output-delaunay.ppm output.obj
    
    # make tilts from original and bumpmap images.
    ./tools tilt input.ppm output-base input-bump.ppm
    
    # make (pseudo) lowpoly and get match.ppm and .obj file.
    ./tools lpoly input.ppm output-match.ppm output.obj
    
    # list matches.
    ./tools match input-matchbase.ppm output-base input-to-bematched.ppm ref-to-be-matched.ppm ref-matchbase.ppm

# ライブラリとしての使い方
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> input;
    
    #include "enlarge.hh"
    enlarger2ex<float> enlarger;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> enlarged(enlarger.enlarge2(input, enlarger2ex<float>::ENLARGE_QUAD));
    
    enlarger2exds<float> enlargerds;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> enlargedds(enlargerds.enlarge2ds(input, enlarger2exds<float>::ENLARGE_BOTH));
    
    #include "edgedetect.hh"
    edgedetect<float> detect;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> edgecollect(detect.detect(input, edgedetect<float>::COLLECT_QUAD));
    
    #include "fisheye.hh"
    PseudoBump<float> bump;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> bumpped(bump.getPseudoBump(input));
    
    #include "tilt.hh"
    tilter<float> tilt;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> tilted(tilt.tilt(input, bumpped, 0, 8, .95));
    
    #include "scancontext.hh"
    lowFreq<float> lf;
    matchPartialPartial<float> statmatch;
    std::vector<Eigen::Matrix<float, 3, 1> > shape(lf.getLowFreq(input, 300));
    std::vector<Eigen::Matrix<float, 3, 1> > shape2(lf.getLowFreq(input2, 300));
    statmatch.init(shape0, .85, .25, .125, .3);
    std::vector<match_t<float> > matches(statmatch.match(shape1, 20));
    // match operations.
    
    // If you need, please scope with namespace block.
    // but include guard definition may harms.

# デモ
http://services.limpid-intensity.info/ にサンプルがあります。  
画像をアップロードした後、リンク先のページをブックマークしてください。その後、計算終了までに空いていて数分かかります。  
N.B. 5分に一回バッチが回ります。アップロード順ではなく sha256 の名前順です。  
N.B. アップロードのサイズを 20Mo までに制限してあります。

# Tips
enlarge と enlarge-ds は DFT 時の半整数空間から擬似的にとってきています。  
collect は単純に DFT 微分の後、abs をとってきています。  
bump は F=&infin; と y軸方向への偏光フィルタを仮定しています。  
match は入力されるファイルがバンプマップであることを仮定しています。また、内部で深度の比率を調節してください。  
match は z 軸方向まで含めて合致する部分を探します。  
match3d は入力にバンプマップと .obj ファイルを仮定しています。  
match は片方が稠密な頂点、もう片方が lowPoly された頂点で有る入力を仮定していますが、現行の実装では、lowFreq の出力がよろしくない結果を返します。

# バグ
PseudoBump はもっともらしいバンプマップを返しますが、正しくない場合があります。

# その他のダウンロードサイト
* https://ja.osdn.net/projects/goki-check/
* https://www.sourceforge.net/projects/gokicheck/
* https://konbu.sakura.ne.jp/files/goki_check_cc-1.01-lack-rotate-stable11.tar.gz
* http://files.limpid-intensity.info/goki_check_cc-1.01-lack-rotate-stable11.tar.gz
