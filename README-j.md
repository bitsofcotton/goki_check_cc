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
* * rstp    : ぼやけ具合を検出する際に使用される実際の画像の幅もしくは高さです。
* * vmax    : pseudoBumpVec で返される点の最大数です。(若干上下します)
* tilt.hh
* * z_ratio : [0,1] から [0,z_atio] への線形写像。
* scancontext.hh
* * matchPartialPartial::ndiv    : 合致する角度の分割数です
* * matchPartialPartial::thresh  : 合致を集めてくる際の平行なベクトルのための誤差です。
* * matchPartialPartial::thresht : 合致を集めてくる際の平行なベクトルの長さの比率に対する誤差です。
* * matchPartialPartial::threshp : 検出される最小の合致する点の数です。
* * matchPartialPartial::threshr : 検出される最小の画像倍率です。
* * matchPartialPartial::threshs : 合致が同じかどうか判定する閾値です。

# 文脈
写真の後でのピント調整プログラムに刺激されました。
また、この分野に関して、様々な(これと異なる)付帯条件での先行がたくさんありました。
(例えば、たくさんのカメラによる画像を使用するものや、あらかじめレイヤ毎に用意しておくもの、球状に膨らませるもの、
 動画から生成するものなどです。ただし、正しい公式ページがどれなのかわからないために、リンクは差し控えさせていただきます。)

# 状態
回転のある場合の全体の一致の実装をしています。

# 使い方
    make tools
    
    # enlarge
    ./tools enlarge input.ppm output.ppm
    
    # detect edges with absolute differential with dft value.
    ./tools collect input.ppm output.ppm
    
    # make 2d to 3d pseudo bumpmap
    ./tools bump input.ppm output.ppm output-delaunay.ppm output.obj
    
    # make tilts from original and bumpmap images.
    ./tools tilt input.ppm output-base input-bump.ppm
    
    # make (pseudo) lowpoly and get match.ppm and .obj file.
    ./tools lpoly input.ppm output-match.ppm output.obj
    
    # list matches.
    ./tools match input-matchbase.ppm output-base input-tobematched.ppm
    
    # list matches 2d - 3d.
    ./tools match3d input-matchbase.ppm output-base input-tobematched.obj
    
# ライブラリとしての使い方
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> input;
    
    #include "enlarge.hh"
    enlarger2ex<float> enlarger;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> enlarged(enlarger.compute(input, enlarger.ENLARGE_BOTH));
    
    enlarger2ex<float> detect;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> edgecollect(detect.compute(input, detect.COLLECT_BOTH));
    
    #include "fisheye.hh"
    PseudoBump<float> bump;
    std::vector<Eigen::Matrix<float, 3, 1> > points;
    std::vector<Eigen::Matrix<int,   3, 1> > delaunay;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> bumps;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> bumpd(bump.getPseudoBumpVec(input, points, delaunay, bumps));
    
    #include "tilt.hh"
    tilter<float> tilt;
    Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> tilted(tilt.tilt(input, bumpped, 0, 8, .95));
    
    #include "scancontext.hh"
    lowFreq<float> lf;
    matchPartialPartial<float> statmatch;
    std::vector<Eigen::Matrix<float, 3, 1> > shape0(lf.getLowFreq(input, 300));
    std::vector<Eigen::Matrix<float, 3, 1> > shape1(lf.getLowFreq(input2, 300));
    std::vector<match_t<float> > matches(statmatch.match(shape0, shape1));
    // match operations.
    
    // If you need, please scope with namespace block.
    // but include guard definition may harms.

# デモ
https://services.limpid-intensity.info/ にサンプルがあります。  
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
matchPartialPartial クラスは全体に対して安定な合致を返しますが、実際に必要なのは、単連結な部分に対して安定で、その他の部分に対して関連しない合致です。

# その他のダウンロードサイト
* https://ja.osdn.net/projects/goki-check/
* https://www.sourceforge.net/projects/gokicheck/
* https://konbu.sakura.ne.jp/files/goki_check_cc-1.01-lack-rotate-release3.tar.gz
* https://files.limpid-intensity.info/goki_check_cc-1.01-lack-rotate-release3.tar.gz
