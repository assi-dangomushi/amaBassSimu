#!/usr/bin/octave --no-init-file

## Copyright (C) 2022 Amanogawa Audio Labo
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

# AmanogawaBassSimulation Version 0.02 (2022/10/02)
versionNo = "AmaBassSimu ver0.02";

pkg load signal;


function g = calcAndDisplay(g)
  isMultiMic = strcmp(g.isMultiMic, "Y");

  figNum = 1;
  figure(figNum); clf; hold on; figNum = figNum + 1;
  plotRoom(g);
  
  # calc ir    ira: mic1 irb(1:16): mic16
  [ira, irb, iraSumSp, irbSumSp] = getAllIr(g);

  # calc h 
  [w, ha, hb, hMean] = calcH(g, iraSumSp, irbSumSp);

  ################################### plot #############################3
  # plot freq
  figure(figNum); clf; hold on; figNum = figNum + 1;
  plotFreq(g, w, ha, hb, hMean, "non-Eq");

  # plot time domain
  figNum = plotTimedomain(figNum, g, iraSumSp, irbSumSp, "non-Eq");
  uiwait(msgbox('Operation Completed','Success'));

  # Eq
  figNum2 = figNum;
  while (true)
    switch menu("Eq", {"edit Eq", "plot Total Eq", "plot Individual Eq"})
    case 1 # edit Eq
      g = editstructure(g, {"nSmooth", "eqLength", "irFreqs"});
    case 2 # total Eq
      figNum = figNum2;
      eq1 = mkEq(g, w, hMean);
      figure(figNum); clf; hold on; figNum = figNum + 1;
      eqPlot(g, eq1);
      title("Total Eq-Curve", "FontSize", 25);
      iraSumSp2 = fftconv(eq1, iraSumSp);
      irbSumSp2 = {};
      if isMultiMic
        for micNum = 1:16
          irbSumSp2{micNum} = fftconv(eq1, irbSumSp{micNum});
        endfor
      endif  
      [w, ha2, hb2, hMean2] = calcH(g, iraSumSp2, irbSumSp2);
      # plot freq
      figure(figNum); clf; hold on; figNum = figNum + 1;
      plotFreq(g, w, ha2, hb2, hMean2, "Total Equalized");

      # plot time domain
      figNum = plotTimedomain(figNum, g, iraSumSp2, irbSumSp2, "Total Equalized");
      uiwait(msgbox('Operation Completed','Success'));
    case 3 # Individual Eq
      figNum = figNum2;
      inEq = mkInEq(g, ira, irb);
      figure(figNum); clf; hold on; figNum = figNum + 1;
      eqPlot(g, inEq);
      title("Individual Eq-Curve", "FontSize", 25);
      [iraSumSp2, irbSumSp2] = applyIndividualEq(g, inEq, ira, irb);
      [w, ha2, hb2, hMean2] = calcH(g, iraSumSp2, irbSumSp2);
      # plot freq
      figure(figNum); clf; hold on; figNum = figNum + 1;
      plotFreq(g, w, ha2, hb2, hMean2, "Individual Equalized");
      # plot time domain
      figNum = plotTimedomain(figNum, g, iraSumSp2, irbSumSp2, "Individual Equalized");
      uiwait(msgbox('Operation Completed','Success'));
    otherwise
      break;
    endswitch
  endwhile
endfunction

function ir = calcir(roomSize, roomRef, spPosition, micPosition, n = 50, sec =1)
  p = sprintf("%f ", [roomSize, roomRef, spPosition, micPosition, n, sec]);
  fp = popen(["./calcir ", p], "r");
  ir = fread(fp, "float64"); 
  pclose(fp);
endfunction

function micPositions = calcMicPositions(g)
  micPositions = zeros(16,3);
  for micNum = 1:16
    micPositions(micNum, 1) = mod(micNum-1, 4) * (g.mic.area/3) + (g.mic.position(1) - g.mic.area/2);
    micPositions(micNum, 2) = mod(fix((micNum-1)/4), 4) * (g.mic.area/3) + (g.mic.position(2) - g.mic.area/2);
    micPositions(micNum, 3) = g.mic.position(3);
  endfor
endfunction


function spDelays = calcSpDelay(g)
  spDistance = zeros(1, g.spNum);
  spDelays = zeros(g.spNum, 1);
  for n = 1:g.spNum
    spDistance(n) = norm(g.sp.positions(n,:) .- g.mic.position);
    spDelays(n) = - (floor((spDistance(n) / 340) * g.fs));
  endfor
  spDelays = spDelays - min(spDelays);
endfunction


# plot freq
function plotFreq(g, w, ha, hb, hMean, titleStr)
  w2 = w(w <= g.fRange(2));
  if strcmp(g.fScale, "linear")
    fp = @plot;
  elseif strcmp(g.fScale, "log")
    fp = @semilogx;
  else
    return
  endif
  
  if strcmp(g.isMultiMic, "Y")
    for micNum = 1:16
      fp(w2, mag2db(abs(hb{micNum}(1:length(w2)))), "Color", [0.7, 0.7, 0.7]);
    endfor
    fp(w2, mag2db(hMean)(1:length(w2)), "Color", [0, 0, 0], "LineWidth", 1.5);
  endif
  fp(w2, mag2db(abs(ha)(1:length(w2))), "1");
  axis([g.fRange(1), g.fRange(2), g.dRange(1), g.dRange(2)]);
  grid minor on;
  title(titleStr, "FontSize", 25);
endfunction

# plot Room
function plotRoom(g)
  
  isMultiMic = strcmp(g.isMultiMic, "Y");
  
  hold on;

  #speaker
  for n = 1:g.spNum
    plot3 (g.sp.positions(n,1), g.sp.positions(n,2), g.sp.positions(n,3), '4s', 'markerSize', 10);
  endfor

  axis('equal', 'ij');
  grid on;
  
  #mic 
  plot3(g.mic.position(1), g.mic.position(2), g.mic.position(3), 'x', 'markerSize', 10);
  if isMultiMic
    micPositions = calcMicPositions(g);
    for n = 1:16
      plot3(micPositions(n,1), micPositions(n,2), micPositions(n,3), 'o', 'markerSize', 5, "Color", [0.7, 0.7, 0.7]);
    endfor
  endif

  #room
  plot3([0, g.room.width], [0, 0], [0, 0], '3' , 'LineWidth', 2);
  plot3([0, g.room.width], [g.room.depth, g.room.depth], [0, 0], '3', 'LineWidth', 2);
  plot3([0, 0], [0, g.room.depth], [0, 0], '3', 'LineWidth', 2);
  plot3([g.room.width, g.room.width], [0, g.room.depth], [0, 0], '3', 'LineWidth', 2);
  plot3([0, g.room.width], [0, 0], [g.room.hight, g.room.hight], '1', 'LineWidth', 2 );
  plot3([0, g.room.width], [g.room.depth, g.room.depth], [g.room.hight, g.room.hight], '5');
  plot3([0, 0], [0, g.room.depth], [g.room.hight, g.room.hight], '5' );
  plot3([g.room.width, g.room.width], [0, g.room.depth], [g.room.hight, g.room.hight], '5' );
  plot3([0, 0], [0, 0], [0, g.room.hight], '1', 'LineWidth', 2);
  plot3([g.room.width, g.room.width], [0, 0], [0, g.room.hight], '1', 'LineWidth', 2);
  plot3([0, 0], [g.room.depth, g.room.depth], [0, g.room.hight], '5');
  plot3([g.room.width, g.room.width], [g.room.depth, g.room.depth], [0, g.room.hight], '5');
endfunction


function eqPlot(g, eq1)
    
  if strcmp(g.fScale, "linear")
    fp = @plot;
  elseif strcmp(g.fScale, "log")
    fp = @semilogx;
  else
    return
  endif

  eq2 = {};
  if iscell(eq1) 
    eq2 = eq1;
  else
    eq2{1} = eq1;
  endif

  for eqNum = 1:numel(eq2)
    [h, w] = freqz(eq2{eqNum}, 1, 8192*64, g.fs);
    w2 = w(w <= g.fRange(2));
    fp(w2, mag2db(abs(h)(1:length(w2))), "3");
  endfor
  axis([g.fRange(1), g.fRange(2), g.dRange(1), g.dRange(2)]);
  grid minor on;
endfunction

# calc h
function  [w, ha, hb, hMean] = calcH(g, iraSumSp, irbSumSp)
  if length(iraSumSp) > 8192*64
    iraSumSp = iraSumSp(1:8192*64)
  endif
  [h, w] = freqz(iraSumSp, 1, 8192*64, g.fs);
  ha = h;
  hb = {};
  if strcmp(g.isMultiMic, "Y")
    hMean = zeros(length(w),1);
    for micNum = 1:16
      if length(irbSumSp{micNum}) > 8192*64
        irbSumSp{micNum} = irbSumSp{micNum}(1:8192*64);
      endif
      [h, w] = freqz(irbSumSp{micNum}, 1, 8192*64, g.fs);
      hb{micNum} = h;
      hMean = hMean .+ abs(hb{micNum});
    endfor
    hMean = hMean / 16;
  else
    hMean = ha;
  endif
endfunction

function [ira, irb, iraSumSp, irbSumSp] = getAllIr(g)
  switch g.taMode
    case "auto"
      delay = calcSpDelay(g);
    case "none"
      delay = zeros(g.spNum, 1);
    case "manual"
      delay = g.sp.delays;
    otherwise  
      questdlg("taMode is  [auto|manual|none]");
      return
  endswitch
  maxDelay = max(delay);
  
  ira = {};
  for spNum = 1:g.spNum
    ir = calcir([g.room.width, g.room.depth, g.room.hight],
                   [g.room.refLeft, g.room.refRight, g.room.refFront, g.room.refBack, g.room.refBotom, g.room.refTop],
                   [g.sp.positions(spNum,:)],
                   g.mic.position, g.refNum, g.irDuration);
    ir = [zeros(delay(spNum),1); ir; zeros(maxDelay - delay(spNum), 1)] * g.sp.gains(spNum);
    ira{spNum} = ir;
  endfor  
  iraSumSp = ira{1};
  if g.spNum >= 2
    for spNum = 2:g.spNum
      iraSumSp = iraSumSp .+ ira{spNum};
    endfor
    iraSumSp = iraSumSp / g.spNum;
  endif

  irb = {}; irbSumSp = {};
  if strcmp(g.isMultiMic, "Y")
    micPositions = calcMicPositions(g);
    for spNum = 1:g.spNum
      for micNum = 1:16
        ir = calcir([g.room.width, g.room.depth, g.room.hight],
                   [g.room.refLeft, g.room.refRight, g.room.refFront, g.room.refBack, g.room.refBotom, g.room.refTop],
                   [g.sp.positions(spNum,:)],
                   micPositions(micNum,:), g.refNum, g.irDuration);
        ir = [zeros(delay(spNum),1); ir; zeros(maxDelay - delay(spNum), 1)] * g.sp.gains(spNum);
        irb{spNum, micNum} = ir;
      endfor
    endfor
    for micNum = 1:16
      irbSumSp{micNum} = irb{1, micNum};
      if g.spNum >= 2
        for spNum = 2:g.spNum
          irbSumSp{micNum} = irbSumSp{micNum} .+ irb{spNum, micNum};
        endfor
        irbSumSp{micNum} = irbSumSp{micNum} / g.spNum;
      endif
    endfor
  endif
endfunction

function eq = mkEq(g, w, h)
  w20 = f2w(w, 20);
  w500 = f2w(w, 500);
  hs = fsmooth(h, g.nSmooth);
  hs(1:w20) = hs(w20);
  hs(w500:end) = hs(w500);
  hs = 1 ./ hs;
  eq = mifft(hs);
  eq = fftshift(eq);
  eq = l2m(eq);
  eq = ircut(eq, 0, g.eqLength);
endfunction

function inEq = mkInEq(g, ira, irb)
  inEq = {};
  if strcmp(g.isMultiMic, "Y")
    for spNum = 1:g.spNum
      hb = zeros(8192*64, 1);
      for micNum = 1:16
        [h, w] = freqz(irb{spNum, micNum}, 1, 8192*64, g.fs);
        hb = hb .+ abs(h);
      endfor
      hb = hb / 16;
      inEq{spNum} = mkEq(g, w, hb);
    endfor
  else
    for spNum = 1:g.spNum
      [h, w] = freqz(ira{spNum}, 1, 8192*64, g.fs);
      inEq{spNum} = mkEq(g, w, h);
    endfor
  endif
endfunction

function figNum = plotTimedomain(figNum, g, iraSumSp, irbSumSp, titleStr)
  for freq = g.irFreqs
    m = morlet2(freq, g.fs);
    plotLength = min([length(m)*7, length(iraSumSp)]);
    figure(figNum); clf; hold on; figNum = figNum + 1;
    if strcmp(g.isMultiMic, "Y")
      for micNum = 1:16
        plot(fftconv(irbSumSp{micNum}(1:plotLength),m)(1:end-length(m)), "Color", [0.7, 0.7, 0.7]);
      endfor
    endif
    plot(fftconv(iraSumSp(1:plotLength),m)(1:end-length(m)), '1');   
    title([num2str(freq), " Hz ", titleStr], "FontSize", 25);
  endfor    
endfunction

function [iraSumSp2, irbSumSp2] = applyIndividualEq(g, inEq, ira, irb)
  irb2 = {};
  for spNum = 1:g.spNum
    if strcmp(g.isMultiMic, "Y")
      for micNum = 1:16
        irb2{spNum, micNum} = fftconv(inEq{spNum}, irb{spNum, micNum});
      endfor
    endif
    ira2{spNum} = fftconv(inEq{spNum}, ira{spNum});
  endfor

  iraSumSp2 = ira2{1};
  if g.spNum >= 2
    for spNum = 2:g.spNum
      iraSumSp2 = iraSumSp2 .+ ira2{spNum};
    endfor
    iraSumSp2 = iraSumSp2 / g.spNum;
  endif
  
  irbSumSp2 = {};
  if strcmp(g.isMultiMic, "Y")  
    for micNum = 1:16
      irbSumSp2{micNum} = irb2{1, micNum};
      if g.spNum >= 2
        for spNum = 2:g.spNum
          irbSumSp2{micNum} = irbSumSp2{micNum} .+ irb2{spNum, micNum};
        endfor
        irbSumSp2{micNum} = irbSumSp2{micNum} / g.spNum;
      endif
    endfor
  endif

endfunction

#################################################################

#hの絶対値を1/nオクターブで平滑化
#ガウス加重を使って平滑化
#h,wはfreqzで計算されたもの
#角度情報は失われる
function hs=fsmooth(h,n)
  # hを絶対に変更
  h = abs(h);

  # wを作る （end は範囲を拡張）
  w = (0:1/length(h):1)';
  w(end) = 3; 
  # fs/2を補完
  h(end+1) = h(end);
 
  # wlog を作る #
  
  # 1octをlogでn1等分
  n1 = 1024;
  t = 2 ^ (1/n1);
  sf = w(2)/2;
  ef = 2;
  a = log(ef/sf)/log(t);

  x = 0:ceil(a);
  wlog = (sf * t .^ x)';

  # interp1 で wlog 対応のh2を作る
  h2 = interp1(w,h,wlog);

  # 長さ1/n oct のローパスフィルタ（ガウス）を作る
  n2 = ceil(n1 / n / 2) * 2 + 1; # n2は奇数
  b = 4;    
  low = gausswin(n2); 
  low = db2mag(-frmax(low))*low;

  #ローパスフィルタを畳み込んで平滑化
  h3 = fftconv(h2,low);
  
  #同じ位置を切り出す
  h3 = h3((n2-1) / 2 + 1 : end);
  h3 = h3 (1 : length(h2));

  # w に対応したhsに戻す w(1)は範囲外なので除外して戻す
  hs = interp1(wlog, h3, w(2:end-1));
  hs = [hs(1);hs];  
endfunction

#freqzで計算したw(n)で周波数fを超えない最大のnを返す
function n = f2w(w,f)
 n=floor(f/w(2))+1;
 if (n > length(w))
   n=length(w);
 endif
endfunction

#インパルスレスポンスaを前後pre、postサンプルで切り出して
#窓関数wをかけて返す。 @hamming など
#窓関数を省略した場合には6次のmado1関数を使用
function b = ircut (a, pre, post,w)
  if !(exist("w"))
    w=@(x)mado1(x,6);
  endif

  if (size(a)(2)>1)
    a=a';
  endif

  [m1 m2] = max (abs(a));
  b = a(m2-pre:m2+post-1);
  m = hmado(length(b),pre+1,w);
  b = b .* m;
endfunction

#y=x^tの形の窓をかける
function a = mado1(l,t)
  a=[ones(1,l)];
  for n=1:l
    a(n)=((n-l/2)/l)^t;
  endfor
  a=-a/(0.5^t)+1;
  a=a';
endfunction

#xの最大値(dB)とその周波数を返す
function [m,f]=frmax(x,fs);
  fs1=pi*2;
  if (exist("fs"))
    fs1=fs;
  endif
  [h,w]=freqz(x,1,8192*16,fs1);
  h=mag2db(abs(h));
  [m,b]=max(h);
  f=w(b);
endfunction

# 構造体s1の要素をGUIで編集して新しい構造体s2を返す
function  s2 = editstructure(s1, f = fieldnames(s1)', n = 14)
  prompt = f;
  s2 = s1;
  #デフォルト値を文字列に変換
  defaults = cell(0);
  for i = prompt
    j = cell2mat(i);
    if isnumeric (s1.(j))
      #defaults{end+1} = num2str(s1.(j),16);
      defaults{end+1} = num2str(s1.(j),"%.16g  ");
    else
      defaults{end+1} = s1.(j);
    endif
  endfor
  colrow = [];
  for i = prompt
    colrow = [colrow; 1, n];
  endfor
  dims = inputdlg (prompt, ["Edit ", inputname(1)], colrow, defaults);
  if isempty(dims)
    return;
  endif  
  # 元が数値なら数値に変換
  n = 1; 
  for  i = prompt
    j = cell2mat(i);
    if ~(ischar (s1.(j)))
      s2.(j) = str2num(dims{n});
    else
      s2.(j) = dims{n};
    endif
    n = n + 1;
  endfor
endfunction

# 配列をGUIで編集
function a2  = editarray(a1)
  arrayname = inputname(1);
  s = struct;
  n = size(a1)(1);
  for j = 1:n
    s.([arrayname, sprintf("%02d", j)]) = a1(j,:);
  endfor
  s = editstructure(s);
  for j = 1:n
    a2(j,:) = s.([arrayname, sprintf("%02d", j)]);
  endfor
endfunction

# freqzで計算したhにエイリアシング成分を足して逆fftする。
function b=mifft(a)
  b=flipud(a);
  b=conj(b);
  c=abs(a(end));
  b=[a;c;b(1:end-1)];
  b=real(ifft(b));
endfunction

#　直線位相（に限らないけど）フィルタを最小位相フィルタに変換
function m = l2m(l)
  if (size(l)(1) ~= 1)
    fl=1;
    l=l';
  else
    fl=0;
  endif
  n=length(l);
  l2=[l,zeros(1,n*2)];
  [y m]=rceps(l2);
  m=m(1:n);
  if (fl==1)
    m=m';
  endif
endfunction

#長さL,中心tの非対称窓
#w:窓関数 @hamming など
function a = hmado(L,t,w)
  if (t != 1)
    w1=w(2*(t-1));
    w1=w1(1:t-1);
    w2=w(2*(L-t));
    w2=w2(L-t:end);
    a=[w1;w2];
  else
    w2=w(2*(L-t));
    w2=w2(L-t:end);
    a=w2;
  endif
endfunction

function Y = morlet2(f,fs)
  [Y, X] = morlet (-4, 4, floor(8*fs*5/(2*pi)/f));
endfunction

########################################
# main 
########################################
feet = 0.3048;
maxSpNum = 50;

# Global Structure "g"
# default value
if !(exist("g"))
  g = struct();
  g.versioNo = versionNo;
  g.refNum = 150;
  g.isMultiMic = "Y"; # "Y" or "N" 
  g.spNum = 4; # 1..maxSpNum
  g.taMode = "none"; # auto manual none
  g.irFreqs = [50, 70, 100]; # time domain
  g.irDuration = 5; 
  g.fScale = "linear"; # log | linear
  g.fRange = [1, 100];
  g.dRange = [-30, 25];
  g.nSmooth = 6; # 1/6 Oct Smoothing 
  g.eqLength = 8192;
  g.fs = 96000; # Fixed!

  g.room = struct();
  g.room.width = 20*feet;
  g.room.depth = 24*feet;
  g.room.hight = 9*feet;
  g.room.refLeft = 0.95;
  g.room.refRight = 0.95;
  g.room.refFront = 0.95;
  g.room.refBack = 0.95;
  g.room.refBotom = 0.95;
  g.room.refTop = 0.95;

  g.mic = struct();
  g.mic.position = [g.room.width /2, g.room.depth/2, 1.1];
  g.mic.area = 6 * feet; 
  
  g.sp = struct;
  g.sp.positions = zeros(maxSpNum,3); # W, D, H
  g.sp.gains = ones(maxSpNum,1);
  g.sp.delays = zeros(maxSpNum, 1);  
endif

while(true)
  s = { "global",
        "Edit Room",
        "Mic positon",
        "Speaker position",
        "Speaker gain",
        "Time Alignment",
        "Load",
        "Save",
        "Plot room",
        "Calc and Display",
      };
  switch menu(versionNo, s)
    case 1 #global
      g = editstructure(g, {"refNum", "isMultiMic", "spNum", "taMode", "irFreqs", "fScale", "fRange", ...
                            "dRange", "nSmooth", "eqLength", "irDuration", });
    case 2 #room
      g.room = editstructure(g.room);
    case 3 #mic
      g.mic = editstructure(g.mic);
    case 4 #Sp position
      g.sp.positions((1:g.spNum),:) = editarray(g.sp.positions((1:g.spNum),:));
    case 5 # SP gain
      g.sp.gains(1:g.spNum) = editarray(g.sp.gains(1:g.spNum));
    case 6 #TA
      while(true)
        switch menu("TA menu", {"calc delay", "clear", "edit delay"})
          case 1 #calc delay
            g.sp.delays(1:g.spNum) = calcSpDelay(g);
          case 2 #clear delay
            g.sp.delays(1:g.spNum) = 0;
          case 3 #edit delay
            g.sp.delays(1:g.spNum) = editarray(g.sp.delays(1:g.spNum));
          otherwise
            break;
        endswitch
      endwhile
    case 7 #load
      [FNAME, FPATH, FLTIDX] = uigetfile ("./data/*.wooferconf");
      if FNAME ~= 0
        load ([FPATH, FNAME]);
      endif
    case 8 #save
      [FNAME, FPATH, FLTIDX] = uiputfile ("./data/*.wooferconf");
      if FNAME ~= 0
        save ([FPATH, FNAME], "g") ;
      endif
    case 9 # Plot Room
      figure(1); clf; hold on;
      plotRoom(g); 
      uiwait(msgbox('Operation Completed','Success'));
    case 10 #calc and display
      g = calcAndDisplay(g);
      #uiwait(msgbox('Operation Completed','Success'));
    otherwise
      break; 
  endswitch
endwhile

