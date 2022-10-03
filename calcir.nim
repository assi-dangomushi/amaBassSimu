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

# calcir.nim version 0.12 (2022/09/05))
# calcir room_W room_D room_H roomref_left roomref_right roomref_front roomref_back roomref_botom roomref_top sp_W sp_D sp_H mic_W mic_D mic_H numOfRef sec

import math, os, strutils

const
  Fs: int = 96000 # サンプリング周波数
  Ss: float64 = 340.0 #音速

type Pos3 = array[3, float64]

let
  v = commandLineParams()
  roomSize: array[3, float64] = [v[0].parseFloat, v[1].parseFloat, v[2].parseFloat]
  roomRef: array[6, float64] = [v[3].parseFloat, v[4].parseFloat, v[5].parseFloat, v[6].parseFloat, v[7].parseFloat, v[8].parseFloat] #左右前後下上
  spPosition: Pos3 = [v[9].parseFloat, v[10].parseFloat, v[11].parseFloat]
  micPosition: Pos3 = [v[12].parseFloat, v[13].parseFloat, v[14].parseFloat]
  n: int = (v[15].parseFloat).int
  sec: float = v[16].parseFloat
  L = (sec * Fs.float).int

func getDistance(p1, p2: Pos3): float64 =
  return sqrt((p1[0] - p2[0])^2 + (p1[1] - p2[1])^2 + (p1[2] - p2[2])^2)

proc getMirrorPosition(roomNo: array[3, int]): Pos3 =
  for i in 0..2:
    result[i] = roomNo[i].float64 * roomSize[i] + spPosition[i] + abs((roomNo[i] mod 2)).float64 * (roomSize[i] - 2 * spPosition[i])

var
  left = newSeq[float64](n+1)
  right = newSeq[float64](n+1)
  front = newSeq[float64](n+1)
  back = newSeq[float64](n+1)
  botom = newSeq[float64](n+1)
  top = newSeq[float64](n+1)
for i in 0..<(n+1):
  left[i] = roomRef[0]^i
  right[i] = roomRef[1]^i
  front[i] = roomRef[2]^i
  back[i] = roomRef[3]^i
  botom[i] = roomRef[4]^i
  top[i] = roomRef[5]^i

proc getRefCoef(roomNo: array[3, int]): float64 =
  let
    leftRef: float64 = left[abs(((roomNo[0] + 2*n) div 2) - n)]
    rightRef: float64 = right[abs(((roomNo[0] + 2*n + 1) div 2) - n)]
    frontRef: float64 = front[abs(((roomNo[1] + 2*n) div 2) - n)]
    backRef: float64 = back[abs(((roomNo[1] + 2*n + 1) div 2) - n)]
    botomRef: float64 = botom[abs(((roomNo[2] + 2*n) div 2) - n)]
    topRef: float64 = top[abs(((roomNo[2] + 2*n + 1) div 2) - n)]
  return leftRef * rightRef * frontRef * backRef * botomRef * topRef

proc calcir(): seq[float64] =
  var
    ir = newSeq[float64](L)
    mp: Pos3
    d, r, y: float64
    x: int

  for n1 in -n..n:
    for n2 in -(n-abs(n1))..(n-abs(n1)):
      for n3 in -(n-abs(n1)-abs(n2))..(n-abs(n1)-abs(n2)):
        mp = getMirrorPosition([n1, n2, n3])
        d = getDistance(mp, micPosition)
        r = getRefCoef([n1, n2, n3])
        x = (d / Ss * Fs.float64).int
        y = 1 / d * r
        if (x < L):
          ir[x] = ir[x] + y
  return ir

var ir = calcir()
discard stdout.writeBuffer(ir[0].addr, sizeof(float64)*L)
