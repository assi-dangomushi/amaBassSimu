# amaBassSimu

multi Bass Simulation program

only Rectangle Room

plot Freq and TimeDomain

Equalizer: FIR, minimal Phase


## install

### install octave and octave-signal pkg

ubuntu 22.04

sudo apt install octave-signal

### install choosenim

sudo apt install curl  
curl https://nim-lang.org/choosenim/init.sh -sSf | sh

### install amaBassSimu

git clone https://github.com/assi-dangomushi/amaBassSimu.git 

### Compile calcir

cd amaBassSimu  
nim c -d:release --opt:speed calcir.nim

## Usage

./amaBassSimu.m

