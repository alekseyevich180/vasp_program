#!/bin/bash

# awk を使用して格子定数を抽出
awk '
  NR == 3 {
    a = sqrt($1^2 + $2^2 + $3^2)
  }

  NR == 4 {
    b = sqrt($1^2 + $2^2 + $3^2)
  }

  NR == 5 {
    c = sqrt($1^2 + $2^2 + $3^2)
    printf "Lattice Constants (a, b, c): %f %f %f\n", a, b, c
    exit
  }
' "$1"
