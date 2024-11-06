#!/bin/bash

# awk ���g�p���Ċi�q�萔�𒊏o
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
