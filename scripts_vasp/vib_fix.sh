#!/bin/bash

target_numbers=("$@")

add_fix.sh 1000

# POSCARファイルのパス
poscar_file="POSCAR"

# 新しいファイルのパス
output_file="POSCAR_new"

# 通し番号を初期化
atom_index=1

# POSCARファイルを読み込み、新しいファイルに通し番号を追加して保存
while IFS= read -r line; do
  stripped_line=$(echo "$line" | awk '{$1=$1};1')  # 行の先頭と末尾の空白文字を取り除く

  if [[ "$stripped_line" == *" F F F" ]]; then
    echo "${stripped_line} ${atom_index}" >> "$output_file"
    ((atom_index++))
  else
    echo "$line" >> "$output_file"
  fi
done < "$poscar_file"

mv POSCAR_new POSCAR

while IFS= read -r line; do
    if [[ $line == *"F F F"* ]]; then
        line_parts=($line)
        last_part=${line_parts[-1]}
        for target_number in "${target_numbers[@]}"; do
            if [[ $last_part == $target_number ]]; then
                line=${line//"F F F"/"T T T"}
                break
            fi
        done
    fi
    echo "$line"
done < POSCAR > POSCAR_new

mv POSCAR_new POSCAR