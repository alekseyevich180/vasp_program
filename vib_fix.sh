#!/bin/bash

target_numbers=("$@")

add_fix.sh 1000

# POSCAR�t�@�C���̃p�X
poscar_file="POSCAR"

# �V�����t�@�C���̃p�X
output_file="POSCAR_new"

# �ʂ��ԍ���������
atom_index=1

# POSCAR�t�@�C����ǂݍ��݁A�V�����t�@�C���ɒʂ��ԍ���ǉ����ĕۑ�
while IFS= read -r line; do
  stripped_line=$(echo "$line" | awk '{$1=$1};1')  # �s�̐擪�Ɩ����̋󔒕�������菜��

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