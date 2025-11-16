#!/usr/bin/env bash

read -p "Enter the number of molecules: " molecule_num
if ! [[ "$molecule_num" =~ ^[0-9]+$ ]]; then
  echo "Error: Please enter a valid number." >&2
  exit 1
fi
echo "Using NUM = ${molecule_num}."

parent_directory=""

first_pdb=$(ls "${parent_directory}/Model01"/*.pdb 2>/dev/null | head -n 1)
if [ -z "$first_pdb" ]; then
    echo "Error" >&2
    exit 1
fi

dg_count=$(grep -o " DG " "$first_pdb" | wc -l)
g_count=$(grep -o " G " "$first_pdb" | wc -l)

if [ "$dg_count" -gt 0 ]; then
    calculated_mw=$(echo "scale=2; ($dg_count / 22) * 347" | bc)
    echo "MW = $calculated_mw"
elif [ "$g_count" -gt 0 ]; then
    calculated_mw=$(echo "scale=2; ($g_count / 22) * 363" | bc)
    echo "MW = $calculated_mw"
else
    echo "Warning: No guanine found."
    calculated_mw=5000
fi

shopt -s nullglob
mtz_candidates=( "${parent_directory}"/*.mtz )
shopt -u nullglob
if (( ${#mtz_candidates[@]} == 0 )); then
  echo "Error: No .mtz file found in ${parent_directory}" >&2
  exit 1
elif (( ${#mtz_candidates[@]} == 1 )); then
  mtz_file="${mtz_candidates[0]}"
else
  mtz_file="$(ls -t "${parent_directory}"/*.mtz | head -n1)"
  echo "Warning: Multiple .mtz files found. Using the most recent: ${mtz_file}" >&2
fi
echo "Reflection file to use: ${mtz_file}"

mkdir -p "${parent_directory}/Phaser-MR"
file_path="${parent_directory}/Phaser-MR"
cd "${file_path}" || exit 1

pdb_directory_1="${file_path}/Model01"
finish_directory_1="${pdb_directory_1}/finish"
mkdir -p "${pdb_directory_1}" "${finish_directory_1}"
cp "${parent_directory}/Model01"/*.pdb "${pdb_directory_1}" 2>/dev/null || true

for pdb_file_1 in "${pdb_directory_1}"/*.pdb; do
  [ -f "${pdb_file_1}" ] || continue
  base_name_1=$(basename "${pdb_file_1}" .pdb)
  output_directory="phaser-${base_name_1}"
  mkdir -p "${output_directory}"
  phaser <<EOF
TITLe ${base_name_1}
MODE MR_AUTO
HKLIn ${mtz_file}
ENSEmble ${base_name_1} PDB ${pdb_file_1} IDENtity 80
COMPosition NUCLeic MW ${calculated_mw} NUM ${molecule_num}
SEARch ENSEmble ${base_name_1} NUM ${molecule_num}
EOF
  [ -f PHASER.sol ]   && mv PHASER.sol   "${output_directory}/"
  [ -f PHASER.1.mtz ] && mv PHASER.1.mtz "${output_directory}/"
  [ -f PHASER.1.pdb ] && mv PHASER.1.pdb "${output_directory}/"
  cp "${pdb_file_1}" "${finish_directory_1}/"
done

llg_extracted="./extracted_data_LGG.txt"
llg_output_file="./TopLLG.txt"
: > "$llg_extracted"; : > "$llg_output_file"
for folder_path in ./phaser-*/; do
  [ -d "$folder_path" ] || continue
  phaser_sol="${folder_path}PHASER.sol"
  [ -f "$phaser_sol" ] || continue
  folder_name=$(awk 'NR==1{print $2; exit}' "$phaser_sol")
  solu_set_line=$(awk '/SOLU SET/ {print; exit}' "$phaser_sol")
  [ -n "$solu_set_line" ] || continue
  result=$(echo "$solu_set_line" | awk '{
    sep="";
    for(i=1;i<=NF;i++){
      if($i ~ /LLG=/){
        gsub(/[^0-9.]/,"",$i);
        printf "%s%s", sep, $i;
        sep=",";
      }
    }
    printf "\n"
  }')
  [ -n "$result" ] || continue
  printf '\n%s\n%s\n' "$folder_name" "$result" >> "$llg_extracted"
  max_value=$(echo "$result" | awk -F',' '{max=$1; for(i=2;i<=NF;i++) if($i>max) max=$i; print max}')
  printf '\n%s,%s\n' "$folder_name" "$max_value" >> "$llg_output_file"
done

tfz_extracted="./extracted_data_TFZ.txt"
tfz_output_file="./TopTFZ.txt"
: > "$tfz_extracted"; : > "$tfz_output_file"
for folder_path in ./phaser-*/; do
  [ -d "$folder_path" ] || continue
  phaser_sol="${folder_path}PHASER.sol"
  [ -f "$phaser_sol" ] || continue
  folder_name=$(awk 'NR==1{print $2; exit}' "$phaser_sol")
  solu_set_line=$(awk '/SOLU SET/ {print; exit}' "$phaser_sol")
  [ -n "$solu_set_line" ] || continue
  result=$(echo "$solu_set_line" | awk '{
    sep="";
    for(i=1;i<=NF;i++){
      if($i ~ /TFZ=/){
        gsub(/[^0-9.]/,"",$i);
        printf "%s%s", sep, $i;
        sep=",";
      }
    }
    printf "\n"
  }')
  [ -n "$result" ] || continue
  printf '\n%s\n%s\n' "$folder_name" "$result" >> "$tfz_extracted"
  max_value=$(echo "$result" | awk -F',' '{max=$1; for(i=2;i<=NF;i++) if($i>max) max=$i; print max}')
  printf '\n%s,%s\n' "$folder_name" "$max_value" >> "$tfz_output_file"
done

{
  printf 'Model TopLLG TopTFZ\n'
  awk -F',' 'NR==FNR{llg[$1]=$2; next} {print $1, llg[$1], $2}' OFS=' ' TopLLG.txt TopTFZ.txt
} > "${file_path}/results.txt"

echo finish
exit 0
