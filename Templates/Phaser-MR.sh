#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
parent_directory="$(pwd)"
models_dir="${parent_directory}/Models"
work_dir="${parent_directory}/Phaser-MR"
finish_dir="${work_dir}/Models-Finish"
info_file="${parent_directory}/Info.txt"
mkdir -p "${work_dir}" "${finish_dir}"
if [[ ! -f "${info_file}" ]]; then
  echo "Error: ${info_file} not found. Please create it with lines like '--MW 1' and '--NUM 1'." >&2
  exit 1
fi
mw_value=""
num_value=""
while IFS= read -r line; do
  line="${line#"${line%%[![:space:]]*}"}"
  line="${line%"${line##*[![:space:]]}"}"
  [[ -z "${line}" ]] && continue
  [[ "${line}" =~ ^# ]] && continue
  key="$(awk '{print $1}' <<< "${line}")"
  val="$(awk '{print $2}' <<< "${line}")"
  case "${key}" in
    --MW|--mw|--Mw|--mW)
      mw_value="${val}"
      ;;
    --NUM|--num|--Num|--nUm|--nuM)
      num_value="${val}"
      ;;
  esac
done < "${info_file}"
if [[ -z "${mw_value}" ]] || [[ ! "${mw_value}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
  echo "Error: --MW must be provided in ${info_file} as a positive number (e.g., '--MW 1' or '--MW 1000')." >&2
  exit 1
fi
if [[ -z "${num_value}" ]] || [[ ! "${num_value}" =~ ^[1-9][0-9]*$ ]]; then
  echo "Error: --NUM must be provided in ${info_file} as an integer >= 1 (e.g., '--NUM 1')." >&2
  exit 1
fi
echo "Loaded composition parameters: MW=${mw_value}, NUM=${num_value}"
mtz_candidates=( "${parent_directory}"/*.mtz )
if (( ${#mtz_candidates[@]} == 0 )); then
  echo "Error: No .mtz file found in ${parent_directory}" >&2
  exit 1
fi
mtz_file="$(ls -t "${parent_directory}"/*.mtz | head -n1)"
echo "Reflection file to use: ${mtz_file}"
cd "${work_dir}"
pdb_files=( "${models_dir}"/*.pdb )
if (( ${#pdb_files[@]} == 0 )); then
  echo "Error: No .pdb files found in ${models_dir}" >&2
  exit 1
fi
for pdb_path in "${pdb_files[@]}"; do
  filename_noext="$(basename "${pdb_path}" .pdb)"
  base="${filename_noext%%_*}"   # Use part before the first underscore
  outdir="${work_dir}/phaser-${base}"
  mkdir -p "${outdir}"
  echo "==> Running Phaser for: ${base}"
  phaser <<EOF
TITLe ${base}
MODE MR_AUTO
HKLIn ${mtz_file}
ENSEmble ${base} PDB ${pdb_path} IDENtity 80
COMPosition NUCLeic MW ${mw_value} NUM ${num_value}
SEARch ENSEmble ${base} NUM 1
EOF
  for f in PHASER.sol PHASER.1.mtz PHASER.1.pdb; do
    [[ -f "$f" ]] && mv "$f" "${outdir}/"
  done
  mv "${pdb_path}" "${finish_dir}/"
done
llg_extracted="${work_dir}/extracted_data_LGG.txt"
llg_output="${work_dir}/TopLLG.txt"
tfz_extracted="${work_dir}/extracted_data_TFZ.txt"
tfz_output="${work_dir}/TopTFZ.txt"
: > "${llg_extracted}"; : > "${llg_output}"
: > "${tfz_extracted}"; : > "${tfz_output}"
for folder_path in "${work_dir}"/phaser-*/; do
  phaser_sol="${folder_path}PHASER.sol"
  [[ -f "${phaser_sol}" ]] || continue
  folder_name_raw="$(awk 'NR==1{print $2; exit}' "${phaser_sol}")"
  folder_name="${folder_name_raw%%_*}"
  solu_set_line="$(awk '/SOLU SET/ {print; exit}' "${phaser_sol}")"
  [[ -n "${solu_set_line}" ]] || continue
  llg_line="$(echo "${solu_set_line}" | awk '{
    sep="";
    for(i=1;i<=NF;i++){
      if($i ~ /LLG=/){
        gsub(/[^0-9.\-]/,"",$i);
        printf "%s%s", sep, $i; sep=","
      }
    }
    printf "\n"
  }')"
  if [[ -n "${llg_line}" ]]; then
    printf '\n%s\n%s\n' "${folder_name}" "${llg_line}" >> "${llg_extracted}"
    llg_max="$(echo "${llg_line}" | awk -F',' '{max=$1; for(i=2;i<=NF;i++) if($i>max) max=$i; print max}')"
    printf '%s,%s\n' "${folder_name}" "${llg_max}" >> "${llg_output}"
  fi
  tfz_line="$(echo "${solu_set_line}" | awk '{
    sep="";
    for(i=1;i<=NF;i++){
      if($i ~ /TFZ=/){
        gsub(/[^0-9.\-]/,"",$i);
        printf "%s%s", sep, $i; sep=","
      }
    }
    printf "\n"
  }')"
  if [[ -n "${tfz_line}" ]]; then
    printf '\n%s\n%s\n' "${folder_name}" "${tfz_line}" >> "${tfz_extracted}"
    tfz_max="$(echo "${tfz_line}" | awk -F',' '{max=$1; for(i=2;i<=NF;i++) if($i>max) max=$i; print max}')"
    printf '%s,%s\n' "${folder_name}" "${tfz_max}" >> "${tfz_output}"
  fi
done
{
  printf 'Model TopLLG TopTFZ\n'
  awk -F',' 'NR==FNR{llg[$1]=$2; next} {print $1, llg[$1], $2}' OFS=' ' "${llg_output}" "${tfz_output}"
} > "${work_dir}/results.txt"
echo "Done: Results have been saved under ${work_dir}"
echo " - Each model: phaser-<model_name>/"
echo " - Processed models: ${finish_dir}/"
echo " - Aggregated results: ${work_dir}/results.txt"
