#!/usr/bin/env bash

# bash strict mode - snakemake.readthedocs.io/en/stable/project_info/faq.html
set -euo pipefail
IFS=$'\n\t'

usage_error() {
  echo "usage: $0 -i key_file -r user@host pwiz_image raw_file mzml_file" >&2
  exit 1
}

while getopts ":i:r:" opt; do
  case ${opt} in
    i ) key_file="$OPTARG"
      ;;
    r ) remote_addr="$OPTARG"
      ;;
    \? ) usage_error
      ;;
    : ) usage_error
      ;;
  esac
done

shift "$((OPTIND-1))"
if [[ $# -ne 3 ]]
then usage_error
fi

pwiz_image=$1
raw_file_path=$2
mzml_file_path=$3

if [[ -z "$key_file" || -z "$remote_addr" || \
      -z "$pwiz_image" || -z "$raw_file_path" || -z "$mzml_file_path" ]]
then usage_error
fi

if [[ "$raw_file_path" != /* ]]
then echo "raw file must be an absolute path" >&2; exit 1;
fi

if [[ "$mzml_file_path" != /* ]]
then echo "mzML file must be an absolute path" >&2; exit 1;
fi

raw_file_name=$(basename "$raw_file_path")
mzml_file_name=$(basename "$mzml_file_path")

echo "msconvert starting"

echo "checking ProteoWizard Docker image: '$pwiz_image'"
if [[ "$(docker images -q $pwiz_image 2>/dev/null)" == "" ]]
then docker pull -q $pwiz_image
fi

# The ProteoWizard Docker image needs a writable output directory and
# won't accept links, so we copy the input file to a temp directory.
TWD=$(mktemp -d)
echo "preparing temp directory: '$TWD'"
trap "{ rm -rf $TWD; }" EXIT
chmod a+w "$TWD"

tmp_raw_file="$TWD/$raw_file_name"
echo "creating temp raw file: '$tmp_raw_file'"
rsync -e "ssh -i $key_file" "$remote_addr:$raw_file_path" "$tmp_raw_file"

echo "converting raw file: '$raw_file_name'"
# use options -z and --gzip to compress output mzML file
# use --singleThreaded for reasons discussed in: https://github.com/ProteoWizard/pwiz/issues/684
docker run --rm -e WINEDEBUG=-all -v "$TWD":/data $pwiz_image \
  wine msconvert --mzML -z --gzip --singleThreaded "$raw_file_name"

echo "checking mzML file: '$mzml_file_name'"
tmp_mzml_file="$TWD/$mzml_file_name"
gunzip -c "$tmp_mzml_file" | xmllint --stream --noout -

echo "transferring mzML file '$mzml_file_name'"
mzml_dir=$(dirname "$mzml_file_path")
rsync -I -e "ssh -i $key_file" --rsync-path="mkdir -p $mzml_dir && rsync" \
  "$tmp_mzml_file" "$remote_addr:$mzml_file_path"

echo "msconvert done"
