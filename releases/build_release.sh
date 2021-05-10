#!/usr/bin/env bash

# version number of the build
version_number="0.1.0"

# folders to include in release build
folders=(\
    "../docs" \
    "../examples" \
    "../src" \
    "../utilities")

# specific files to include in release build
files=("LICENSE" \
    "README.md" \
    "changelog.md" \
    "EXCEED-DM-docs.md" \
    "fobos" \
    "install.md")

# create the folder
mkdir "v${version_number}"

# move all the necessary folders
for fol in "${folders[@]}"; do
    cp -r "$fol" "v${version_number}/"
done

# move all the necessary files
for file in "${files[@]}"; do
    cp "../${file}" "v${version_number}/"
done

# tar the files
cd "./v${version_number}"
tar -czvf "../EXCEED-DM-v${version_number}.tar.gz" .

cd ../
# remove the version folder (no reason to save all that)
rm -r "./v${version_number}"

