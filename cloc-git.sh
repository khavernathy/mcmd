#!/usr/bin/env bash

theurl="https://github.com/khavernathy/mcmd"

git clone --depth 1 $theurl temp-linecount-repo &&
  printf "('temp-linecount-repo' will be deleted automatically)\n\n\n" &&
  cloc temp-linecount-repo &&
  rm -rf temp-linecount-repo
