#!/bin/bash

install_dir=/usr/bin
cd $install_dir

git clone --branch wsl2 https://github.com/pbs-assess/gfiscam

cd $install_dir/gfiscam

sed -i 's#/home/cgrandin#/usr/bin#g' src/main/iscam.tpl

make dist
