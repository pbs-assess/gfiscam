#!/bin/bash

install_dir=/usr/bin
iscam_dir=$install_dir/gfiscam
cd $install_dir

rm -rf $iscam_dir
git clone --branch wsl2 https://github.com/pbs-assess/gfiscam

cd $install_dir/gfiscam

sed -i 's#/home/cgrandin#/usr/bin#g' src/main/iscam.tpl

make dist
