#!/bin/bash

# Must be root to install iscam, or call from a Dockerfile because
# they install as root

install_dir=/usr/bin
iscam_dir=$install_dir/gfiscam
cd $install_dir

rm -rf $iscam_dir
git clone --branch wsl2 https://github.com/pbs-assess/gfiscam

cd $install_dir/gfiscam

# Change include statement
sed -i 's#/home/cgrandin#/usr/bin#g' src/main/iscam.tpl

make dist

cd $install_dir/gfiscam
cp docker/i $iscam_dir
chmode 777 i
sed -i 's#$ISCAM_BASE#/usr/bin/gfiscam#g' $iscam_dir/i
# Needed to make non-root users (who are in /etc/sudoers with full
# root privileges) access to create directories during the build.
# Look up "NFS root squash".
# https://superuser.com/questions/387111/sudo-make-install-permission-denied
chmod -R o+rw .

