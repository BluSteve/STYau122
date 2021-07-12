#!/bin/bash
mkdir ~/.jdks
cd ~/Downloads || exit
wget "https://github.com/AdoptOpenJDK/openjdk11-binaries/releases/download/jdk-11.0.11%2B9/OpenJDK11U-jdk_x64_linux_hotspot_11.0.11_9.tar.gz"
tar -xzf OpenJDK11U-jdk_x64_linux_hotspot_11.0.11_9.tar.gz -C ~/.jdks
export PATH=~/.jdks/jdk-11.0.11+9/bin:$PATH
echo "export PATH=~/.jdks/jdk-11.0.11+9/bin:$PATH" >> ~/.bashrc
rm OpenJDK11U-jdk_x64_linux_hotspot_11.0.11_9.tar.gz

sudo apt install libgfortran5 htop
java -version