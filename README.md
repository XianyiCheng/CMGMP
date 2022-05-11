
install 
```   
sudo apt-get install build-essential cmake pkg-config git

sudo apt-get install libglpk-dev libeigen3-dev libassimp-dev libccd-dev libfcl-dev libboost-regex-dev libboost-system-dev libopenscenegraph-dev libbullet-dev libtinyxml2-dev liburdfdom-dev libxi-dev libxmu-dev freeglut3-dev
    
```
install qhull

```
git clone https://github.com/qhull/qhull.git
cd qhull
cd build
cmake ..
make
ctest
sudo make install
```

Update the external libraries
```
git submodule update --init --recursive --progress
```

Use cmake to build this project in the project root folder

```
mkdir build
cd build
cmake ..
make
```

Check out examples in the build folder

```
./examples/[name of the example]
```
