# BMA
mkdir build

cd build

cmake ..

cd tutorials/pathtracer

make -j9

cd ../../

./embree_pathtracer
