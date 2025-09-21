# BMA

This is the implementation code for the paper "BMA: Accelerating BVH Traversal in Ray Tracing via Probabilistic Sampling".

You can modify the file 'pathtrace_device.cpp' to adjust the convolution-like operation.

The following are the steps for building and executing:

mkdir build

cd build

cmake ..

cd tutorials/pathtracer

make -j9

cd ../../

./embree_pathtracer


