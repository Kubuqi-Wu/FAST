% Build all mex files

cd utils/mergeCells 
mex mergeCells.c
cd ../intersect
mex plusScalarFast.c
cd ../forwardMapping
mex forwardMapping.c
cd ../../