#------------------------------------------------------------------------------
cd ${0%/*} || exit 1
clear
./Allclean
##  MESHING ##
#------------------------------------------------------------------------------	
	echo "BlockMesh"
	blockMesh >> logMeshing
#------------------------------------------------------------------------------
	echo "surfaceFeatureExtract"
	surfaceFeatureExtract >> logMeshing
#------------------------------------------------------------------------------
	echo "snappyHexMesh"
	snappyHexMesh -overwrite >> logMeshing
#------------------------------------------------------------------------------
	echo "extrudeMesh"
	extrudeMesh >> logMeshing
#------------------------------------------------------------------------------
	echo "createPatch"	
	createPatch -overwrite  >> logMeshing
#------------------------------------------------------------------------------
	changeDictionary  >> logMeshing
	
##  SIMULATION ##
#------------------------------------------------------------------------------
	echo "DecomposePar"
	decomposePar
        echo "PimpleFoam"
	mpirun -np 3 pimpleFoam -parallel >> logSimulation
	echo "ReconstructPar"
	reconstructPar
	cd postProcessing/forces1/0
	mv forces.dat forces1.dat
	cd ../../forces2/0
	mv forces.dat forces2.dat
	cd ../../forces3/0
	mv forces.dat forces3.dat
#------------------------------------------------------------------------------
##  OVER  ##
#------------------------------------------------------------------------------
