/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "SVMMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include <fstream>
using namespace std;

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SVMMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        SVMMotionSolver,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * /
Foam::label Foam::SVMMotionSolver::samePoint(const point& a, const point& b) const
{
	scalar x = a.x()-b.x();
	if (x > 1e-8 || x < -1e-8) return 0;

	scalar y = a.y()-b.y();
	if (y > 1e-8 || y < -1e-8) return 0;

	scalar z = a.z()-b.z();	
	if (z > 1e-8 || z < -1e-8) return 0;

	return 1;
}

void Foam::SVMMotionSolver::makeParallelControlIDs()
{
    // Points that are neither on moving nor on static patches
    // will be marked with 0
    labelList markedPoints(mesh().nPoints(), 0);

	// Mark all points on coupled patches with -2
	// when they appear in moving patches, remark these points with -3
	// when they appear in static patches, remark these points with -4
	const polyBoundaryMesh& bMesh = mesh().boundaryMesh();
	forAll (bMesh, patchI)
	{
		const polyPatch& pp = bMesh[patchI];
		if (pp.coupled())
		{
			const labelList& mp = pp.meshPoints();
			forAll (mp, i)
			{
				markedPoints[mp[i]] = -2;
			}
		}
	}    

    // Mark all points on moving patches with 1
    label nMovingPoints = 0;

    forAll (movingPatches_, patchI)
    {
        // Find the patch in boundary
        label patchIndex =
            mesh().boundaryMesh().findPatchID(movingPatches_[patchI]);

        if (patchIndex < 0)
        {
	        continue;
	        //commented by gaoxiang
	        //while in parallel, some boundary patches may not in this region	        
            /*FatalErrorIn("void RBFMotionSolver::makeControlIDs()")
                << "Patch " << movingPatches_[patchI] << " not found.  "
                << "valid patch names: " << mesh().boundaryMesh().names()
                << abort(FatalError);
            */
        }

        const labelList& mp = mesh().boundaryMesh()[patchIndex].meshPoints();

        forAll (mp, i)
        {
	        if (markedPoints[mp[i]] == 0)
	        {
            	markedPoints[mp[i]] = 1;
            	nMovingPoints++;
        	}
        	else if (markedPoints[mp[i]] == -2)
        	{
	        	markedPoints[mp[i]] = -3;
	        	nMovingPoints++;
        	}
        }
    }

    // Mark moving points and select control points from moving patches
    movingIDs_.setSize(nMovingPoints);    

    const pointField& points = mesh().points();

    // Re-use counter to count moving points
    // Note: the control points also hold static points in the second part
    // of the list if static patches are included in the RBF
    // HJ, 24/Mar/2011
    nMovingPoints = 0;

    // Count moving points first
    forAll (markedPoints, i)
    {
        if (markedPoints[i] == 1)
        {
            movingIDs_[nMovingPoints] = i;
            nMovingPoints++;
        }
    }
    
    label nLocalMovingPoints = nMovingPoints;
    forAll (markedPoints, i)
    {
        if (markedPoints[i] == -3)
        {
            movingIDs_[nMovingPoints] = i;
            nMovingPoints++;
        }
    }

    movingIDs_.setSize(nMovingPoints);

	Info<< "Total points on moving boundaries: " << nMovingPoints << endl;

    // Actual location of moving points will be set later on request
    // HJ, 19/Dec/2008
    movingPoints_.setSize(nMovingPoints, vector::zero);

    // Mark all points on static patches with -1
    label nStaticPoints = 0;

    forAll (staticPatches_, patchI)
    {
        // Find the patch in boundary
        label patchIndex =
            mesh().boundaryMesh().findPatchID(staticPatches_[patchI]);

        if (patchIndex < 0)
        {
	        continue;
	        //commented by gaoxiang
	        //while in parallel, some boundary patches may not in this region
            /*FatalErrorIn("void RBFMotionSolver::makeControlPoints()")
                << "Patch " << staticPatches_[patchI] << " not found.  "
                << "valid patch names: " << mesh().boundaryMesh().names()
                << abort(FatalError);
            */
        }

        const labelList& mp = mesh().boundaryMesh()[patchIndex].meshPoints();

        forAll (mp, i)
        {
            if (markedPoints[mp[i]] == 0)
	        {
            	markedPoints[mp[i]] = -1;
            	nStaticPoints++;
        	}
        	else if (markedPoints[mp[i]] == -2)
        	{
	        	markedPoints[mp[i]] = -4;
	        	nStaticPoints++;
        	}
        }
    }
    
    staticIDs_.setSize(nStaticPoints);

    // Re-use counter
    nStaticPoints = 0;

    // Count total number of control points
    forAll (markedPoints, i)
    {
        if (markedPoints[i] == -1)
        {
            staticIDs_[nStaticPoints] = i;
            nStaticPoints++;
        }
    }
    
    label nLocalStaticPoints = nStaticPoints;
    forAll (markedPoints, i)
    {
        if (markedPoints[i] == -4)
        {
            staticIDs_[nStaticPoints] = i;
            nStaticPoints++;
        }
    }

    staticIDs_.setSize(nStaticPoints);

    Info<< "Total points on static boundaries: " << nStaticPoints << endl;

    if (includeStaticPatches_)
    {
	    Info << "include static patches" << endl;
    }
    else
    {
	    Info << "not include static patches" << endl;
    }

    // Control IDs also potentially include points on static patches
    // controlIDs_ in fact does not used in parallel running
    controlIDs_.setSize(movingIDs_.size() + staticIDs_.size());
	label nLocalRealBoundaryPoints = 0;
    forAll(movingIDs_, i)
    {
	    controlIDs_[nLocalRealBoundaryPoints] = movingIDs_[i];
	    nLocalRealBoundaryPoints++;
    }

    if (includeStaticPatches_)
    {
	    forAll(staticIDs_, i)
	    {
		    controlIDs_[nLocalRealBoundaryPoints] = staticIDs_[i];
		    nLocalRealBoundaryPoints++;
	    }
    }
    controlIDs_.setSize(nLocalRealBoundaryPoints);

    Info<< "Selected " << nLocalRealBoundaryPoints
        << " total Local real boundary control points" << endl;


    // Pick up all internal points, include processor boundary patches
    //control points + internal points = all points in this domain
    internalIDs_.setSize(points.size());
    internalPoints_.setSize(points.size());

    // Count internal points
    label nInternalPoints = 0;

    forAll (markedPoints, i)
    {
        if (markedPoints[i] == 0 || markedPoints[i] == -2)
        {
            // Grab internal point
            internalIDs_[nInternalPoints] = i;
            internalPoints_[nInternalPoints] = points[i];
            nInternalPoints++;
        }
        //TODO sperate coupled patch points
    }

    Info << "Number of internal points: " << nInternalPoints << endl;

    /*if (Pstream::myProcNo() == 0)
    {
        std::cout << "markedPoints[3155]=" << markedPoints[3155] << std::endl;
    }
    else
    {
        std::cout << "markedPoints[4599]=" << markedPoints[4599] << std::endl;
    }*/

    // Resize the lists
    internalIDs_.setSize(nInternalPoints);
    internalPoints_.setSize(nInternalPoints);

    newInternalPoints_.setSize(nInternalPoints);

    // get global boundary infomation
    //moving points
    List<vectorField> globalMovingPoints;
    globalMovingPoints.setSize(Pstream::nProcs());
    vectorField& myMovingPoints = globalMovingPoints[Pstream::myProcNo()];

    myMovingPoints.setSize(nLocalMovingPoints);
    for (int i = 0; i < nLocalMovingPoints; i++)
    {
	    myMovingPoints[i] = points[movingIDs_[i]];
    }    
    
    Pstream::gatherList(globalMovingPoints);
    Pstream::scatterList(globalMovingPoints);
    
	//static points
    List<vectorField> globalStaticPoints;
    globalStaticPoints.setSize(Pstream::nProcs());
    vectorField& myStaticPoints = globalStaticPoints[Pstream::myProcNo()];

    myStaticPoints.setSize(nLocalStaticPoints);
    for (int i = 0; i < nLocalStaticPoints; i++)
    {
	    myStaticPoints[i] = points[staticIDs_[i]];
    }    
    
    Pstream::gatherList(globalStaticPoints);
    Pstream::scatterList(globalStaticPoints);
    
	//shared points
	List<vectorField> sharedMovingPoints;
    sharedMovingPoints.setSize(Pstream::nProcs());
    vectorField& mySharedMovingPoints = sharedMovingPoints[Pstream::myProcNo()];

    mySharedMovingPoints.setSize(nMovingPoints-nLocalMovingPoints);
    for (int i = nLocalMovingPoints, j = 0; i < nMovingPoints; i++, j++)
    {
	    mySharedMovingPoints[j] = points[movingIDs_[i]];
    }    
    
    Pstream::gatherList(sharedMovingPoints);
    Pstream::scatterList(sharedMovingPoints);

	// sharedStaticPoints
    List<vectorField> sharedStaticPoints;
    sharedStaticPoints.setSize(Pstream::nProcs());
    vectorField& mySharedStaticPoints = sharedStaticPoints[Pstream::myProcNo()];

    mySharedStaticPoints.setSize(nStaticPoints-nLocalStaticPoints);
    for (int i = nLocalStaticPoints, j = 0; i < nStaticPoints; i++, j++)
    {
	    mySharedStaticPoints[j] = points[staticIDs_[i]];
    }    
    
    Pstream::gatherList(sharedStaticPoints);
    Pstream::scatterList(sharedStaticPoints);

    //pick out reduplicate point
	//a point in sharedMovingPoints may also appear in sharedStaticPoints on different processors
    vectorField finalSharedMovingPoints;
    vectorField finalSharedStaticPoints;
    label nFinalSharedMovingPoints = 0;
    label nFinalSharedStaticPoints = 0;
    for (int i = 0; i < Pstream::nProcs(); i++)
    {
	    nFinalSharedMovingPoints += sharedMovingPoints[i].size();
	    nFinalSharedStaticPoints += sharedStaticPoints[i].size();
    }
    finalSharedMovingPoints.setSize(nFinalSharedMovingPoints);
    finalSharedStaticPoints.setSize(nFinalSharedStaticPoints);

    //insure moving point saved first
    nFinalSharedMovingPoints = nFinalSharedStaticPoints = 0;

    for (int i = 0; i < Pstream::nProcs(); i++)
    {  
	    const vectorField& iMovingPoints = sharedMovingPoints[i];
	    label tmpMovingLen = nFinalSharedMovingPoints;
	    
	    forAll (iMovingPoints, j)
	    {
		    label flag = 1;
		    for (int k = 0; k < tmpMovingLen; k++)
		    {
			    if (samePoint(iMovingPoints[j], finalSharedMovingPoints[k]))
			    {
				    flag = 0;
				    break;
			    }
		    }
		    if (flag)
		    {
			    finalSharedMovingPoints[nFinalSharedMovingPoints] = iMovingPoints[j];
			    nFinalSharedMovingPoints++;
		    }
	    }
    }
    
    for (int i = 0; i < Pstream::nProcs(); i++)
    {	    
	    const vectorField& iStaticPoints = sharedStaticPoints[i];
	    label tmpStaticLen = nFinalSharedStaticPoints;
	    	    	    
	    forAll (iStaticPoints, j)
	    {
		    label flag = 1;
		    for (int k = 0; k < tmpStaticLen; k++)
		    {
			    if (samePoint(iStaticPoints[j], finalSharedStaticPoints[k]))
			    {
				    flag = 0;
				    break;
			    }
		    }
		    if (flag)
		    {
			    for (int k = 0; k < nFinalSharedMovingPoints; k++)
			    {
				    if (samePoint(iStaticPoints[j], finalSharedMovingPoints[k]))
				    {
					    flag = 0;
					    break;
				    }
			    }
		    }
		    if (flag)
		    {
			    finalSharedStaticPoints[nFinalSharedStaticPoints] = iStaticPoints[j];
			    nFinalSharedStaticPoints++;
		    }
	    }
    }
    finalSharedMovingPoints.setSize(nFinalSharedMovingPoints);
    finalSharedStaticPoints.setSize(nFinalSharedStaticPoints);
	
	label nRealControlPoints = nFinalSharedMovingPoints + nFinalSharedStaticPoints;
	forAll(globalMovingPoints, i)
	{
		nRealControlPoints += globalMovingPoints[i].size() + globalStaticPoints[i].size();
	}
	
	std::cout << "processor " << Pstream::myProcNo() << ", total boundary points=" << nRealControlPoints << std::endl;

    controlPoints_.setSize(nRealControlPoints);
    nRealControlPoints = 0;

	/*              Data layout schematic
	  |----------------controlPoints_---------------------|
	  |---------motion_(Dynamic)-----|---Static-----------|
	*/   
    for (int i = 0; i < Pstream::nProcs(); i++)
    {	    
	    const vectorField& tmpPoints = globalMovingPoints[i];
	    forAll(tmpPoints, j)
	    {
		    /*if (tmpPoints[j].x() > minX && tmpPoints[j].x() < maxX &&
		        tmpPoints[j].y() > minY && tmpPoints[j].y() < maxY &&
		        tmpPoints[j].z() > minZ && tmpPoints[j].z() < maxZ)*/
		    {
		        controlPoints_[nRealControlPoints] = tmpPoints[j];
		        nRealControlPoints++;
		    }
	    }
	}

	for (int i = 0; i < nFinalSharedMovingPoints; i++)
	{
		/*if (finalSharedMovingPoints[i].x() > minX && finalSharedMovingPoints[i].x() < maxX &&
		    finalSharedMovingPoints[i].y() > minY && finalSharedMovingPoints[i].y() < maxY &&
		    finalSharedMovingPoints[i].z() > minZ && finalSharedMovingPoints[i].z() < maxZ)*/
		{
	        controlPoints_[nRealControlPoints] = finalSharedMovingPoints[i];
	        nRealControlPoints++;
		}
	}
	nRealMovingPoints_ = nRealControlPoints; //moving part of the entire boundary point

	if (includeStaticPatches_)
	{		
		for (int i = 0; i < Pstream::nProcs(); i++)
		{
		    const vectorField& tmpPoints = globalStaticPoints[i];
		    forAll(tmpPoints, j)
		    {
			    /*if (tmpPoints[j].x() > minX && tmpPoints[j].x() < maxX &&
			        tmpPoints[j].y() > minY && tmpPoints[j].y() < maxY &&
			        tmpPoints[j].z() > minZ && tmpPoints[j].z() < maxZ)*/
			    {
			        controlPoints_[nRealControlPoints] = tmpPoints[j];
			        nRealControlPoints++;
			    }
		    }
	    }
	    
		for (int i = 0; i < nFinalSharedStaticPoints; i++)
		{
			/*if (finalSharedStaticPoints[i].x() > minX && finalSharedStaticPoints[i].x() < maxX &&
			    finalSharedStaticPoints[i].y() > minY && finalSharedStaticPoints[i].y() < maxY &&
			    finalSharedStaticPoints[i].z() > minZ && finalSharedStaticPoints[i].z() < maxZ)*/
			{
		        controlPoints_[nRealControlPoints] = finalSharedStaticPoints[i];
		        nRealControlPoints++;
			}
		}
	}

    //Resize
    controlPoints_.setSize(nRealControlPoints);
    std::cout << "processor " << Pstream::myProcNo() << ", nRealControlPoints=" << nRealControlPoints << std::endl;
    

    motion_.setSize(nRealControlPoints, vector::zero);
    
    xSpace = Malloc(struct svm_node, nRealControlPoints*4);

    label j = 0;
    forAll (controlPoints_, i)
    {        
        xSpace[j  ].index = 1;
        xSpace[j  ].value = controlPoints_[i].x();
        xSpace[j+1].index = 2;
        xSpace[j+1].value = controlPoints_[i].y();
        xSpace[j+2].index = 3;
        xSpace[j+2].value = controlPoints_[i].z();
        xSpace[j+3].index = -1;
        j += 4;
    }

    // SVM Problem set
    xSVMProblem.l = nRealControlPoints;
    xSVMProblem.y = Malloc(double, xSVMProblem.l);
    xSVMProblem.x = Malloc(struct svm_node *, xSVMProblem.l);
    
    ySVMProblem.l = nRealControlPoints;
    ySVMProblem.y = Malloc(double, ySVMProblem.l);
    ySVMProblem.x = Malloc(struct svm_node *, ySVMProblem.l);
    
    zSVMProblem.l = nRealControlPoints;
    zSVMProblem.y = Malloc(double, zSVMProblem.l);
    zSVMProblem.x = Malloc(struct svm_node *, zSVMProblem.l);

	j = 0;
    for (label i = 0; i < nRealControlPoints; i++)
    {
	    xSVMProblem.x[i] = &xSpace[j];
	    ySVMProblem.x[i] = &xSpace[j];
	    zSVMProblem.x[i] = &xSpace[j];
	    j += 4;
    }
    
    // initial moving points' coordinates
    statPoints_.setSize(nRealMovingPoints_, vector::zero);
    forAll (statPoints_, i)
    {
        statPoints_[i] = controlPoints_[i];
    }

    localMotion_.setSize(movingIDs_.size(), vector::zero);
    localStatPoints_.setSize(movingIDs_.size(), vector::zero);
    forAll (movingIDs_, i)
    {
	    localStatPoints_[i] = points[movingIDs_[i]];
    }

    // Inital SVM
    xSVMModel = NULL;
    ySVMModel = NULL,
    zSVMModel = NULL;

    param.svm_type = svm_type_; //3=EPSILON_SVR, 4=NU_SVR 
    if (svm_type_ == EPSILON_SVR)
    {
        Info << "Choose EPSILON_SVR" << endl;
    }
    else if (svm_type_ == NU_SVR)
    {
        Info << "Choose NU_SVR" << endl;
    }
    param.kernel_type = kernel_type_; //2=RBF
    param.degree = degree_;  // for polynomial kernel function
    param.gamma = gamma_;  // 1/num_features, for RBF/sigmoid function
    param.coef0 = coef0_;  // for polynomial/sigmoid function
    param.nu = nu_;        //for nu-SVR
    param.cache_size = cache_size_;
    param.C = C_;
    param.eps = eps_;  // stop condition
    param.p = p_; // for espilon_SVR
    param.shrinking = shrinking_;  //use heuristic method or not
    param.probability = probability_;
    param.nr_weight = 0;
    param.weight_label = NULL;
    param.weight = NULL;


	//use 2 processors for test
    /*std::fstream fp;
    if (Pstream::myProcNo() == 0)
    {
    	fp.open("p0.txt", ios::out);
	}
	else
	{
		fp.open("p1.txt", ios::out);
	}
	fp << controlPoints_.size() << std::endl;
    forAll (controlPoints_, i)
    {
	    fp << controlPoints_[i][0] << " " << controlPoints_[i][1] << " " << controlPoints_[i][2] << std::endl;
    }
    fp.close();*/
}

void Foam::SVMMotionSolver::makeControlIDs()
{
    // Points that are neither on moving nor on static patches
    // will be marked with 0
    labelList markedPoints(mesh().nPoints(), 0);

    // Mark all points on moving patches with 1
    label nMovingPoints = 0;

    forAll (movingPatches_, patchI)
    {
        // Find the patch in boundary
        label patchIndex =
            mesh().boundaryMesh().findPatchID(movingPatches_[patchI]);

        if (patchIndex < 0)
        {
            FatalErrorIn("void SVMMotionSolver::makeControlIDs()")
                << "Patch " << movingPatches_[patchI] << " not found.  "
                << "valid patch names: " << mesh().boundaryMesh().names()
                << abort(FatalError);
        }

        const labelList& mp = mesh().boundaryMesh()[patchIndex].meshPoints();

        forAll (mp, i)
        {
            markedPoints[mp[i]] = 1;
            nMovingPoints++;
        }
    }

    // Mark moving points and select control points from moving patches
    movingIDs_.setSize(nMovingPoints);

    Info<< "Total points on moving boundaries: " << nMovingPoints << endl;

    const vectorField& points = mesh().points();

    // Re-use counter to count moving points
    // Note: the control points also hold static points in the second part
    // of the list if static patches are included in the SVM
    nMovingPoints = 0;

    // Count moving points first
    forAll (markedPoints, i)
    {
        if (markedPoints[i] == 1)
        {
            // Grab internal point
            movingIDs_[nMovingPoints] = i;
            nMovingPoints++;
        }
    }

    movingIDs_.setSize(nMovingPoints);

    // Actual location of moving points will be set later on request
    movingPoints_.setSize(nMovingPoints, vector::zero);

    // initial moving points' coordinates
    statPoints_.setSize(nMovingPoints, vector::zero);
    forAll (movingIDs_, i)
    {
        statPoints_[i] = points[movingIDs_[i]];
    }

    // Mark all points on static patches with -1
    label nStaticPoints = 0;

    forAll (staticPatches_, patchI)
    {
        // Find the patch in boundary
        label patchIndex =
            mesh().boundaryMesh().findPatchID(staticPatches_[patchI]);

        if (patchIndex < 0)
        {
            FatalErrorIn("void SVMMotionSolver::makeControlPoints()")
                << "Patch " << staticPatches_[patchI] << " not found.  "
                << "valid patch names: " << mesh().boundaryMesh().names()
                << abort(FatalError);
        }

        const labelList& mp = mesh().boundaryMesh()[patchIndex].meshPoints();

        forAll (mp, i)
        {
	        if (markedPoints[mp[i]] == 0)
	        {
            	markedPoints[mp[i]] = -1;
            	nStaticPoints++;
        	}
        }
    }

    Info<< "Total points on static boundaries: " << nStaticPoints << endl;
    staticIDs_.setSize(nStaticPoints);

    // Re-use counter
    nStaticPoints = 0;

    // Count total number of control points
    forAll (markedPoints, i)
    {
        if (markedPoints[i] == -1)
        {
            staticIDs_[nStaticPoints] = i;
            nStaticPoints++;
        }
    }

    staticIDs_.setSize(nStaticPoints);

    // Control IDs also potentially include points on static patches
    controlIDs_.setSize(movingIDs_.size() + staticIDs_.size());
    
    label nControlPoints = 0;

    forAll (movingIDs_, i)
    {
	    // Pick point as control point
        controlIDs_[nControlPoints] = movingIDs_[i];
        nControlPoints++;
    }

    Info<< "Selected " << nControlPoints
        << " control points on moving boundaries" << endl;

    if (includeStaticPatches_)
    {
	    forAll (staticIDs_, i)
    	{
	    	// Pick point as control point
        	controlIDs_[nControlPoints] = staticIDs_[i];
        	nControlPoints++;
    	}

        Info<< "Selected " << nControlPoints
            << " total control points" << endl;
    }

    // Resize control IDs
    controlIDs_.setSize(nControlPoints);

    // Pick up point locations
    controlPoints_.setSize(nControlPoints);
    motion_.setSize(nControlPoints, vector::zero);
    xSpace = Malloc(struct svm_node, nControlPoints*4);

    label j = 0;
    forAll (controlIDs_, i)
    {
        controlPoints_[i] = points[controlIDs_[i]];
        xSpace[j  ].index = 1;
        xSpace[j  ].value = controlPoints_[i].x();
        xSpace[j+1].index = 2;
        xSpace[j+1].value = controlPoints_[i].y();
        xSpace[j+2].index = 3;
        xSpace[j+2].value = controlPoints_[i].z();
        xSpace[j+3].index = -1;
        j += 4;
    }

    // SVM Problem set
    xSVMProblem.l = nControlPoints;
    xSVMProblem.y = Malloc(double, xSVMProblem.l);
    xSVMProblem.x = Malloc(struct svm_node *, xSVMProblem.l);
    
    ySVMProblem.l = nControlPoints;
    ySVMProblem.y = Malloc(double, ySVMProblem.l);
    ySVMProblem.x = Malloc(struct svm_node *, ySVMProblem.l);
    
    zSVMProblem.l = nControlPoints;
    zSVMProblem.y = Malloc(double, zSVMProblem.l);
    zSVMProblem.x = Malloc(struct svm_node *, zSVMProblem.l);

	j = 0;
    for (label i = 0; i < nControlPoints; i++)
    {
	    xSVMProblem.x[i] = &xSpace[j];
	    ySVMProblem.x[i] = &xSpace[j];
	    zSVMProblem.x[i] = &xSpace[j];
	    j += 4;
    }
    
    // Pick up all internal points
    internalIDs_.setSize(points.size());
    internalPoints_.setSize(points.size());

    // Count internal points
    label nInternalPoints = 0;

    forAll (markedPoints, i)
    {
        if (markedPoints[i] == 0)
        {
            // Grab internal point
            internalIDs_[nInternalPoints] = i;
            internalPoints_[nInternalPoints] = points[i];
            nInternalPoints++;
        }
    }

    Info << "Number of internal points: " << nInternalPoints << endl;

    // Resize the lists
    internalIDs_.setSize(nInternalPoints);
    internalPoints_.setSize(nInternalPoints);

    newInternalPoints_.setSize(nInternalPoints);

    // Inital SVM
    xSVMModel = NULL;
    ySVMModel = NULL,
    zSVMModel = NULL;

    param.svm_type = svm_type_; //3=EPSILON_SVR, 4=NU_SVR 
    if (svm_type_ == EPSILON_SVR)
    {
        Info << "Choose EPSILON_SVR" << endl;
    }
    else if (svm_type_ == NU_SVR)
    {
        Info << "Choose NU_SVR" << endl;
    }
    param.kernel_type = kernel_type_; //2=RBF
    param.degree = degree_;  // for polynomial kernel function
    param.gamma = gamma_;  // 1/num_features, for RBF/sigmoid function
    param.coef0 = coef0_;  // for polynomial/sigmoid function
    param.nu = nu_;        //for nu-SVR
    param.cache_size = cache_size_;
    param.C = C_;
    param.eps = eps_;  // stop condition
    param.p = p_; // for espilon_SVR
    param.shrinking = shrinking_;  //use heuristic method or not
    param.probability = probability_;
    param.nr_weight = 0;
    param.weight_label = NULL;
    param.weight = NULL;
}


void Foam::SVMMotionSolver::setMovingPoints() const
{
    const vectorField& points = mesh().points();

    // Set moving points
    forAll (movingIDs_, i)
    {
        movingPoints_[i] = points[movingIDs_[i]];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SVMMotionSolver::SVMMotionSolver
(
    const polyMesh& mesh,
    Istream&
)
:
    motionSolver(mesh),
    movingPatches_(lookup("movingPatches")),
    staticPatches_(lookup("staticPatches")),
    includeStaticPatches_(lookup("includeStaticPatches")),
    frozenRegression_(lookup("frozenRegression")),
    movingIDs_(0),
    movingPoints_(0),
    staticIDs_(0),
    controlIDs_(0),
    controlPoints_(0),
    internalIDs_(0),
    internalPoints_(0),
    motion_(0),
    rotationAmplitude_(readScalar(lookup("rotationAmplitude"))),
    rotationFrequency_(readScalar(lookup("rotationFrequency"))),
    translationAmplitude_(lookup("translationAmplitude")),
    translationFrequency_(lookup("translationFrequency")),
    initialRotationOrigin_(lookup("initialRotationOrigin")),
    statPoints_(0),
    newInternalPoints_(0),
    sumCtrlPoints2D_(0),
    sumCtrlPoints3D_(0),

    svm_type_(readScalar(lookup("svm_type"))),
    kernel_type_(readScalar(lookup("kernel_type"))),
    degree_(readScalar(lookup("degree"))),  // for polynomial kernel function
    gamma_(readScalar(lookup("gamma"))),  // 1/num_features, for RBF/sigmoid function
    coef0_(readScalar(lookup("coef0"))),  // for polynomial/sigmoid function
    nu_(readScalar(lookup("nu"))),
    cache_size_(readScalar(lookup("cache_size"))),
    C_(readScalar(lookup("C"))),
    eps_(readScalar(lookup("eps"))),  // stop condition
    p_(readScalar(lookup("p"))), // for espilon_svr
    shrinking_(readScalar(lookup("shrinking"))),  //use heuristic method or not
    probability_(0)
{
	if (Pstream::parRun())
	{
		makeParallelControlIDs();
	}
	else
	{
		makeControlIDs();
	}    
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::SVMMotionSolver::~SVMMotionSolver()
{
	svm_destroy_param(&param);
	free(xSVMProblem.x);
	free(xSVMProblem.y);
	free(ySVMProblem.x);
	free(ySVMProblem.y);
	free(zSVMProblem.x);
	free(zSVMProblem.y);
	free(xSpace);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::SVMMotionSolver::setMotion(const vectorField& m)
{
    if (Pstream::parRun())
	{
		if (m.size() != nRealMovingPoints_)
	    {
	        FatalErrorIn
	        (
	            "void RBFMotionSolver::setMotion(const vectorField& m)"
	        )   << "Incorrect size of motion points: m = " << m.size()
	            << " nRealMovingPoints_ = " << nRealMovingPoints_
	            << abort(FatalError);
	    }
	}
	else
	{
	    if (m.size() != movingIDs_.size())
	    {
	        FatalErrorIn
	        (
	            "void RBFMotionSolver::setMotion(const vectorField& m)"
	        )   << "Incorrect size of motion points: m = " << m.size()
	            << " movingIDs = " << movingIDs_.size()
	            << abort(FatalError);
	    }
	}

    // Motion of static points is zero and moving points are first
    // in the list.

    if (!frozenRegression_)
    {
        // Set control points
        const vectorField& points = mesh().points();

        forAll(internalIDs_, i)
        {
	        internalPoints_[i] = points[internalIDs_[i]];
        }
        if (Pstream::parRun())
        {
	        label j = 0;
	        for (int i = 0; i < nRealMovingPoints_; i++)
	        {
		        controlPoints_[i] += motion_[i];
               	xSpace[j  ].value = controlPoints_[i].x();
            	xSpace[j+1].value = controlPoints_[i].y();
            	xSpace[j+2].value = controlPoints_[i].z();
            	j += 4;
	        }
        }
        else
        {
            label j = 0;
        	forAll (controlIDs_, i)
        	{
            	controlPoints_[i] = points[controlIDs_[i]];
            	xSpace[j  ].value = controlPoints_[i].x();
            	xSpace[j+1].value = controlPoints_[i].y();
            	xSpace[j+2].value = controlPoints_[i].z();
            	j += 4;
        	}
        }
    }
    
    forAll (m, i)
    {
        motion_[i] = m[i];
    }    

    forAll (motion_, i)
    {
        xSVMProblem.y[i] = motion_[i].x();
        ySVMProblem.y[i] = motion_[i].y();
        zSVMProblem.y[i] = motion_[i].z();
    }    
}


const Foam::vectorField& Foam::SVMMotionSolver::movingPoints() const
{
    // Update moving points based on current mesh
    setMovingPoints();

    return movingPoints_;
}


Foam::tmp<Foam::vectorField> Foam::SVMMotionSolver::curPoints() const
{
    // Prepare new points: same as old point
    tmp<vectorField> tcurPoints
    (
        new vectorField(mesh().nPoints(), vector::zero)
    );
    vectorField& curPoints = tcurPoints();

    // Add motion to existing points

    // 1. Insert prescribed motion of moving points
    if (Pstream::parRun())
    {
	    forAll (movingIDs_, i)
	    {
	        curPoints[movingIDs_[i]] = localMotion_[i]; //localMotion_ is the moving points of this domain
	    }
    }
    else
    {
	    forAll (movingIDs_, i)
	    {
	        curPoints[movingIDs_[i]] = motion_[i];
	    }
    }
    
    // 2. Insert zero motion of static points
    /*forAll (staticIDs_, i)
    {
        curPoints[staticIDs_[i]] = vector::zero;
    }*/

    // 3. Insert SVM regression motion
    forAll (internalIDs_, i)
    {
        curPoints[internalIDs_[i]] = newInternalPoints_[i];
    }

    // 4. Add old point positions
    curPoints += mesh().points();

    twoDCorrectPoints(tcurPoints());

    return tcurPoints;
}


void Foam::SVMMotionSolver::solve()
{
	// Move the moving boundary
#   include "kinematicModel.H"
//#   include "kinematicModel_for_M6.H"
    setMotion(motionBoundary);

    // Call interpolation
    SVMRegression();
}


void Foam::SVMMotionSolver::updateMesh(const mapPolyMesh&)
{
    // Recalculate control point IDs
    makeControlIDs();
}

void Foam::SVMMotionSolver::SVMRegression()
{
	struct svm_node *x;
	const char *error_msg;
	
	error_msg = svm_check_parameter(&xSVMProblem, &param);
	if(error_msg)
	{
        printf("ERROR: %s\n", error_msg);
		FatalErrorIn
    	(
        	"void Foam::SVMMotionSolver::SVMRegression()"
    	)   << abort(FatalError);
	}    
	xSVMModel = svm_train(&xSVMProblem, &param);
	
	error_msg = svm_check_parameter(&ySVMProblem, &param);
	if(error_msg)
	{
        printf("ERROR: %s\n", error_msg);
		FatalErrorIn
    	(
        	"void Foam::SVMMotionSolver::SVMRegression()"
    	)   << abort(FatalError);
	}
	ySVMModel = svm_train(&ySVMProblem, &param);
	
	error_msg = svm_check_parameter(&zSVMProblem, &param);
	if(error_msg)
	{
        printf("ERROR: %s\n", error_msg);
		FatalErrorIn
    	(
        	"void Foam::SVMMotionSolver::SVMRegression()"
    	)   << abort(FatalError);
	}
	zSVMModel = svm_train(&zSVMProblem, &param);

	x = Malloc(struct svm_node, 4);
	forAll (internalPoints_, i)
	{
		x[0].index = 1;
		x[0].value = internalPoints_[i].x();
		x[1].index = 2;
		x[1].value = internalPoints_[i].y();
		x[2].index = 3;
		x[2].value = internalPoints_[i].z();
		x[3].index = -1;
		
		double xPredict = svm_predict(xSVMModel, x);
		double yPredict = svm_predict(ySVMModel, x);
		double zPredict = svm_predict(zSVMModel, x);
		Vector<double> tmpPoint(xPredict, yPredict, zPredict);
		newInternalPoints_[i] = tmpPoint;
	}
    
    sumCtrlPoints2D_ += xSVMModel->l + ySVMModel->l;
    sumCtrlPoints3D_ += xSVMModel->l + ySVMModel->l + zSVMModel->l;

    Info << "this step 2D average controlPoints =" << (xSVMModel->l + ySVMModel->l)/2 << endl;
    Info << "this step 3D average controlPoints =" << (xSVMModel->l + ySVMModel->l + zSVMModel->l)/3 << endl;

    Info << "2D sum control points of all steps = " << sumCtrlPoints2D_ 
         << ", average control points = " 
         << (sumCtrlPoints2D_/2.0)/mesh().time().timeIndex() << endl;

    Info << "3D sum control points of all steps = " << sumCtrlPoints3D_ 
         << ", average control points = " 
         << (sumCtrlPoints3D_/3.0)/mesh().time().timeIndex() << endl;


	//print support vector
/*	fstream foutput;
    foutput.open("controlpointsX.plt", ios::out);    
    foutput << "variables = \"x\", \"y\", \"z\"" << std::endl;
    foutput << "zone n = " << xSVMModel->l << ", e = " << 1
          		<< ", f = fepoint, et = lineseg" << std::endl;
    for (int i = 0; i < xSVMModel->l; i++)
    {
	    struct svm_node *tmpSV = xSVMModel->SV[i];
	    foutput << tmpSV[0].value << " " << tmpSV[1].value << " " << tmpSV[2].value << std::endl;
    }
	foutput << 1 << " " << 1 << std::endl;
    foutput.close();

    foutput.open("controlpointsY.plt", ios::out);    
    foutput << "variables = \"x\", \"y\", \"z\"" << std::endl;
    foutput << "zone n = " << ySVMModel->l << ", e = " << 1
          		<< ", f = fepoint, et = lineseg" << std::endl;
    for (int i = 0; i < ySVMModel->l; i++)
    {
	    struct svm_node *tmpSV = ySVMModel->SV[i];
	    foutput << tmpSV[0].value << " " << tmpSV[1].value << " " << tmpSV[2].value << std::endl;
    }
	foutput << 1 << " " << 1 << std::endl;
    foutput.close();

    foutput.open("controlpointsZ.plt", ios::out);    
    foutput << "variables = \"x\", \"y\", \"z\"" << std::endl;
    foutput << "zone n = " << zSVMModel->l << ", e = " << 1
          		<< ", f = fepoint, et = lineseg" << std::endl;
    for (int i = 0; i < zSVMModel->l; i++)
    {
	    struct svm_node *tmpSV = zSVMModel->SV[i];
	    foutput << tmpSV[0].value << " " << tmpSV[1].value << " " << tmpSV[2].value << std::endl;
    }
	foutput << 1 << " " << 1 << std::endl;
    foutput.close();*/

    svm_free_and_destroy_model(&xSVMModel);
    svm_free_and_destroy_model(&ySVMModel);
    svm_free_and_destroy_model(&zSVMModel);
	free(x);
}

// ************************************************************************* //
