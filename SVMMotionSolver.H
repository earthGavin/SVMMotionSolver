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

Class
    SVMMotionSolver

Description
    Support Vector Machine regression motion solver

Author
    Gao Xiang, College of Computer, NUDT China.  All rights reserved.

SourceFiles
    SVMMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef SVMMotionSolver_H
#define SVMMotionSolver_H

#include "motionSolver.H"
#include "polyMesh.H"

#include "svm.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class SVMMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class SVMMotionSolver
:
    public motionSolver
{
    // Private data

        //- Moving patch names
        wordList movingPatches_;

        //- Static patch names
        wordList staticPatches_;

        //- Include zero motion of static patches in SVM regression
        Switch includeStaticPatches_;

        //- Frozen regression
        Switch frozenRegression_;

        //- Moving point IDs
        labelList movingIDs_;

        //- Moving points on the boundary
        mutable vectorField movingPoints_;

        //- Static point IDs
        labelList staticIDs_;

        //- Control point IDs
        labelList controlIDs_;

        //- Control points on the boundary
        vectorField controlPoints_;

        //- Internal point IDs
        labelList internalIDs_;

        //- Internal points
        vectorField internalPoints_;

        //- Prescribed motion on control points
        vectorField motion_;

        //- Prescribed motion parameters
        scalar rotationAmplitude_;
        scalar rotationFrequency_;
        vector translationAmplitude_;
        vector translationFrequency_;
        vector initialRotationOrigin_;

        //- Stationary points, i.e. inital moving points
        vectorField statPoints_;

        //- New internal points
        vectorField newInternalPoints_;

        //- parallel data        
        label nRealMovingPoints_;               
        pointField localStatPoints_;
        vectorField localMotion_;

        //- SVM training problem
        struct svm_problem xSVMProblem;
        struct svm_problem ySVMProblem;
        struct svm_problem zSVMProblem;

        //- SVM model
        struct svm_model *xSVMModel;
        struct svm_model *ySVMModel;
        struct svm_model *zSVMModel;

        //- SVM parameters
        struct svm_parameter param;
        label svm_type_;


        label kernel_type_;
        label degree_; /* for poly */
        scalar gamma_;   /* for poly/rbf/sigmoid */
        scalar coef0_;   /* for poly/sigmoid */

        /* these are for training only */
        scalar cache_size_; /* in MB */
        scalar eps_; /* stopping criteria */
        scalar C_;   /* for C_SVC, EPSILON_SVR and NU_SVR */
        scalar nu_;  /* for NU_SVC, ONE_CLASS, and NU_SVR */
        scalar p_;   /* for EPSILON_SVR */
        label shrinking_;  /* use the shrinking heuristics */
        label probability_; /* do probability estimates */

        //- SVM training data space
        struct svm_node *xSpace;

        //-sum nSV for all time steps
        label sumCtrlPoints2D_;
        label sumCtrlPoints3D_;



    // Private Member Functions

        //- Disallow default bitwise copy construct
        SVMMotionSolver(const SVMMotionSolver&);

        //- Disallow default bitwise assignment
        void operator=(const SVMMotionSolver&);


        //- Make control point IDs.  Constructor helper
        void makeControlIDs();

		//-Make control point IDs for parallel running.  Constructor helper
        void makeParallelControlIDs();

        //-check if two points are the same
        label samePoint(const point& a, const point& b) const;

        //- Set location of points
        void setMovingPoints() const;

        //- SVM regression compute newInternalPionts_
        void SVMRegression();


public:

    //- Runtime type information
    TypeName("SVMMotionSolver");


    // Constructors

        //- Construct from polyMesh
        SVMMotionSolver
        (
            const polyMesh& mesh,
            Istream& msData
        );


    // Destructor

        virtual ~SVMMotionSolver();


    // Member Functions

        //- Return moving point indices in mesh points
        const labelList& movingIDs() const
        {
            return movingIDs_;
        }

        //- Return moving points
        const vectorField& movingPoints() const;

        //- Set motion of moving points
        void setMotion(const vectorField&);


        //- Return point location obtained from the current motion field
        virtual tmp<vectorField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update the mesh corresponding to given map
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
