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
    myRBF

Description
    Radial basis function motion solver

Author
    Frank Bos, TU Delft.  All rights reserved.

SourceFiles
    myRBF.C

\*---------------------------------------------------------------------------*/

#ifndef myRBF_H
#define myRBF_H

#include "motionSolver.H"
#include "polyMesh.H"
#include "RBFInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class myRBF Declaration
\*---------------------------------------------------------------------------*/

class myRBF
:
    public motionSolver
{
    // Private data

        //- Moving patch names
        wordList movingPatches_;

        //- Static patch names
        wordList staticPatches_;

        //- Coarsening ratio
        label coarseningRatio_;

        //- Include zero motion of static patches in RBF interpolation
        Switch includeStaticPatches_;

        //- Frozen interpolation
        Switch frozenInterpolation_;

        //- Moving point IDs
        labelList movingIDs_;

		labelList movingPointIDs_;//Inserted by Sandeep 25/11/2014

        //- Moving points on the boundary
        mutable vectorField movingPoints_;

        //- Static point IDs
        labelList staticIDs_;

        //- Control point IDs
        labelList controlIDs_;

		labelList markedPoints;

        //- Control points on the boundary
        mutable vectorField controlPoints_;

        //- Internal point IDs
        labelList internalIDs_;

        //- Internal points
        vectorField internalPoints_;

        //- Prescribed motion on control points
        vectorField motion_;

        //- RBF interpolation
        RBFInterpolation interpolation_;

		// My declarations
		label count;
		//vectorField motionOfControl;
		//vectorField globalMotion;
		//vectorField globalMotionSum;
		


    // Private Member Functions

        //- Disallow default bitwise copy construct
        myRBF(const myRBF&);

        //- Disallow default bitwise assignment
        void operator=(const myRBF&);


        //- Make control point IDs.  Constructor helper
        void makeControlIDs();

        //- Set location of points
        void setMovingPoints() const;


public:

    //- Runtime type information
    TypeName("myRBF");


    // Constructors

        //- Construct from polyMesh
        myRBF
        (
            const polyMesh& mesh,
            Istream& msData
        );


    // Destructor

        virtual ~myRBF();


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
        virtual tmp<pointField> curPoints() const;

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
