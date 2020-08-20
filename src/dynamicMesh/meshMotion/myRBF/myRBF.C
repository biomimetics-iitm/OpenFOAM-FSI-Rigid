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

#include "myRBF.H"
#include "addToRunTimeSelectionTable.H"
#include <iostream>


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(myRBF, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        myRBF,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::myRBF::makeControlIDs()
{

//-----------------------------------------------------------------------
		markedPoints.setSize(mesh().nPoints());
		forAll (markedPoints, i)
			{
					markedPoints[i]=0;					
			}
		const pointZone& zone = mesh().pointZones()[0]; //Assuming there is only one point zone for moving mesh
		Info << "Name: "<<zone.name()<<endl;
		Info<<"Size: "<<zone.size()<<endl;

			forAll (zone, i)
			{
				Info<<"Points: "<<zone[i]<<endl;
  			    			
			}	

//--------------------------------------------------------------------------


	// Points that are neither on moving nor on static patches
    // will be marked with 0
    //labelList markedPoints(mesh().nPoints(), 0);

    // Mark all points on moving patches with 1
    label nMovingPoints = 0;

    forAll (movingPatches_, patchI)
    {
        // Find the patch in boundary
        label patchIndex =
            mesh().boundaryMesh().findPatchID(movingPatches_[patchI]);

	
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
            // Grab internal point
            movingIDs_[nMovingPoints] = i;
            nMovingPoints++;
        }
    }

    movingIDs_.setSize(nMovingPoints);

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

        const labelList& mp = mesh().boundaryMesh()[patchIndex].meshPoints();

        forAll (mp, i)
        {
            markedPoints[mp[i]] = -1;
            nStaticPoints++;
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
    // HJ, 24/Mar/2011
    controlIDs_.setSize(movingIDs_.size() + staticIDs_.size());
    motion_.setSize(controlIDs_.size(), vector::zero);

    label nControlPoints = 0;

    forAll (movingPatches_, patchI)
    {
        // Find the patch in boundary
        label patchIndex =
            mesh().boundaryMesh().findPatchID(movingPatches_[patchI]);

        const labelList& mp = mesh().boundaryMesh()[patchIndex].meshPoints();

        for
        (
            label pickedPoint = 0;
            pickedPoint < mp.size();
            pickedPoint += 1
        )
        {
            // Pick point as control point
            controlIDs_[nControlPoints] = mp[pickedPoint]; 			
		
			// Mark the point as picked
            markedPoints[mp[pickedPoint]] = 2;
            nControlPoints++;
        }
    }

    Info<< "Selected " << nControlPoints
        << " control points on moving boundaries" << endl;

    // Resize control IDs
    controlIDs_.setSize(nControlPoints);

   /* // Pick up point locations
    controlPoints_.setSize(nControlPoints);

    // Set control points
    forAll (controlIDs_, i)
    {
        controlPoints_[i] = points[controlIDs_[i]];
		Info<<"controlPoints_[i]:"<<controlPoints_[i]<<endl;
    }*/

//////////////////////////////////////////////////////////////////////////////

/*

	  // Pick up point locations
     controlPoints_.setSize(zone.size(),vector::zero);

    // Set control points 
     forAll(points, i)
      {
		  if (markedPoints[i] == 2) 
	      {controlPoints_[zone.whichPoint(i)] = points[i];}
      }

	vectorField controlPoints_Sum = controlPoints_;
	reduce(controlPoints_Sum,sumOp<vectorField>());



	forAll(controlPoints_, i)
	{	
		if ((controlPoints_Sum[i]).x() < 0)
		{reduce((controlPoints_[i]).x(), minOp<scalar>());}
		else
		{reduce((controlPoints_[i]).x(), maxOp<scalar>());}
		if ((controlPoints_Sum[i]).y() < 0)
		{reduce((controlPoints_[i]).y(), minOp<scalar>());}
		else
		{reduce((controlPoints_[i]).y(), maxOp<scalar>());}
		if ((controlPoints_Sum[i]).z() < 0)
		{reduce((controlPoints_[i]).z(), minOp<scalar>());}
		else
		{reduce((controlPoints_[i]).z(), maxOp<scalar>());}

	}


	*/




	vectorField tempControlPts;

	tempControlPts.setSize(zone.size(),vector::zero);

	forAll(points, i)
      {
		  if (markedPoints[i] == 2) 
	      {tempControlPts[zone.whichPoint(i)] = points[i];}
      }

	vectorField tempControlPtsSum = tempControlPts;
	reduce(tempControlPtsSum,sumOp<vectorField>());



	forAll(tempControlPts, i)
	{	
		if ((tempControlPtsSum[i]).x() < 0)
		{reduce((tempControlPts[i]).x(), minOp<scalar>());}
		else
		{reduce((tempControlPts[i]).x(), maxOp<scalar>());}
		if ((tempControlPtsSum[i]).y() < 0)
		{reduce((tempControlPts[i]).y(), minOp<scalar>());}
		else
		{reduce((tempControlPts[i]).y(), maxOp<scalar>());}
		if ((tempControlPtsSum[i]).z() < 0)
		{reduce((tempControlPts[i]).z(), minOp<scalar>());}
		else
		{reduce((tempControlPts[i]).z(), maxOp<scalar>());}

	}

	
	count=0;

	for (label i =0; i<zone.size(); i= i+coarseningRatio_ )
	{	
		count++;
	}

	controlPoints_.setSize(count,vector::zero);

	label j=0;

	for (label i =0; i<zone.size(); i= i+coarseningRatio_ )
	{	
		controlPoints_[j] =  tempControlPts[i];
		j++;
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

}


void Foam::myRBF::setMovingPoints() const
{
    const pointField& points = mesh().points();

    // Set moving points
    forAll (movingIDs_, i)
    {
        movingPoints_[i] = points[movingIDs_[i]];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myRBF::myRBF
(
    const polyMesh& mesh,
    Istream&
)
:
    motionSolver(mesh),
    movingPatches_(lookup("movingPatches")),
    staticPatches_(lookup("staticPatches")),
	//movingPointIDs_(lookup("movingPoints")),
    coarseningRatio_(readLabel(lookup("coarseningRatio"))),
    includeStaticPatches_(lookup("includeStaticPatches")),
    frozenInterpolation_(lookup("frozenInterpolation")),
    movingIDs_(0),
    movingPoints_(0),
    staticIDs_(0),
    controlIDs_(0),
    controlPoints_(0),
    internalIDs_(0),
    internalPoints_(0),
    motion_(0),
    interpolation_
    (
        subDict("interpolation"),
        controlPoints_,
        internalPoints_
    )
	
{
    makeControlIDs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::myRBF::~myRBF()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::myRBF::setMotion(const vectorField& m)
{
    if (m.size() != movingIDs_.size())
    {
        FatalErrorIn
        (
            "void myRBF::setMotion(const vectorField& m)"
        )   << "Incorrect size of motion points: m = " << m.size()
            << " movingIDs = " << movingIDs_.size()
            << abort(FatalError);
    }

    // Motion of static points is zero and moving points are first
    // in the list.  HJ, 24/Mar/2011
    motion_ = vector::zero;

    forAll (m, i)
    {
        motion_[i] = m[i];
    }

}


const Foam::vectorField& Foam::myRBF::movingPoints() const
{
    // Update moving points based on current mesh
    setMovingPoints();

    return movingPoints_;
}


Foam::tmp<Foam::pointField> Foam::myRBF::curPoints() const
{
    // Prepare new points: same as old point
    tmp<pointField> tcurPoints
    (
        new vectorField(mesh().nPoints(), vector::zero)
    );
    pointField& curPoints = tcurPoints();

    // Add motion to existing points

    // 1. Insert prescribed motion of moving points
    forAll (movingIDs_, i)
    {
        curPoints[movingIDs_[i]] = motion_[i];
    }

    // 2. Insert zero motion of static points
    forAll (staticIDs_, i)
    {
        curPoints[staticIDs_[i]] = vector::zero;
    }

     vectorField localMotion = curPoints;
    const pointZone& zone = mesh().pointZones()[0];
    // put local U into globalU
    //const label pointStart = mesh.Points().start();
    vectorField globalMotion(zone.size(), vector::zero);


    // Set control points and globalMotion
     forAll(localMotion, i)
      {
		  if (markedPoints[i] == 2) 
	      {globalMotion[zone.whichPoint(i)] = localMotion[i];}
      }


	vectorField globalMotionSum = globalMotion;

	reduce(globalMotionSum,sumOp<vectorField>());

	forAll(globalMotion, i)
	{	
		if ((globalMotionSum[i]).x() < 0)
		{reduce((globalMotion[i]).x(), minOp<scalar>());}
		else
		{reduce((globalMotion[i]).x(), maxOp<scalar>());}

		if ((globalMotionSum[i]).y() < 0)
		{reduce((globalMotion[i]).y(), minOp<scalar>());}
		else
		{reduce((globalMotion[i]).y(), maxOp<scalar>());}
		
		if ((globalMotionSum[i]).z() < 0)
		{reduce((globalMotion[i]).z(), minOp<scalar>());}
		else
		{reduce((globalMotion[i]).z(), maxOp<scalar>());}

	}


	//vectorField motionOfControl = globalMotion;

	

	vectorField motionOfControl;

	motionOfControl.setSize(count,vector::zero);

	label j=0;

	for (label i =0; i<zone.size(); i= i+coarseningRatio_ )
	{	
		motionOfControl[j] = globalMotion[i];
		j++;
	}
	
   /* // Set motion of control
    vectorField motionOfControl(controlIDs_.size());

    // 2. Capture positions of control points
    forAll (controlIDs_, i)
    {
        motionOfControl[i] = curPoints[controlIDs_[i]];
		Info<< "motionOfControl[i]: "<<motionOfControl[i]<<endl;
    }

	//Sandeep testing.................................................................................
	Info<< "Total points:" <<mesh().nPoints() <<endl;

	// testing ends....................................................................................*/

    // Call interpolation
    vectorField interpolatedMotion =
        interpolation_.interpolate(motionOfControl);



    // 3. Insert RBF interpolated motion
    forAll (internalIDs_, i)
    {
        curPoints[internalIDs_[i]] = interpolatedMotion[i];
    }
	
	// 4. Add old point positions
    curPoints += mesh().points();

	twoDCorrectPoints(tcurPoints());

	return tcurPoints;
}


void Foam::myRBF::solve()
{}


void Foam::myRBF::updateMesh(const mapPolyMesh&)
{
    // Recalculate control point IDs
    makeControlIDs();
}

