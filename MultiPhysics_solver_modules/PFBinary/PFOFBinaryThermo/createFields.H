  //! Creating phase-field with the option to write
   volScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
      
  //! Creating chemical potential field with the option to write
    volScalarField mu
    (
        IOobject
        (
            "mu",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    //! Creating orientation field with the option to write
    volScalarField theta
    (
        IOobject
        (
            "theta",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
 // The input properties are imported
 #include "readTransportProperties.H"
 #include "readTDB.H"

    //! use "decomposePar -cellDist"
    //! to get cells belonging to domains
    
    /*volScalarField cellDist
    (
        IOobject
        (
            "cellDist",
            "0",
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );*/
    
    //! get a field of global cell indices
    /*volScalarField cellGlobal = 0.0*cellDist;
    
    forAll(cellDist, i)
    {
        cellGlobal[i] = i;
        Info << i << " " << cellGlobal[i] << " " << cellDist[i] << endl;
    }*/
    
    // There may be better approaches for doing the above
    
    /*
    List<List<label> > processCellToGlobalAddr;
    List<label> globalCellToProcessAddr;

    if (Pstream::parRun())
    {
	processCellToGlobalAddr.resize
	(
	    Pstream::nProcs()
	);
        
	//read local cell addressing
	labelIOList localCellProcAddr
	(
	    IOobject
	    (
		"cellProcAddressing",
		mesh.facesInstance(),
		mesh.meshSubDir,
		mesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	    )
	);
	
	processCellToGlobalAddr[Pstream::myProcNo()] = localCellProcAddr;
	
	//send local cell addressing to master process
	if (Pstream::master())
	{
	    for (label jSlave=Pstream::firstSlave(); jSlave<=Pstream::lastSlave(); jSlave++)
	    {
		IPstream fromSlave(Pstream::scheduled, jSlave);
		label nSlaveCells = 0;
		fromSlave >> nSlaveCells;
		processCellToGlobalAddr[jSlave].resize(nSlaveCells);
		labelList& slaveCellProcAddr = processCellToGlobalAddr[jSlave];
		forAll(slaveCellProcAddr, iCell)
		{
		    fromSlave >> slaveCellProcAddr[iCell];
		}
	    }
	}
	else
	{
	    OPstream toMaster (Pstream::scheduled, Pstream::masterNo());
	    toMaster << localCellProcAddr.size();
	    
	    forAll(localCellProcAddr, iCell)
	    {
		toMaster << localCellProcAddr[iCell];
	    }
	}
	
	//redistribute cell addressing to slave processes
	if (Pstream::master())
	{
	    for (label jSlave=Pstream::firstSlave(); jSlave<=Pstream::lastSlave(); jSlave++)
	    {
		OPstream toSlave (Pstream::scheduled, jSlave);
		forAll(processCellToGlobalAddr, iProcess)
		{
		    const labelList& thisProcessAddr = processCellToGlobalAddr[iProcess];
		    const label nCells = thisProcessAddr.size();
		    toSlave << nCells;
		    forAll(thisProcessAddr, jCell)
		    {
			toSlave << thisProcessAddr[jCell];
		    }
		}
	    }
	}
	else
	{
	    IPstream fromMaster(Pstream::scheduled, Pstream::masterNo());
	    forAll(processCellToGlobalAddr, iProcess)
	    {
		labelList& thisProcessAddr = processCellToGlobalAddr[iProcess];
		label nCells = 0;
		fromMaster >> nCells;
		thisProcessAddr.resize(nCells);
		forAll(thisProcessAddr, jCell)
		{
		    fromMaster >> thisProcessAddr[jCell];
		}
	    }
	}

	forAll(processCellToGlobalAddr, jProc)
	{
	    const labelList& jProcessAddr = processCellToGlobalAddr[jProc];
	    forAll(jProcessAddr, iCell)
	    {
		label iGlobalCell = jProcessAddr[iCell];
		globalCellToProcessAddr[iGlobalCell] = iCell;
	    }
	}
    }
    */
    
    
