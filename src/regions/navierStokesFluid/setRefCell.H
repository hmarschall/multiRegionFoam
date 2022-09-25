labelList whichProcHasRef(Pstream::nProcs(), 0);

// Find reference cell
if (p_().needReference())
{
    point refPointi(mesh().solutionDict().subDict("PISO").lookup("pRefPoint"));
    label refCelli = mesh().findCell(refPointi);
    whichProcHasRef[Pstream::myProcNo()] = (refCelli >= 0 ? 1 : 0);
    Pstream::gatherList(whichProcHasRef);

    if (Pstream::master())
    {
        label sumHasRef = sum(whichProcHasRef);
        if (sumHasRef == 0)
        {
            FatalErrorIn("setRefCell")
                << "Unable to set reference cell for field "
                    << p_().name()
                    << nl << "    Reference point pdRefPoint"
                    << " found on " << sumHasRef << " domains (should be one)"
                    << nl << exit(FatalError);
        } else if (sumHasRef > 1)
        {
            // On the off chance, the point was accepted by more than
            // one processor, just pick the first one
            label procWithRef = -1;
            forAll(whichProcHasRef, pi)
            {
                if (whichProcHasRef[pi])
                {
                    procWithRef = pi;
                }
                break;
            }
            whichProcHasRef = labelList(Pstream::nProcs(), 0);
            whichProcHasRef[procWithRef] = 1;
            Pstream::scatterList(whichProcHasRef);
        }
    }

    if (whichProcHasRef[Pstream::myProcNo()])
    {
        pRefCell_ = refCelli;
    } else {
        pRefCell_ = -1;
    }
}

