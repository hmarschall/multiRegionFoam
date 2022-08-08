// Find reference cell
if (p_().needReference())
{
    point refPointi(mesh().solutionDict().subDict("blockSolver").lookup("pRefPoint"));
    label refCelli = mesh().findCell(refPointi);
    whichProcHasRef_[Pstream::myProcNo()] = (refCelli >= 0 ? 1 : 0);
    Pstream::gatherList(whichProcHasRef_);

    if (Pstream::master())
    {
        label sumHasRef = sum(whichProcHasRef_);
        if (sumHasRef == 0)
        {
            FatalError
                << "Unable to set reference cell for field "
                    << p_().name()
                    << nl << "    Reference point pRefPoint"
                    << " found on " << sumHasRef << " domains (should be one)"
                    << nl << exit(FatalError);
        } else if (sumHasRef > 1)
        {
            // On the off chance, the point was accepted by more than
            // one processor, just pick the first one
            label procWithRef = -1;
            forAll(whichProcHasRef_, pi)
            {
                if (whichProcHasRef_[pi])
                {
                    procWithRef = pi;
                }
                break;
            }
            whichProcHasRef_ = labelList(Pstream::nProcs(), 0);
            whichProcHasRef_[procWithRef] = 1;
            Pstream::scatterList(whichProcHasRef_);
        }
    }

    if (whichProcHasRef_[Pstream::myProcNo()])
    {
        pRefCell_ = refCelli;
    } else {
        pRefCell_ = -1;
    }
}
