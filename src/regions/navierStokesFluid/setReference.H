// if (pNeedRef && Pstream::master())
if (pNeedRef && procHasRef)
{
    pEqn.source()[pRefCell_] +=
        pEqn.diag()[pRefCell_]*pRefValue_;

    pEqn.diag()[pRefCell_] +=
        pEqn.diag()[pRefCell_];
}
