#include "fvMesh.H"
#include "fvMatrix.H"


namespace Foam
{

template<class T>
fvMatrix<T>& regionType::getCoupledEqn
(
    word name
)
{
//    notImplemented;        
}

template<>
fvMatrix<scalar>& regionType::getCoupledEqn
(
    word name
)
{
    if (fvScalarMatrices.found(name))
    {
        return *fvScalarMatrices[name];
    }
}

template<>
fvMatrix<vector>& regionType::getCoupledEqn
(
    word name
)
{
    if (fvVectorMatrices.found(name))
    {
        return *fvVectorMatrices[name];
    }
}

}
