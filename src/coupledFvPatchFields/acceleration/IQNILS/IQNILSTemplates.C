/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  Field<Type> ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::IQNILS<Type>::updateIQNILS(Field<Type> &curFld)
{
    notImplemented
    (
        "IQNILSTemplates.C\n"
        "void IQNILS<Type>::updateIQNILS\n"
        "(\n"
        "    Field<Type> &curFld\n"
        ")\n"
        "not implemented"
    );
}

template<>
inline void Foam::IQNILS<Foam::scalar>::updateIQNILS(Field<scalar> &curFld)
{
    label cols = V_.size();

    RectangularMatrix<scalar> R(cols, cols, 0.0);
    RectangularMatrix<scalar> C(cols, 1);
    RectangularMatrix<scalar> Rcolsum(1, cols);

    DynamicList<Field<scalar> > Q;

    //- QR decomposition with modified Gram-Schmidt procedure
    for (label i = 0; i < cols; i++)
    {
        Q.append(V_[cols-1-i]);
    }

    for (label i = 0; i < cols; i++)
    {
        R[i][i] = Foam::sqrt(sum(Q[i] * Q[i]));
        Q[i] /= R[i][i];

        // Orthogonalize columns to the right of column i
        for (label j = i+1; j < cols; j++)
        {
            R[i][j] = sum(Q[i] * Q[j]);
            Q[j] -= R[i][j]*Q[i];
        }

        // Project minus the residual vector on the Q
        C[i][0] = sum
            (
                Q[i]
              * (
                    - this->resFld_
                )
            );
    }

    // Solve the upper triangular system
    for (label j = 0; j < cols; j++)
    {
        Rcolsum[0][j] = 0.0;

        for (label i = 0; i < j+1; i++)
        {
            Rcolsum[0][j] += cmptMag(R[i][j]);
        }
    }

    scalar epsilon = 1.0E-10*max(Rcolsum);

    for (label i = 0; i < cols; i++)
    {
        if (cmptMag(R[i][i]) > epsilon)
        {
            for (label j = i + 1; j < cols; j++)
            {
                R[i][j] /= R[i][i];
            }

            C[i][0] /= R[i][i];
            R[i][i] = 1.0;
        }
    }

    for (label j = cols-1; j >= 0; j--)
    {
        if (cmptMag(R[j][j]) > epsilon)
        {
            for (label i = 0; i < j; i++)
            {
                C[i][0] -= C[j][0]*R[i][j];
            }
        }
        else
        {
            C[j][0] = 0.0;
        }
    }

    //- Update curent field
    for (label i = 0; i < cols; i++)
    {
        curFld += W_[i]*C[cols-1-i][0];
    }
}


template<>
inline void Foam::IQNILS<Foam::vector>::updateIQNILS(Field<vector> &curFld)
{
    label cols = V_.size();

    RectangularMatrix<scalar> R(cols, cols, 0.0);
    RectangularMatrix<scalar> C(cols, 1);
    RectangularMatrix<scalar> Rcolsum(1, cols);

    DynamicList<Field<vector> > Q;

    //- QR decomposition with modified Gram-Schmidt procedure
    for (label i = 0; i < cols; i++)
    {
        Q.append(V_[cols-1-i]);
    }

    for (label i = 0; i < cols; i++)
    {
        R[i][i] = Foam::sqrt(sum(Q[i] & Q[i]));
        Q[i] /= R[i][i];

        // Orthogonalize columns to the right of column i
        for (label j = i+1; j < cols; j++)
        {
            R[i][j] = sum(Q[i] & Q[j]);
            Q[j] -= R[i][j]*Q[i];
        }

        // Project minus the residual vector on the Q
        C[i][0] = sum
            (
                Q[i]
              & (
                    - this->resFld_
                )
            );
    }

    // Solve the upper triangular system
    for (label j = 0; j < cols; j++)
    {
        Rcolsum[0][j] = 0.0;

        for (label i = 0; i < j+1; i++)
        {
            Rcolsum[0][j] += cmptMag(R[i][j]);
        }
    }

    scalar epsilon = 1.0E-10*max(Rcolsum);

    for (label i = 0; i < cols; i++)
    {
        if (cmptMag(R[i][i]) > epsilon)
        {
            for (label j = i + 1; j < cols; j++)
            {
                R[i][j] /= R[i][i];
            }

            C[i][0] /= R[i][i];
            R[i][i] = 1.0;
        }
    }

    for (label j = cols-1; j >= 0; j--)
    {
        if (cmptMag(R[j][j]) > epsilon)
        {
            for (label i = 0; i < j; i++)
            {
                C[i][0] -= C[j][0]*R[i][j];
            }
        }
        else
        {
            C[j][0] = 0.0;
        }
    }

    //- Update curent field
    for (label i = 0; i < cols; i++)
    {
        curFld += W_[i]*C[cols-1-i][0];
    }
}
// ************************************************************************* //
