inline scalar Cv
(
    const scalar p,
    const scalar T
) const
{
    #ifdef __clang__
    volatile const scalar cp = Cp(p, T);
    volatile const scalar cpmcv = equationOfState::CpMCv(p, T);
    return cp - cpmcv;
    #else
    return Cp(p, T) - equationOfState::CpMCv(p, T);
    #endif
}

inline scalar Es
(
    const scalar p,
    const scalar T
) const
{
    #ifdef __clang__
    volatile const scalar hs = Hs(p, T);
    volatile const scalar rho = equationOfState::rho(p, T);
    return hs - p/rho;
    #else
    return Hs(p, T) - p/equationOfState::rho(p, T);
    #endif
}

inline scalar Ea
(
    const scalar p,
    const scalar T
) const
{
    #ifdef __clang__
    volatile const scalar ha = Ha(p, T);
    volatile const scalar rho = equationOfState::rho(p, T);
    return ha - p/rho;
    #else
    return Ha(p, T) - p/equationOfState::rho(p, T);
    #endif
}
