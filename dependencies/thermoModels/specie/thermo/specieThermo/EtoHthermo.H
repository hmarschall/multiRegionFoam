inline scalar Cp
(
    const scalar p,
    const scalar T
) const
{
    #ifdef __clang__
    volatile const scalar cv = Cv(p, T);
    volatile const scalar cpmcv = equationOfState::CpMCv(p, T);
    return cv + cpmcv;
    #else
    return Cv(p, T) + equationOfState::CpMCv(p, T);
    #endif
}

inline scalar Hs
(
    const scalar p,
    const scalar T
) const
{
    #ifdef __clang__
    volatile const scalar es = Es(p, T);
    volatile const scalar rho = equationOfState::rho(p, T);
    return es + p/rho;
    #else
    return Es(p, T) + p/equationOfState::rho(p, T);
    #endif
}

inline scalar Ha
(
    const scalar p,
    const scalar T
) const
{
    #ifdef __clang__
    volatile const scalar ea = Ea(p, T);
    volatile const scalar rho = equationOfState::rho(p, T);
    return ea + p/rho;
    #else
    return Ea(p, T) + p/equationOfState::rho(p, T);
    #endif
}
