
int initialize P(( double, double));

:Begin:
:Function:       initialize
:Pattern:        InitializeSeff[r0_Real, a0_Real]
:Arguments:      { r0, a0 }
:ArgumentTypes:  { Real, Real }
:ReturnType:     Integer
:End:

:Evaluate: InitializeSeff::usage = "InitializeSeff[r0, a0] initializes efective source variables."

double PsiPm P(( int, double, double));

:Begin:
:Function:       PsiPm
:Pattern:        PsiPm[m_Integer, r_Real, theta_Real]
:Arguments:      { m, r, theta }
:ArgumentTypes:  { Integer, Real, Real }
:ReturnType:     Real
:End:

:Evaluate: PsiPm::usage = "PsiPm[m, r, theta] gives the m modes of the singular field."


double Seffm P(( int, double, double));

:Begin:
:Function:       Seffm
:Pattern:        Seffm[ m_Integer, r_Real, theta_Real]
:Arguments:      { m, r, theta }
:ArgumentTypes:  { Integer, Real, Real }
:ReturnType:     Real
:End:

:Evaluate: Seffm::usage = "Seffm[m, r, theta] gives the m modes of the effective source."