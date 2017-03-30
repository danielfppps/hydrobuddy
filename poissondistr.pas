{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier

Contributors:
    * Sergey Bochkanov (ALGLIB project). Translation from C to
      pseudocode.

See subroutines comments for additional copyrights.

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the 
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses

>>> END OF LICENSE >>>
*************************************************************************)
unit poissondistr;
interface
uses Math, Sysutils, Ap, gammafunc, normaldistr, igammaf;

function PoissonDistribution(k : AlglibInteger; m : Double):Double;
function PoissonCDistribution(k : AlglibInteger; m : Double):Double;
function InvPoissonDistribution(k : AlglibInteger; y : Double):Double;

implementation

(*************************************************************************
Poisson distribution

Returns the sum of the first k+1 terms of the Poisson
distribution:

  k         j
  --   -m  m
  >   e    --
  --       j!
 j=0

The terms are not summed directly; instead the incomplete
gamma integral is employed, according to the relation

y = pdtr( k, m ) = igamc( k+1, m ).

The arguments must both be positive.
ACCURACY:

See incomplete gamma function

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*************************************************************************)
function PoissonDistribution(k : AlglibInteger; m : Double):Double;
begin
    Assert((k>=0) and AP_FP_Greater(m,0), 'Domain error in PoissonDistribution');
    Result := IncompleteGammaC(k+1, m);
end;


(*************************************************************************
Complemented Poisson distribution

Returns the sum of the terms k+1 to infinity of the Poisson
distribution:

 inf.       j
  --   -m  m
  >   e    --
  --       j!
 j=k+1

The terms are not summed directly; instead the incomplete
gamma integral is employed, according to the formula

y = pdtrc( k, m ) = igam( k+1, m ).

The arguments must both be positive.

ACCURACY:

See incomplete gamma function

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*************************************************************************)
function PoissonCDistribution(k : AlglibInteger; m : Double):Double;
begin
    Assert((k>=0) and AP_FP_Greater(m,0), 'Domain error in PoissonDistributionC');
    Result := IncompleteGamma(k+1, m);
end;


(*************************************************************************
Inverse Poisson distribution

Finds the Poisson variable x such that the integral
from 0 to x of the Poisson density is equal to the
given probability y.

This is accomplished using the inverse gamma integral
function and the relation

   m = igami( k+1, y ).

ACCURACY:

See inverse incomplete gamma function

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*************************************************************************)
function InvPoissonDistribution(k : AlglibInteger; y : Double):Double;
begin
    Assert((k>=0) and AP_FP_Greater_Eq(y,0) and AP_FP_Less(y,1), 'Domain error in InvPoissonDistribution');
    Result := InvIncompleteGammaC(k+1, y);
end;


end.