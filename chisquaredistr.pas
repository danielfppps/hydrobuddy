{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier

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
unit chisquaredistr;
interface
uses Math, Sysutils, Ap, gammafunc, normaldistr, igammaf;

function ChiSquareDistribution(v : Double; x : Double):Double;
function ChiSquareCDistribution(v : Double; x : Double):Double;
function InvChiSquareDistribution(v : Double; y : Double):Double;

implementation

(*************************************************************************
Chi-square distribution

Returns the area under the left hand tail (from 0 to x)
of the Chi square probability density function with
v degrees of freedom.


                                  x
                                   -
                       1          | |  v/2-1  -t/2
 P( x | v )   =   -----------     |   t      e     dt
                   v/2  -       | |
                  2    | (v/2)   -
                                  0

where x is the Chi-square variable.

The incomplete gamma integral is used, according to the
formula

y = chdtr( v, x ) = igam( v/2.0, x/2.0 ).

The arguments must both be positive.

ACCURACY:

See incomplete gamma function


Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************)
function ChiSquareDistribution(v : Double; x : Double):Double;
begin
    Assert(AP_FP_Greater_Eq(x,0) and AP_FP_Greater_Eq(v,1), 'Domain error in ChiSquareDistribution');
    Result := IncompleteGamma(v/Double(2.0), x/Double(2.0));
end;


(*************************************************************************
Complemented Chi-square distribution

Returns the area under the right hand tail (from x to
infinity) of the Chi square probability density function
with v degrees of freedom:

                                 inf.
                                   -
                       1          | |  v/2-1  -t/2
 P( x | v )   =   -----------     |   t      e     dt
                   v/2  -       | |
                  2    | (v/2)   -
                                  x

where x is the Chi-square variable.

The incomplete gamma integral is used, according to the
formula

y = chdtr( v, x ) = igamc( v/2.0, x/2.0 ).

The arguments must both be positive.

ACCURACY:

See incomplete gamma function

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************)
function ChiSquareCDistribution(v : Double; x : Double):Double;
begin
    Assert(AP_FP_Greater_Eq(x,0) and AP_FP_Greater_Eq(v,1), 'Domain error in ChiSquareDistributionC');
    Result := IncompleteGammaC(v/Double(2.0), x/Double(2.0));
end;


(*************************************************************************
Inverse of complemented Chi-square distribution

Finds the Chi-square argument x such that the integral
from x to infinity of the Chi-square density is equal
to the given cumulative probability y.

This is accomplished using the inverse gamma integral
function and the relation

   x/2 = igami( df/2, y );

ACCURACY:

See inverse incomplete gamma function


Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************)
function InvChiSquareDistribution(v : Double; y : Double):Double;
begin
    Assert(AP_FP_Greater_Eq(y,0) and AP_FP_Less_Eq(y,1) and AP_FP_Greater_Eq(v,1), 'Domain error in InvChiSquareDistribution');
    Result := 2*InvIncompleteGammaC(Double(0.5)*v, y);
end;


end.