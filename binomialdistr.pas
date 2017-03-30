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
unit binomialdistr;
interface
uses Math, Sysutils, Ap, gammafunc, normaldistr, ibetaf, nearunityunit;

function BinomialDistribution(k : AlglibInteger;
     n : AlglibInteger;
     p : Double):Double;
function BinomialCDistribution(k : AlglibInteger;
     n : AlglibInteger;
     p : Double):Double;
function InvBinomialDistribution(k : AlglibInteger;
     n : AlglibInteger;
     y : Double):Double;

implementation

(*************************************************************************
Binomial distribution

Returns the sum of the terms 0 through k of the Binomial
probability density:

  k
  --  ( n )   j      n-j
  >   (   )  p  (1-p)
  --  ( j )
 j=0

The terms are not summed directly; instead the incomplete
beta integral is employed, according to the formula

y = bdtr( k, n, p ) = incbet( n-k, k+1, 1-p ).

The arguments must be positive, with p ranging from 0 to 1.

ACCURACY:

Tested at random points (a,b,p), with p between 0 and 1.

              a,b                     Relative error:
arithmetic  domain     # trials      peak         rms
 For p between 0.001 and 1:
   IEEE     0,100       100000      4.3e-15     2.6e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*************************************************************************)
function BinomialDistribution(k : AlglibInteger;
     n : AlglibInteger;
     p : Double):Double;
var
    dk : Double;
    dn : Double;
begin
    Assert(AP_FP_Greater_Eq(p,0) and AP_FP_Less_Eq(p,1), 'Domain error in BinomialDistribution');
    Assert((k>=-1) and (k<=n), 'Domain error in BinomialDistribution');
    if k=-1 then
    begin
        Result := 0;
        Exit;
    end;
    if k=n then
    begin
        Result := 1;
        Exit;
    end;
    dn := n-k;
    if k=0 then
    begin
        dk := power(Double(1.0)-p, dn);
    end
    else
    begin
        dk := k+1;
        dk := IncompleteBeta(dn, dk, Double(1.0)-p);
    end;
    Result := dk;
end;


(*************************************************************************
Complemented binomial distribution

Returns the sum of the terms k+1 through n of the Binomial
probability density:

  n
  --  ( n )   j      n-j
  >   (   )  p  (1-p)
  --  ( j )
 j=k+1

The terms are not summed directly; instead the incomplete
beta integral is employed, according to the formula

y = bdtrc( k, n, p ) = incbet( k+1, n-k, p ).

The arguments must be positive, with p ranging from 0 to 1.

ACCURACY:

Tested at random points (a,b,p).

              a,b                     Relative error:
arithmetic  domain     # trials      peak         rms
 For p between 0.001 and 1:
   IEEE     0,100       100000      6.7e-15     8.2e-16
 For p between 0 and .001:
   IEEE     0,100       100000      1.5e-13     2.7e-15

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*************************************************************************)
function BinomialCDistribution(k : AlglibInteger;
     n : AlglibInteger;
     p : Double):Double;
var
    dk : Double;
    dn : Double;
begin
    Assert(AP_FP_Greater_Eq(p,0) and AP_FP_Less_Eq(p,1), 'Domain error in BinomialDistributionC');
    Assert((k>=-1) and (k<=n), 'Domain error in BinomialDistributionC');
    if k=-1 then
    begin
        Result := 1;
        Exit;
    end;
    if k=n then
    begin
        Result := 0;
        Exit;
    end;
    dn := n-k;
    if k=0 then
    begin
        if AP_FP_Less(p,Double(0.01)) then
        begin
            dk := -expm1(dn*log1p(-p));
        end
        else
        begin
            dk := Double(1.0)-power(Double(1.0)-p, dn);
        end;
    end
    else
    begin
        dk := k+1;
        dk := IncompleteBeta(dk, dn, p);
    end;
    Result := dk;
end;


(*************************************************************************
Inverse binomial distribution

Finds the event probability p such that the sum of the
terms 0 through k of the Binomial probability density
is equal to the given cumulative probability y.

This is accomplished using the inverse beta integral
function and the relation

1 - p = incbi( n-k, k+1, y ).

ACCURACY:

Tested at random points (a,b,p).

              a,b                     Relative error:
arithmetic  domain     # trials      peak         rms
 For p between 0.001 and 1:
   IEEE     0,100       100000      2.3e-14     6.4e-16
   IEEE     0,10000     100000      6.6e-12     1.2e-13
 For p between 10^-6 and 0.001:
   IEEE     0,100       100000      2.0e-12     1.3e-14
   IEEE     0,10000     100000      1.5e-12     3.2e-14

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*************************************************************************)
function InvBinomialDistribution(k : AlglibInteger;
     n : AlglibInteger;
     y : Double):Double;
var
    dk : Double;
    dn : Double;
    p : Double;
begin
    Assert((k>=0) and (k<n), 'Domain error in InvBinomialDistribution');
    dn := n-k;
    if k=0 then
    begin
        if AP_FP_Greater(y,Double(0.8)) then
        begin
            p := -expm1(log1p(y-Double(1.0))/dn);
        end
        else
        begin
            p := Double(1.0)-power(y, Double(1.0)/dn);
        end;
    end
    else
    begin
        dk := k+1;
        p := IncompleteBeta(dn, dk, Double(0.5));
        if AP_FP_Greater(p,Double(0.5)) then
        begin
            p := InvIncompleteBeta(dk, dn, Double(1.0)-y);
        end
        else
        begin
            p := Double(1.0)-InvIncompleteBeta(dn, dk, y);
        end;
    end;
    Result := p;
end;


end.