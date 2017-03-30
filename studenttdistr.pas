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
unit studenttdistr;
interface
uses Math, Sysutils, Ap, gammafunc, normaldistr, ibetaf;

function StudentTDistribution(k : AlglibInteger; t : Double):Double;
function InvStudentTDistribution(k : AlglibInteger; p : Double):Double;

implementation

(*************************************************************************
Student's t distribution

Computes the integral from minus infinity to t of the Student
t distribution with integer k > 0 degrees of freedom:

                                     t
                                     -
                                    | |
             -                      |         2   -(k+1)/2
            | ( (k+1)/2 )           |  (     x   )
      ----------------------        |  ( 1 + --- )        dx
                    -               |  (      k  )
      sqrt( k pi ) | ( k/2 )        |
                                  | |
                                   -
                                  -inf.

Relation to incomplete beta integral:

       1 - stdtr(k,t) = 0.5 * incbet( k/2, 1/2, z )
where
       z = k/(k + t**2).

For t < -2, this is the method of computation.  For higher t,
a direct method is derived from integration by parts.
Since the function is symmetric about t=0, the area under the
right tail of the density is found by calling the function
with -t instead of t.

ACCURACY:

Tested at random 1 <= k <= 25.  The "domain" refers to t.
                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE     -100,-2      50000       5.9e-15     1.4e-15
   IEEE     -2,100      500000       2.7e-15     4.9e-17

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*************************************************************************)
function StudentTDistribution(k : AlglibInteger; t : Double):Double;
var
    x : Double;
    rk : Double;
    z : Double;
    f : Double;
    tz : Double;
    p : Double;
    xsqk : Double;
    j : AlglibInteger;
begin
    Assert(k>0, 'Domain error in StudentTDistribution');
    if AP_FP_Eq(t,0) then
    begin
        Result := Double(0.5);
        Exit;
    end;
    if AP_FP_Less(t,-Double(2.0)) then
    begin
        rk := k;
        z := rk/(rk+t*t);
        Result := Double(0.5)*IncompleteBeta(Double(0.5)*rk, Double(0.5), z);
        Exit;
    end;
    if AP_FP_Less(t,0) then
    begin
        x := -t;
    end
    else
    begin
        x := t;
    end;
    rk := k;
    z := Double(1.0)+x*x/rk;
    if k mod 2<>0 then
    begin
        xsqk := x/sqrt(rk);
        p := arctan(xsqk);
        if k>1 then
        begin
            f := Double(1.0);
            tz := Double(1.0);
            j := 3;
            while (j<=k-2) and AP_FP_Greater(tz/f,MachineEpsilon) do
            begin
                tz := tz*((j-1)/(z*j));
                f := f+tz;
                j := j+2;
            end;
            p := p+f*xsqk/z;
        end;
        p := p*Double(2.0)/PI;
    end
    else
    begin
        f := Double(1.0);
        tz := Double(1.0);
        j := 2;
        while (j<=k-2) and AP_FP_Greater(tz/f,MachineEpsilon) do
        begin
            tz := tz*((j-1)/(z*j));
            f := f+tz;
            j := j+2;
        end;
        p := f*x/sqrt(z*rk);
    end;
    if AP_FP_Less(t,0) then
    begin
        p := -p;
    end;
    Result := Double(0.5)+Double(0.5)*p;
end;


(*************************************************************************
Functional inverse of Student's t distribution

Given probability p, finds the argument t such that stdtr(k,t)
is equal to p.

ACCURACY:

Tested at random 1 <= k <= 100.  The "domain" refers to p:
                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE    .001,.999     25000       5.7e-15     8.0e-16
   IEEE    10^-6,.001    25000       2.0e-12     2.9e-14

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*************************************************************************)
function InvStudentTDistribution(k : AlglibInteger; p : Double):Double;
var
    t : Double;
    rk : Double;
    z : Double;
    rflg : AlglibInteger;
begin
    Assert((k>0) and AP_FP_Greater(p,0) and AP_FP_Less(p,1), 'Domain error in InvStudentTDistribution');
    rk := k;
    if AP_FP_Greater(p,Double(0.25)) and AP_FP_Less(p,Double(0.75)) then
    begin
        if AP_FP_Eq(p,Double(0.5)) then
        begin
            Result := 0;
            Exit;
        end;
        z := Double(1.0)-Double(2.0)*p;
        z := InvIncompleteBeta(Double(0.5), Double(0.5)*rk, absReal(z));
        t := sqrt(rk*z/(Double(1.0)-z));
        if AP_FP_Less(p,Double(0.5)) then
        begin
            t := -t;
        end;
        Result := t;
        Exit;
    end;
    rflg := -1;
    if AP_FP_Greater_Eq(p,Double(0.5)) then
    begin
        p := Double(1.0)-p;
        rflg := 1;
    end;
    z := InvIncompleteBeta(Double(0.5)*rk, Double(0.5), Double(2.0)*p);
    if AP_FP_Less(MaxRealNumber*z,rk) then
    begin
        Result := rflg*MaxRealNumber;
        Exit;
    end;
    t := sqrt(rk/z-rk);
    Result := rflg*t;
end;


end.