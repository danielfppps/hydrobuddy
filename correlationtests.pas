{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Copyright (c) 2007, Sergey Bochkanov (ALGLIB project).

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
unit correlationtests;
interface
uses Math, Sysutils, Ap, gammafunc, normaldistr, ibetaf, studenttdistr, correlation;

procedure PearsonCorrelationSignificance(R : Double;
     N : AlglibInteger;
     var BothTails : Double;
     var LeftTail : Double;
     var RightTail : Double);
procedure SpearmanRankCorrelationSignificance(R : Double;
     N : AlglibInteger;
     var BothTails : Double;
     var LeftTail : Double;
     var RightTail : Double);

implementation

function SpearmanTail5(S : Double):Double;forward;
function SpearmanTail6(S : Double):Double;forward;
function SpearmanTail7(S : Double):Double;forward;
function SpearmanTail8(S : Double):Double;forward;
function SpearmanTail9(S : Double):Double;forward;
function SpearmanTail(T : Double; N : AlglibInteger):Double;forward;


(*************************************************************************
Pearson's correlation coefficient significance test

This test checks hypotheses about whether X  and  Y  are  samples  of  two
continuous  distributions  having  zero  correlation  or   whether   their
correlation is non-zero.

The following tests are performed:
    * two-tailed test (null hypothesis - X and Y have zero correlation)
    * left-tailed test (null hypothesis - the correlation  coefficient  is
      greater than or equal to 0)
    * right-tailed test (null hypothesis - the correlation coefficient  is
      less than or equal to 0).

Requirements:
    * the number of elements in each sample is not less than 5
    * normality of distributions of X and Y.

Input parameters:
    R   -   Pearson's correlation coefficient for X and Y
    N   -   number of elements in samples, N>=5.

Output parameters:
    BothTails   -   p-value for two-tailed test.
                    If BothTails is less than the given significance level
                    the null hypothesis is rejected.
    LeftTail    -   p-value for left-tailed test.
                    If LeftTail is less than the given significance level,
                    the null hypothesis is rejected.
    RightTail   -   p-value for right-tailed test.
                    If RightTail is less than the given significance level
                    the null hypothesis is rejected.

  -- ALGLIB --
     Copyright 09.04.2007 by Bochkanov Sergey
*************************************************************************)
procedure PearsonCorrelationSignificance(R : Double;
     N : AlglibInteger;
     var BothTails : Double;
     var LeftTail : Double;
     var RightTail : Double);
var
    T : Double;
    P : Double;
begin
    
    //
    // Some special cases
    //
    if AP_FP_Greater_Eq(R,1) then
    begin
        BothTails := Double(0.0);
        LeftTail := Double(1.0);
        RightTail := Double(0.0);
        Exit;
    end;
    if AP_FP_Less_Eq(R,-1) then
    begin
        BothTails := Double(0.0);
        LeftTail := Double(0.0);
        RightTail := Double(1.0);
        Exit;
    end;
    if N<5 then
    begin
        BothTails := Double(1.0);
        LeftTail := Double(1.0);
        RightTail := Double(1.0);
        Exit;
    end;
    
    //
    // General case
    //
    T := R*Sqrt((N-2)/(1-AP_Sqr(R)));
    P := StudentTDistribution(N-2, T);
    BothTails := 2*Min(P, 1-P);
    LeftTail := P;
    RightTail := 1-P;
end;


(*************************************************************************
Spearman's rank correlation coefficient significance test

This test checks hypotheses about whether X  and  Y  are  samples  of  two
continuous  distributions  having  zero  correlation  or   whether   their
correlation is non-zero.

The following tests are performed:
    * two-tailed test (null hypothesis - X and Y have zero correlation)
    * left-tailed test (null hypothesis - the correlation  coefficient  is
      greater than or equal to 0)
    * right-tailed test (null hypothesis - the correlation coefficient  is
      less than or equal to 0).

Requirements:
    * the number of elements in each sample is not less than 5.

The test is non-parametric and doesn't require distributions X and Y to be
normal.

Input parameters:
    R   -   Spearman's rank correlation coefficient for X and Y
    N   -   number of elements in samples, N>=5.

Output parameters:
    BothTails   -   p-value for two-tailed test.
                    If BothTails is less than the given significance level
                    the null hypothesis is rejected.
    LeftTail    -   p-value for left-tailed test.
                    If LeftTail is less than the given significance level,
                    the null hypothesis is rejected.
    RightTail   -   p-value for right-tailed test.
                    If RightTail is less than the given significance level
                    the null hypothesis is rejected.

  -- ALGLIB --
     Copyright 09.04.2007 by Bochkanov Sergey
*************************************************************************)
procedure SpearmanRankCorrelationSignificance(R : Double;
     N : AlglibInteger;
     var BothTails : Double;
     var LeftTail : Double;
     var RightTail : Double);
var
    T : Double;
    P : Double;
begin
    
    //
    // Special case
    //
    if N<5 then
    begin
        BothTails := Double(1.0);
        LeftTail := Double(1.0);
        RightTail := Double(1.0);
        Exit;
    end;
    
    //
    // General case
    //
    if AP_FP_Greater_Eq(R,1) then
    begin
        T := Double(1.0E10);
    end
    else
    begin
        if AP_FP_Less_Eq(R,-1) then
        begin
            T := -Double(1.0E10);
        end
        else
        begin
            T := R*Sqrt((N-2)/(1-AP_Sqr(R)));
        end;
    end;
    if AP_FP_Less(T,0) then
    begin
        P := SpearmanTail(T, N);
        BothTails := 2*P;
        LeftTail := P;
        RightTail := 1-P;
    end
    else
    begin
        P := SpearmanTail(-T, N);
        BothTails := 2*P;
        LeftTail := 1-P;
        RightTail := P;
    end;
end;


(*************************************************************************
Tail(S, 5)
*************************************************************************)
function SpearmanTail5(S : Double):Double;
begin
    if AP_FP_Less(S,Double(0.000e+00)) then
    begin
        Result := StudentTDistribution(3, -S);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(3.580e+00)) then
    begin
        Result := Double(8.304e-03);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(2.322e+00)) then
    begin
        Result := Double(4.163e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(1.704e+00)) then
    begin
        Result := Double(6.641e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(1.303e+00)) then
    begin
        Result := Double(1.164e-01);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(1.003e+00)) then
    begin
        Result := Double(1.748e-01);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(7.584e-01)) then
    begin
        Result := Double(2.249e-01);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(5.468e-01)) then
    begin
        Result := Double(2.581e-01);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(3.555e-01)) then
    begin
        Result := Double(3.413e-01);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(1.759e-01)) then
    begin
        Result := Double(3.911e-01);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(1.741e-03)) then
    begin
        Result := Double(4.747e-01);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(0.000e+00)) then
    begin
        Result := Double(5.248e-01);
        Exit;
    end;
    Result := 0;
end;


(*************************************************************************
Tail(S, 6)
*************************************************************************)
function SpearmanTail6(S : Double):Double;
begin
    if AP_FP_Less(S,Double(1.001e+00)) then
    begin
        Result := StudentTDistribution(4, -S);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(5.663e+00)) then
    begin
        Result := Double(1.366e-03);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(3.834e+00)) then
    begin
        Result := Double(8.350e-03);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(2.968e+00)) then
    begin
        Result := Double(1.668e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(2.430e+00)) then
    begin
        Result := Double(2.921e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(2.045e+00)) then
    begin
        Result := Double(5.144e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(1.747e+00)) then
    begin
        Result := Double(6.797e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(1.502e+00)) then
    begin
        Result := Double(8.752e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(1.295e+00)) then
    begin
        Result := Double(1.210e-01);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(1.113e+00)) then
    begin
        Result := Double(1.487e-01);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(1.001e+00)) then
    begin
        Result := Double(1.780e-01);
        Exit;
    end;
    Result := 0;
end;


(*************************************************************************
Tail(S, 7)
*************************************************************************)
function SpearmanTail7(S : Double):Double;
begin
    if AP_FP_Less(S,Double(1.001e+00)) then
    begin
        Result := StudentTDistribution(5, -S);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(8.159e+00)) then
    begin
        Result := Double(2.081e-04);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(5.620e+00)) then
    begin
        Result := Double(1.393e-03);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(4.445e+00)) then
    begin
        Result := Double(3.398e-03);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(3.728e+00)) then
    begin
        Result := Double(6.187e-03);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(3.226e+00)) then
    begin
        Result := Double(1.200e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(2.844e+00)) then
    begin
        Result := Double(1.712e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(2.539e+00)) then
    begin
        Result := Double(2.408e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(2.285e+00)) then
    begin
        Result := Double(3.320e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(2.068e+00)) then
    begin
        Result := Double(4.406e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(1.879e+00)) then
    begin
        Result := Double(5.478e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(1.710e+00)) then
    begin
        Result := Double(6.946e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(1.559e+00)) then
    begin
        Result := Double(8.331e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(1.420e+00)) then
    begin
        Result := Double(1.001e-01);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(1.292e+00)) then
    begin
        Result := Double(1.180e-01);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(1.173e+00)) then
    begin
        Result := Double(1.335e-01);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(1.062e+00)) then
    begin
        Result := Double(1.513e-01);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(1.001e+00)) then
    begin
        Result := Double(1.770e-01);
        Exit;
    end;
    Result := 0;
end;


(*************************************************************************
Tail(S, 8)
*************************************************************************)
function SpearmanTail8(S : Double):Double;
begin
    if AP_FP_Less(S,Double(2.001e+00)) then
    begin
        Result := StudentTDistribution(6, -S);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(1.103e+01)) then
    begin
        Result := Double(2.194e-05);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(7.685e+00)) then
    begin
        Result := Double(2.008e-04);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(6.143e+00)) then
    begin
        Result := Double(5.686e-04);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(5.213e+00)) then
    begin
        Result := Double(1.138e-03);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(4.567e+00)) then
    begin
        Result := Double(2.310e-03);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(4.081e+00)) then
    begin
        Result := Double(3.634e-03);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(3.697e+00)) then
    begin
        Result := Double(5.369e-03);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(3.381e+00)) then
    begin
        Result := Double(7.708e-03);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(3.114e+00)) then
    begin
        Result := Double(1.087e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(2.884e+00)) then
    begin
        Result := Double(1.397e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(2.682e+00)) then
    begin
        Result := Double(1.838e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(2.502e+00)) then
    begin
        Result := Double(2.288e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(2.340e+00)) then
    begin
        Result := Double(2.883e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(2.192e+00)) then
    begin
        Result := Double(3.469e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(2.057e+00)) then
    begin
        Result := Double(4.144e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(2.001e+00)) then
    begin
        Result := Double(4.804e-02);
        Exit;
    end;
    Result := 0;
end;


(*************************************************************************
Tail(S, 9)
*************************************************************************)
function SpearmanTail9(S : Double):Double;
begin
    if AP_FP_Less(S,Double(2.001e+00)) then
    begin
        Result := StudentTDistribution(7, -S);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(9.989e+00)) then
    begin
        Result := Double(2.306e-05);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(8.069e+00)) then
    begin
        Result := Double(8.167e-05);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(6.890e+00)) then
    begin
        Result := Double(1.744e-04);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(6.077e+00)) then
    begin
        Result := Double(3.625e-04);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(5.469e+00)) then
    begin
        Result := Double(6.450e-04);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(4.991e+00)) then
    begin
        Result := Double(1.001e-03);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(4.600e+00)) then
    begin
        Result := Double(1.514e-03);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(4.272e+00)) then
    begin
        Result := Double(2.213e-03);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(3.991e+00)) then
    begin
        Result := Double(2.990e-03);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(3.746e+00)) then
    begin
        Result := Double(4.101e-03);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(3.530e+00)) then
    begin
        Result := Double(5.355e-03);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(3.336e+00)) then
    begin
        Result := Double(6.887e-03);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(3.161e+00)) then
    begin
        Result := Double(8.598e-03);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(3.002e+00)) then
    begin
        Result := Double(1.065e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(2.855e+00)) then
    begin
        Result := Double(1.268e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(2.720e+00)) then
    begin
        Result := Double(1.552e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(2.595e+00)) then
    begin
        Result := Double(1.836e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(2.477e+00)) then
    begin
        Result := Double(2.158e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(2.368e+00)) then
    begin
        Result := Double(2.512e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(2.264e+00)) then
    begin
        Result := Double(2.942e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(2.166e+00)) then
    begin
        Result := Double(3.325e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(2.073e+00)) then
    begin
        Result := Double(3.800e-02);
        Exit;
    end;
    if AP_FP_Greater_Eq(S,Double(2.001e+00)) then
    begin
        Result := Double(4.285e-02);
        Exit;
    end;
    Result := 0;
end;


(*************************************************************************
Tail(T,N), accepts T<0
*************************************************************************)
function SpearmanTail(T : Double; N : AlglibInteger):Double;
begin
    if N=5 then
    begin
        Result := SpearmanTail5(-T);
        Exit;
    end;
    if N=6 then
    begin
        Result := SpearmanTail6(-T);
        Exit;
    end;
    if N=7 then
    begin
        Result := SpearmanTail7(-T);
        Exit;
    end;
    if N=8 then
    begin
        Result := SpearmanTail8(-T);
        Exit;
    end;
    if N=9 then
    begin
        Result := SpearmanTail9(-T);
        Exit;
    end;
    Result := StudentTDistribution(N-2, T);
end;


end.