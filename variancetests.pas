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
unit variancetests;
interface
uses Math, Sysutils, Ap, gammafunc, normaldistr, ibetaf, fdistr, igammaf, chisquaredistr;

procedure FTest(const X : TReal1DArray;
     N : AlglibInteger;
     const Y : TReal1DArray;
     M : AlglibInteger;
     var BothTails : Double;
     var LeftTail : Double;
     var RightTail : Double);
procedure OneSampleVarianceTest(const X : TReal1DArray;
     N : AlglibInteger;
     Variance : Double;
     var BothTails : Double;
     var LeftTail : Double;
     var RightTail : Double);

implementation

(*************************************************************************
Two-sample F-test

This test checks three hypotheses about dispersions of the given  samples.
The following tests are performed:
    * two-tailed test (null hypothesis - the dispersions are equal)
    * left-tailed test (null hypothesis  -  the  dispersion  of  the first
      sample is greater than or equal to  the  dispersion  of  the  second
      sample).
    * right-tailed test (null hypothesis - the  dispersion  of  the  first
      sample is less than or equal to the dispersion of the second sample)

The test is based on the following assumptions:
    * the given samples have normal distributions
    * the samples are independent.

Input parameters:
    X   -   sample 1. Array whose index goes from 0 to N-1.
    N   -   sample size.
    Y   -   sample 2. Array whose index goes from 0 to M-1.
    M   -   sample size.

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
     Copyright 19.09.2006 by Bochkanov Sergey
*************************************************************************)
procedure FTest(const X : TReal1DArray;
     N : AlglibInteger;
     const Y : TReal1DArray;
     M : AlglibInteger;
     var BothTails : Double;
     var LeftTail : Double;
     var RightTail : Double);
var
    I : AlglibInteger;
    XMean : Double;
    YMean : Double;
    XVar : Double;
    YVar : Double;
    DF1 : AlglibInteger;
    DF2 : AlglibInteger;
    Stat : Double;
begin
    if (N<=2) or (M<=2) then
    begin
        BothTails := Double(1.0);
        LeftTail := Double(1.0);
        RightTail := Double(1.0);
        Exit;
    end;
    
    //
    // Mean
    //
    XMean := 0;
    I:=0;
    while I<=N-1 do
    begin
        XMean := XMean+X[I];
        Inc(I);
    end;
    XMean := XMean/N;
    YMean := 0;
    I:=0;
    while I<=M-1 do
    begin
        YMean := YMean+Y[I];
        Inc(I);
    end;
    YMean := YMean/M;
    
    //
    // Variance (using corrected two-pass algorithm)
    //
    XVar := 0;
    I:=0;
    while I<=N-1 do
    begin
        XVar := XVar+AP_Sqr(X[I]-XMean);
        Inc(I);
    end;
    XVar := XVar/(N-1);
    YVar := 0;
    I:=0;
    while I<=M-1 do
    begin
        YVar := YVar+AP_Sqr(Y[I]-YMean);
        Inc(I);
    end;
    YVar := YVar/(M-1);
    if AP_FP_Eq(XVar,0) or AP_FP_Eq(YVar,0) then
    begin
        BothTails := Double(1.0);
        LeftTail := Double(1.0);
        RightTail := Double(1.0);
        Exit;
    end;
    
    //
    // Statistic
    //
    DF1 := N-1;
    DF2 := M-1;
    Stat := Min(XVar/YVar, YVar/XVar);
    BothTails := 1-(FDistribution(DF1, DF2, 1/Stat)-FDistribution(DF1, DF2, Stat));
    LeftTail := FDistribution(DF1, DF2, XVar/YVar);
    RightTail := 1-LeftTail;
end;


(*************************************************************************
One-sample chi-square test

This test checks three hypotheses about the dispersion of the given sample
The following tests are performed:
    * two-tailed test (null hypothesis - the dispersion equals  the  given
      number)
    * left-tailed test (null hypothesis - the dispersion is  greater  than
      or equal to the given number)
    * right-tailed test (null hypothesis  -  dispersion is  less  than  or
      equal to the given number).

Test is based on the following assumptions:
    * the given sample has a normal distribution.

Input parameters:
    X           -   sample 1. Array whose index goes from 0 to N-1.
    N           -   size of the sample.
    Variance    -   dispersion value to compare with.

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
     Copyright 19.09.2006 by Bochkanov Sergey
*************************************************************************)
procedure OneSampleVarianceTest(const X : TReal1DArray;
     N : AlglibInteger;
     Variance : Double;
     var BothTails : Double;
     var LeftTail : Double;
     var RightTail : Double);
var
    I : AlglibInteger;
    XMean : Double;
    XVar : Double;
    S : Double;
    Stat : Double;
begin
    if N<=1 then
    begin
        BothTails := Double(1.0);
        LeftTail := Double(1.0);
        RightTail := Double(1.0);
        Exit;
    end;
    
    //
    // Mean
    //
    XMean := 0;
    I:=0;
    while I<=N-1 do
    begin
        XMean := XMean+X[I];
        Inc(I);
    end;
    XMean := XMean/N;
    
    //
    // Variance
    //
    XVar := 0;
    I:=0;
    while I<=N-1 do
    begin
        XVar := XVar+AP_Sqr(X[I]-XMean);
        Inc(I);
    end;
    XVar := XVar/(N-1);
    if AP_FP_Eq(XVar,0) then
    begin
        BothTails := Double(1.0);
        LeftTail := Double(1.0);
        RightTail := Double(1.0);
        Exit;
    end;
    
    //
    // Statistic
    //
    Stat := (N-1)*XVar/Variance;
    S := ChiSquareDistribution(N-1, Stat);
    BothTails := 2*Min(S, 1-S);
    LeftTail := S;
    RightTail := 1-LeftTail;
end;


end.