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
unit stest;
interface
uses Math, Sysutils, Ap, gammafunc, normaldistr, ibetaf, nearunityunit, binomialdistr;

procedure OneSampleSignTest(const X : TReal1DArray;
     N : AlglibInteger;
     Median : Double;
     var BothTails : Double;
     var LeftTail : Double;
     var RightTail : Double);

implementation

(*************************************************************************
Sign test

This test checks three hypotheses about the median of  the  given  sample.
The following tests are performed:
    * two-tailed test (null hypothesis - the median is equal to the  given
      value)
    * left-tailed test (null hypothesis - the median is  greater  than  or
      equal to the given value)
    * right-tailed test (null hypothesis - the  median  is  less  than  or
      equal to the given value)

Requirements:
    * the scale of measurement should be ordinal, interval or ratio  (i.e.
      the test could not be applied to nominal variables).

The test is non-parametric and doesn't require distribution X to be normal

Input parameters:
    X       -   sample. Array whose index goes from 0 to N-1.
    N       -   size of the sample.
    Median  -   assumed median value.

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

While   calculating   p-values   high-precision   binomial    distribution
approximation is used, so significance levels have about 15 exact digits.

  -- ALGLIB --
     Copyright 08.09.2006 by Bochkanov Sergey
*************************************************************************)
procedure OneSampleSignTest(const X : TReal1DArray;
     N : AlglibInteger;
     Median : Double;
     var BothTails : Double;
     var LeftTail : Double;
     var RightTail : Double);
var
    I : AlglibInteger;
    GTCnt : AlglibInteger;
    NECnt : AlglibInteger;
begin
    if N<=1 then
    begin
        BothTails := Double(1.0);
        LeftTail := Double(1.0);
        RightTail := Double(1.0);
        Exit;
    end;
    
    //
    // Calculate:
    // GTCnt - count of x[i]>Median
    // NECnt - count of x[i]<>Median
    //
    GTCnt := 0;
    NECnt := 0;
    I:=0;
    while I<=N-1 do
    begin
        if AP_FP_Greater(X[I],Median) then
        begin
            GTCnt := GTCnt+1;
        end;
        if AP_FP_Neq(X[I],Median) then
        begin
            NECnt := NECnt+1;
        end;
        Inc(I);
    end;
    if NECnt=0 then
    begin
        
        //
        // all x[i] are equal to Median.
        // So we can conclude that Median is a true median :)
        //
        BothTails := Double(0.0);
        LeftTail := Double(0.0);
        RightTail := Double(0.0);
        Exit;
    end;
    BothTails := 2*BinomialDistribution(Min(GTCnt, NECnt-GTCnt), NECnt, Double(0.5));
    LeftTail := BinomialDistribution(GTCnt, NECnt, Double(0.5));
    RightTail := BinomialCDistribution(GTCnt-1, NECnt, Double(0.5));
end;


end.