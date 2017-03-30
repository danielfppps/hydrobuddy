{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Copyright (c) 2005-2009, Sergey Bochkanov (ALGLIB project).

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
unit autogk;
interface
uses Math, Sysutils, Ap, tsort, hblas, reflections, creflections, sblas, ablasf, ablas, ortfac, blas, rotations, hsschur, evd, gammafunc, gq, gkq;

type
(*************************************************************************
Integration report:
* TerminationType = completetion code:
    * -5    non-convergence of Gauss-Kronrod nodes
            calculation subroutine.
    * -1    incorrect parameters were specified
    *  1    OK
* Rep.NFEV countains number of function calculations
* Rep.NIntervals contains number of intervals [a,b]
  was partitioned into.
*************************************************************************)
AutoGKReport = record
    TerminationType : AlglibInteger;
    NFEV : AlglibInteger;
    NIntervals : AlglibInteger;
end;


AutoGKInternalState = record
    A : Double;
    B : Double;
    Eps : Double;
    XWidth : Double;
    X : Double;
    F : Double;
    Info : AlglibInteger;
    R : Double;
    Heap : TReal2DArray;
    HeapSize : AlglibInteger;
    HeapWidth : AlglibInteger;
    HeapUsed : AlglibInteger;
    SumErr : Double;
    SumAbs : Double;
    QN : TReal1DArray;
    WG : TReal1DArray;
    WK : TReal1DArray;
    WR : TReal1DArray;
    N : AlglibInteger;
    RState : RCommState;
end;


(*************************************************************************
This structure stores internal state of the integration algorithm  between
subsequent calls of the AutoGKIteration() subroutine.
*************************************************************************)
AutoGKState = record
    A : Double;
    B : Double;
    Alpha : Double;
    Beta : Double;
    XWidth : Double;
    X : Double;
    XMinusA : Double;
    BMinusX : Double;
    F : Double;
    WrapperMode : AlglibInteger;
    InternalState : AutoGKInternalState;
    RState : RCommState;
    V : Double;
    TerminationType : AlglibInteger;
    NFEV : AlglibInteger;
    NIntervals : AlglibInteger;
end;



procedure AutoGKSmooth(A : Double; B : Double; var State : AutoGKState);
procedure AutoGKSmoothW(A : Double;
     B : Double;
     XWidth : Double;
     var State : AutoGKState);
procedure AutoGKSingular(A : Double;
     B : Double;
     Alpha : Double;
     Beta : Double;
     var State : AutoGKState);
function AutoGKIteration(var State : AutoGKState):Boolean;
procedure AutoGKResults(const State : AutoGKState;
     var V : Double;
     var Rep : AutoGKReport);

implementation

procedure AutoGKInternalPrepare(A : Double;
     B : Double;
     Eps : Double;
     XWidth : Double;
     var State : AutoGKInternalState);forward;
function AutoGKInternalIteration(var State : AutoGKInternalState):Boolean;forward;
procedure MHeapPop(var Heap : TReal2DArray;
     HeapSize : AlglibInteger;
     HeapWidth : AlglibInteger);forward;
procedure MHeapPush(var Heap : TReal2DArray;
     HeapSize : AlglibInteger;
     HeapWidth : AlglibInteger);forward;
procedure MHeapResize(var Heap : TReal2DArray;
     var HeapSize : AlglibInteger;
     NewHeapSize : AlglibInteger;
     HeapWidth : AlglibInteger);forward;


(*************************************************************************
Integration of a smooth function F(x) on a finite interval [a,b].

Fast-convergent algorithm based on a Gauss-Kronrod formula is used. Result
is calculated with accuracy close to the machine precision.

Algorithm works well only with smooth integrands.  It  may  be  used  with
continuous non-smooth integrands, but with  less  performance.

It should never be used with integrands which have integrable singularities
at lower or upper limits - algorithm may crash. Use AutoGKSingular in such
cases.

INPUT PARAMETERS:
    A, B    -   interval boundaries (A<B, A=B or A>B)
    
OUTPUT PARAMETERS
    State   -   structure which stores algorithm state between  subsequent
                calls of AutoGKIteration.  Used for reverse communication.
                This structure should be  passed  to  the  AutoGKIteration
                subroutine.

SEE ALSO
    AutoGKSmoothW, AutoGKSingular, AutoGKIteration, AutoGKResults.
    

  -- ALGLIB --
     Copyright 06.05.2009 by Bochkanov Sergey
*************************************************************************)
procedure AutoGKSmooth(A : Double; B : Double; var State : AutoGKState);
begin
    AutoGKSmoothW(A, B, Double(0.0), State);
end;


(*************************************************************************
Integration of a smooth function F(x) on a finite interval [a,b].

This subroutine is same as AutoGKSmooth(), but it guarantees that interval
[a,b] is partitioned into subintervals which have width at most XWidth.

Subroutine  can  be  used  when  integrating nearly-constant function with
narrow "bumps" (about XWidth wide). If "bumps" are too narrow, AutoGKSmooth
subroutine can overlook them.

INPUT PARAMETERS:
    A, B    -   interval boundaries (A<B, A=B or A>B)

OUTPUT PARAMETERS
    State   -   structure which stores algorithm state between  subsequent
                calls of AutoGKIteration.  Used for reverse communication.
                This structure should be  passed  to  the  AutoGKIteration
                subroutine.

SEE ALSO
    AutoGKSmooth, AutoGKSingular, AutoGKIteration, AutoGKResults.


  -- ALGLIB --
     Copyright 06.05.2009 by Bochkanov Sergey
*************************************************************************)
procedure AutoGKSmoothW(A : Double;
     B : Double;
     XWidth : Double;
     var State : AutoGKState);
begin
    State.WrapperMode := 0;
    State.A := A;
    State.B := B;
    State.XWidth := XWidth;
    SetLength(State.RState.RA, 10+1);
    State.RState.Stage := -1;
end;


(*************************************************************************
Integration on a finite interval [A,B].
Integrand have integrable singularities at A/B.

F(X) must diverge as "(x-A)^alpha" at A, as "(B-x)^beta" at B,  with known
alpha/beta (alpha>-1, beta>-1).  If alpha/beta  are  not known,  estimates
from below can be used (but these estimates should be greater than -1 too).

One  of  alpha/beta variables (or even both alpha/beta) may be equal to 0,
which means than function F(x) is non-singular at A/B. Anyway (singular at
bounds or not), function F(x) is supposed to be continuous on (A,B).

Fast-convergent algorithm based on a Gauss-Kronrod formula is used. Result
is calculated with accuracy close to the machine precision.

INPUT PARAMETERS:
    A, B    -   interval boundaries (A<B, A=B or A>B)
    Alpha   -   power-law coefficient of the F(x) at A,
                Alpha>-1
    Beta    -   power-law coefficient of the F(x) at B,
                Beta>-1

OUTPUT PARAMETERS
    State   -   structure which stores algorithm state between  subsequent
                calls of AutoGKIteration.  Used for reverse communication.
                This structure should be  passed  to  the  AutoGKIteration
                subroutine.

SEE ALSO
    AutoGKSmooth, AutoGKSmoothW, AutoGKIteration, AutoGKResults.


  -- ALGLIB --
     Copyright 06.05.2009 by Bochkanov Sergey
*************************************************************************)
procedure AutoGKSingular(A : Double;
     B : Double;
     Alpha : Double;
     Beta : Double;
     var State : AutoGKState);
begin
    State.WrapperMode := 1;
    State.A := A;
    State.B := B;
    State.Alpha := Alpha;
    State.Beta := Beta;
    State.XWidth := Double(0.0);
    SetLength(State.RState.RA, 10+1);
    State.RState.Stage := -1;
end;


(*************************************************************************
One step of adaptive integration process.

Called after initialization with one of AutoGKXXX subroutines.
See HTML documentation for examples.

Input parameters:
    State   -   structure which stores algorithm state between  calls  and
                which  is  used  for   reverse   communication.   Must  be
                initialized with one of AutoGKXXX subroutines.

If suborutine returned False, iterative proces has converged. If subroutine
returned True, caller should calculate function value State.F  at  State.X
and call AutoGKIteration again.

NOTE:

When integrating "difficult" functions with integrable singularities like

    F(x) = (x-A)^alpha * (B-x)^beta

subroutine may require the value of F at points which are too close to A/B.
Sometimes to calculate integral with high enough precision we  may need to
calculate F(A+delta) when delta is less than machine  epsilon.  In  finite
precision arithmetics A+delta will be effectively equal to A,  so  we  may
find us in situation when  we  are  trying  to  calculate  something  like
1/sqrt(1-1).

To avoid  such  situations,  AutoGKIteration  subroutine  fills  not  only
State.X  field,  but  also   State.XMinusA   (which  equals  to  X-A)  and
State.BMinusX  (which  equals to B-X) fields.  If X is too close to A or B
(X-A<0.001*A, or B-X<0.001*B, for example) use  these  fields  instead  of
State.X


  -- ALGLIB --
     Copyright 07.05.2009 by Bochkanov Sergey
*************************************************************************)
function AutoGKIteration(var State : AutoGKState):Boolean;
var
    S : Double;
    Tmp : Double;
    Eps : Double;
    A : Double;
    B : Double;
    X : Double;
    T : Double;
    Alpha : Double;
    Beta : Double;
    V1 : Double;
    V2 : Double;
label
lbl_5, lbl_0, lbl_6, lbl_3, lbl_9, lbl_1, lbl_10, lbl_11, lbl_2, lbl_12, lbl_7, lbl_rcomm;
begin
    
    //
    // Reverse communication preparations
    // I know it looks ugly, but it works the same way
    // anywhere from C++ to Python.
    //
    // This code initializes locals by:
    // * random values determined during code
    //   generation - on first subroutine call
    // * values from previous call - on subsequent calls
    //
    if State.RState.Stage>=0 then
    begin
        S := State.RState.RA[0];
        Tmp := State.RState.RA[1];
        Eps := State.RState.RA[2];
        A := State.RState.RA[3];
        B := State.RState.RA[4];
        X := State.RState.RA[5];
        T := State.RState.RA[6];
        Alpha := State.RState.RA[7];
        Beta := State.RState.RA[8];
        V1 := State.RState.RA[9];
        V2 := State.RState.RA[10];
    end
    else
    begin
        S := -983;
        Tmp := -989;
        Eps := -834;
        A := 900;
        B := -287;
        X := 364;
        T := 214;
        Alpha := -338;
        Beta := -686;
        V1 := 912;
        V2 := 585;
    end;
    if State.RState.Stage=0 then
    begin
        goto lbl_0;
    end;
    if State.RState.Stage=1 then
    begin
        goto lbl_1;
    end;
    if State.RState.Stage=2 then
    begin
        goto lbl_2;
    end;
    
    //
    // Routine body
    //
    Eps := 0;
    A := State.A;
    B := State.B;
    Alpha := State.Alpha;
    Beta := State.Beta;
    State.TerminationType := -1;
    State.NFEV := 0;
    State.NIntervals := 0;
    
    //
    // smooth function  at a finite interval
    //
    if State.WrapperMode<>0 then
    begin
        goto lbl_3;
    end;
    
    //
    // special case
    //
    if AP_FP_Eq(A,B) then
    begin
        State.TerminationType := 1;
        State.V := 0;
        Result := False;
        Exit;
    end;
    
    //
    // general case
    //
    AutoGKInternalPrepare(A, B, Eps, State.XWidth, State.InternalState);
lbl_5:
    if not AutoGKInternalIteration(State.InternalState) then
    begin
        goto lbl_6;
    end;
    X := State.InternalState.X;
    State.X := X;
    State.XMinusA := X-A;
    State.BMinusX := B-X;
    State.RState.Stage := 0;
    goto lbl_rcomm;
lbl_0:
    State.NFEV := State.NFEV+1;
    State.InternalState.F := State.F;
    goto lbl_5;
lbl_6:
    State.V := State.InternalState.R;
    State.TerminationType := State.InternalState.Info;
    State.NIntervals := State.InternalState.HeapUsed;
    Result := False;
    Exit;
lbl_3:
    
    //
    // function with power-law singularities at the ends of a finite interval
    //
    if State.WrapperMode<>1 then
    begin
        goto lbl_7;
    end;
    
    //
    // test coefficients
    //
    if AP_FP_Less_Eq(Alpha,-1) or AP_FP_Less_Eq(Beta,-1) then
    begin
        State.TerminationType := -1;
        State.V := 0;
        Result := False;
        Exit;
    end;
    
    //
    // special cases
    //
    if AP_FP_Eq(A,B) then
    begin
        State.TerminationType := 1;
        State.V := 0;
        Result := False;
        Exit;
    end;
    
    //
    // reduction to general form
    //
    if AP_FP_Less(A,B) then
    begin
        S := +1;
    end
    else
    begin
        S := -1;
        Tmp := A;
        A := B;
        B := Tmp;
        Tmp := Alpha;
        Alpha := Beta;
        Beta := Tmp;
    end;
    Alpha := Min(Alpha, 0);
    Beta := Min(Beta, 0);
    
    //
    // first, integrate left half of [a,b]:
    //     integral(f(x)dx, a, (b+a)/2) =
    //     = 1/(1+alpha) * integral(t^(-alpha/(1+alpha))*f(a+t^(1/(1+alpha)))dt, 0, (0.5*(b-a))^(1+alpha))
    //
    AutoGKInternalPrepare(0, Power(Double(0.5)*(B-A), 1+Alpha), Eps, State.XWidth, State.InternalState);
lbl_9:
    if not AutoGKInternalIteration(State.InternalState) then
    begin
        goto lbl_10;
    end;
    
    //
    // Fill State.X, State.XMinusA, State.BMinusX.
    // Latter two are filled correctly even if B<A.
    //
    X := State.InternalState.X;
    T := Power(X, 1/(1+Alpha));
    State.X := A+T;
    if AP_FP_Greater(S,0) then
    begin
        State.XMinusA := T;
        State.BMinusX := B-(A+T);
    end
    else
    begin
        State.XMinusA := A+T-B;
        State.BMinusX := -T;
    end;
    State.RState.Stage := 1;
    goto lbl_rcomm;
lbl_1:
    if AP_FP_Neq(Alpha,0) then
    begin
        State.InternalState.F := State.F*Power(X, -Alpha/(1+Alpha))/(1+Alpha);
    end
    else
    begin
        State.InternalState.F := State.F;
    end;
    State.NFEV := State.NFEV+1;
    goto lbl_9;
lbl_10:
    V1 := State.InternalState.R;
    State.NIntervals := State.NIntervals+State.InternalState.HeapUsed;
    
    //
    // then, integrate right half of [a,b]:
    //     integral(f(x)dx, (b+a)/2, b) =
    //     = 1/(1+beta) * integral(t^(-beta/(1+beta))*f(b-t^(1/(1+beta)))dt, 0, (0.5*(b-a))^(1+beta))
    //
    AutoGKInternalPrepare(0, Power(Double(0.5)*(B-A), 1+Beta), Eps, State.XWidth, State.InternalState);
lbl_11:
    if not AutoGKInternalIteration(State.InternalState) then
    begin
        goto lbl_12;
    end;
    
    //
    // Fill State.X, State.XMinusA, State.BMinusX.
    // Latter two are filled correctly (X-A, B-X) even if B<A.
    //
    X := State.InternalState.X;
    T := Power(X, 1/(1+Beta));
    State.X := B-T;
    if AP_FP_Greater(S,0) then
    begin
        State.XMinusA := B-T-A;
        State.BMinusX := T;
    end
    else
    begin
        State.XMinusA := -T;
        State.BMinusX := A-(B-T);
    end;
    State.RState.Stage := 2;
    goto lbl_rcomm;
lbl_2:
    if AP_FP_Neq(Beta,0) then
    begin
        State.InternalState.F := State.F*Power(X, -Beta/(1+Beta))/(1+Beta);
    end
    else
    begin
        State.InternalState.F := State.F;
    end;
    State.NFEV := State.NFEV+1;
    goto lbl_11;
lbl_12:
    V2 := State.InternalState.R;
    State.NIntervals := State.NIntervals+State.InternalState.HeapUsed;
    
    //
    // final result
    //
    State.V := S*(V1+V2);
    State.TerminationType := 1;
    Result := False;
    Exit;
lbl_7:
    Result := False;
    Exit;
    
    //
    // Saving state
    //
lbl_rcomm:
    Result := True;
    State.RState.RA[0] := S;
    State.RState.RA[1] := Tmp;
    State.RState.RA[2] := Eps;
    State.RState.RA[3] := A;
    State.RState.RA[4] := B;
    State.RState.RA[5] := X;
    State.RState.RA[6] := T;
    State.RState.RA[7] := Alpha;
    State.RState.RA[8] := Beta;
    State.RState.RA[9] := V1;
    State.RState.RA[10] := V2;
end;


(*************************************************************************
Adaptive integration results

Called after AutoGKIteration returned False.

Input parameters:
    State   -   algorithm state (used by AutoGKIteration).

Output parameters:
    V       -   integral(f(x)dx,a,b)
    Rep     -   optimization report (see AutoGKReport description)

  -- ALGLIB --
     Copyright 14.11.2007 by Bochkanov Sergey
*************************************************************************)
procedure AutoGKResults(const State : AutoGKState;
     var V : Double;
     var Rep : AutoGKReport);
begin
    V := State.V;
    Rep.TerminationType := State.TerminationType;
    Rep.NFEV := State.NFEV;
    Rep.NIntervals := State.NIntervals;
end;


(*************************************************************************
Internal AutoGK subroutine
eps<0   - error
eps=0   - automatic eps selection

width<0 -   error
width=0 -   no width requirements
*************************************************************************)
procedure AutoGKInternalPrepare(A : Double;
     B : Double;
     Eps : Double;
     XWidth : Double;
     var State : AutoGKInternalState);
begin
    
    //
    // Save settings
    //
    State.A := A;
    State.B := B;
    State.Eps := Eps;
    State.XWidth := XWidth;
    
    //
    // Prepare RComm structure
    //
    SetLength(State.RState.IA, 3+1);
    SetLength(State.RState.RA, 8+1);
    State.RState.Stage := -1;
end;


(*************************************************************************
Internal AutoGK subroutine
*************************************************************************)
function AutoGKInternalIteration(var State : AutoGKInternalState):Boolean;
var
    C1 : Double;
    C2 : Double;
    I : AlglibInteger;
    J : AlglibInteger;
    IntG : Double;
    IntK : Double;
    IntA : Double;
    V : Double;
    TA : Double;
    TB : Double;
    NS : AlglibInteger;
    QEps : Double;
    Info : AlglibInteger;
label
lbl_5, lbl_0, lbl_7, lbl_3, lbl_8, lbl_11, lbl_1, lbl_13, lbl_10, lbl_4, lbl_14, lbl_16, lbl_19, lbl_2, lbl_21, lbl_18, lbl_15, lbl_rcomm;
begin
    
    //
    // Reverse communication preparations
    // I know it looks ugly, but it works the same way
    // anywhere from C++ to Python.
    //
    // This code initializes locals by:
    // * random values determined during code
    //   generation - on first subroutine call
    // * values from previous call - on subsequent calls
    //
    if State.RState.Stage>=0 then
    begin
        I := State.RState.IA[0];
        J := State.RState.IA[1];
        NS := State.RState.IA[2];
        Info := State.RState.IA[3];
        C1 := State.RState.RA[0];
        C2 := State.RState.RA[1];
        IntG := State.RState.RA[2];
        IntK := State.RState.RA[3];
        IntA := State.RState.RA[4];
        V := State.RState.RA[5];
        TA := State.RState.RA[6];
        TB := State.RState.RA[7];
        QEps := State.RState.RA[8];
    end
    else
    begin
        I := 497;
        J := -271;
        NS := -581;
        Info := 745;
        C1 := -533;
        C2 := -77;
        IntG := 678;
        IntK := -293;
        IntA := 316;
        V := 647;
        TA := -756;
        TB := 830;
        QEps := -871;
    end;
    if State.RState.Stage=0 then
    begin
        goto lbl_0;
    end;
    if State.RState.Stage=1 then
    begin
        goto lbl_1;
    end;
    if State.RState.Stage=2 then
    begin
        goto lbl_2;
    end;
    
    //
    // Routine body
    //
    
    //
    // initialize quadratures.
    // use 15-point Gauss-Kronrod formula.
    //
    State.N := 15;
    GKQGenerateGaussLegendre(State.N, Info, State.QN, State.WK, State.WG);
    if Info<0 then
    begin
        State.Info := -5;
        State.R := 0;
        Result := False;
        Exit;
    end;
    SetLength(State.WR, State.N);
    I:=0;
    while I<=State.N-1 do
    begin
        if I=0 then
        begin
            State.WR[I] := Double(0.5)*AbsReal(State.QN[1]-State.QN[0]);
            Inc(I);
            Continue;
        end;
        if I=State.N-1 then
        begin
            State.WR[State.N-1] := Double(0.5)*AbsReal(State.QN[State.N-1]-State.QN[State.N-2]);
            Inc(I);
            Continue;
        end;
        State.WR[I] := Double(0.5)*AbsReal(State.QN[I-1]-State.QN[I+1]);
        Inc(I);
    end;
    
    //
    // special case
    //
    if AP_FP_Eq(State.A,State.B) then
    begin
        State.Info := 1;
        State.R := 0;
        Result := False;
        Exit;
    end;
    
    //
    // test parameters
    //
    if AP_FP_Less(State.Eps,0) or AP_FP_Less(State.XWidth,0) then
    begin
        State.Info := -1;
        State.R := 0;
        Result := False;
        Exit;
    end;
    State.Info := 1;
    if AP_FP_Eq(State.Eps,0) then
    begin
        State.Eps := 1000*MachineEpsilon;
    end;
    
    //
    // First, prepare heap
    // * column 0   -   absolute error
    // * column 1   -   integral of a F(x) (calculated using Kronrod extension nodes)
    // * column 2   -   integral of a |F(x)| (calculated using modified rect. method)
    // * column 3   -   left boundary of a subinterval
    // * column 4   -   right boundary of a subinterval
    //
    if AP_FP_NEq(State.XWidth,0) then
    begin
        goto lbl_3;
    end;
    
    //
    // no maximum width requirements
    // start from one big subinterval
    //
    State.HeapWidth := 5;
    State.HeapSize := 1;
    State.HeapUsed := 1;
    SetLength(State.Heap, State.HeapSize, State.HeapWidth);
    C1 := Double(0.5)*(State.B-State.A);
    C2 := Double(0.5)*(State.B+State.A);
    IntG := 0;
    IntK := 0;
    IntA := 0;
    I := 0;
lbl_5:
    if I>State.N-1 then
    begin
        goto lbl_7;
    end;
    
    //
    // obtain F
    //
    State.X := C1*State.QN[I]+C2;
    State.RState.Stage := 0;
    goto lbl_rcomm;
lbl_0:
    V := State.F;
    
    //
    // Gauss-Kronrod formula
    //
    IntK := IntK+V*State.WK[I];
    if I mod 2=1 then
    begin
        IntG := IntG+V*State.WG[I];
    end;
    
    //
    // Integral |F(x)|
    // Use rectangles method
    //
    IntA := IntA+AbsReal(V)*State.WR[I];
    I := I+1;
    goto lbl_5;
lbl_7:
    IntK := IntK*(State.B-State.A)*Double(0.5);
    IntG := IntG*(State.B-State.A)*Double(0.5);
    IntA := IntA*(State.B-State.A)*Double(0.5);
    State.Heap[0,0] := AbsReal(IntG-IntK);
    State.Heap[0,1] := IntK;
    State.Heap[0,2] := IntA;
    State.Heap[0,3] := State.A;
    State.Heap[0,4] := State.B;
    State.SumErr := State.Heap[0,0];
    State.SumAbs := AbsReal(IntA);
    goto lbl_4;
lbl_3:
    
    //
    // maximum subinterval should be no more than XWidth.
    // so we create Ceil((B-A)/XWidth)+1 small subintervals
    //
    NS := Ceil(AbsReal(State.B-State.A)/State.XWidth)+1;
    State.HeapSize := NS;
    State.HeapUsed := NS;
    State.HeapWidth := 5;
    SetLength(State.Heap, State.HeapSize, State.HeapWidth);
    State.SumErr := 0;
    State.SumAbs := 0;
    J := 0;
lbl_8:
    if J>NS-1 then
    begin
        goto lbl_10;
    end;
    TA := State.A+J*(State.B-State.A)/NS;
    TB := State.A+(J+1)*(State.B-State.A)/NS;
    C1 := Double(0.5)*(TB-TA);
    C2 := Double(0.5)*(TB+TA);
    IntG := 0;
    IntK := 0;
    IntA := 0;
    I := 0;
lbl_11:
    if I>State.N-1 then
    begin
        goto lbl_13;
    end;
    
    //
    // obtain F
    //
    State.X := C1*State.QN[I]+C2;
    State.RState.Stage := 1;
    goto lbl_rcomm;
lbl_1:
    V := State.F;
    
    //
    // Gauss-Kronrod formula
    //
    IntK := IntK+V*State.WK[I];
    if I mod 2=1 then
    begin
        IntG := IntG+V*State.WG[I];
    end;
    
    //
    // Integral |F(x)|
    // Use rectangles method
    //
    IntA := IntA+AbsReal(V)*State.WR[I];
    I := I+1;
    goto lbl_11;
lbl_13:
    IntK := IntK*(TB-TA)*Double(0.5);
    IntG := IntG*(TB-TA)*Double(0.5);
    IntA := IntA*(TB-TA)*Double(0.5);
    State.Heap[J,0] := AbsReal(IntG-IntK);
    State.Heap[J,1] := IntK;
    State.Heap[J,2] := IntA;
    State.Heap[J,3] := TA;
    State.Heap[J,4] := TB;
    State.SumErr := State.SumErr+State.Heap[J,0];
    State.SumAbs := State.SumAbs+AbsReal(IntA);
    J := J+1;
    goto lbl_8;
lbl_10:
lbl_4:
    
    //
    // method iterations
    //
lbl_14:
    if False then
    begin
        goto lbl_15;
    end;
    
    //
    // additional memory if needed
    //
    if State.HeapUsed=State.HeapSize then
    begin
        MHeapResize(State.Heap, State.HeapSize, 4*State.HeapSize, State.HeapWidth);
    end;
    
    //
    // TODO: every 20 iterations recalculate errors/sums
    // TODO: one more criterion to prevent infinite loops with too strict Eps
    //
    if AP_FP_Less_Eq(State.SumErr,State.Eps*State.SumAbs) then
    begin
        State.R := 0;
        J:=0;
        while J<=State.HeapUsed-1 do
        begin
            State.R := State.R+State.Heap[J,1];
            Inc(J);
        end;
        Result := False;
        Exit;
    end;
    
    //
    // Exclude interval with maximum absolute error
    //
    MHeapPop(State.Heap, State.HeapUsed, State.HeapWidth);
    State.SumErr := State.SumErr-State.Heap[State.HeapUsed-1,0];
    State.SumAbs := State.SumAbs-State.Heap[State.HeapUsed-1,2];
    
    //
    // Divide interval, create subintervals
    //
    TA := State.Heap[State.HeapUsed-1,3];
    TB := State.Heap[State.HeapUsed-1,4];
    State.Heap[State.HeapUsed-1,3] := TA;
    State.Heap[State.HeapUsed-1,4] := Double(0.5)*(TA+TB);
    State.Heap[State.HeapUsed,3] := Double(0.5)*(TA+TB);
    State.Heap[State.HeapUsed,4] := TB;
    J := State.HeapUsed-1;
lbl_16:
    if J>State.HeapUsed then
    begin
        goto lbl_18;
    end;
    C1 := Double(0.5)*(State.Heap[J,4]-State.Heap[J,3]);
    C2 := Double(0.5)*(State.Heap[J,4]+State.Heap[J,3]);
    IntG := 0;
    IntK := 0;
    IntA := 0;
    I := 0;
lbl_19:
    if I>State.N-1 then
    begin
        goto lbl_21;
    end;
    
    //
    // F(x)
    //
    State.X := C1*State.QN[I]+C2;
    State.RState.Stage := 2;
    goto lbl_rcomm;
lbl_2:
    V := State.F;
    
    //
    // Gauss-Kronrod formula
    //
    IntK := IntK+V*State.WK[I];
    if I mod 2=1 then
    begin
        IntG := IntG+V*State.WG[I];
    end;
    
    //
    // Integral |F(x)|
    // Use rectangles method
    //
    IntA := IntA+AbsReal(V)*State.WR[I];
    I := I+1;
    goto lbl_19;
lbl_21:
    IntK := IntK*(State.Heap[J,4]-State.Heap[J,3])*Double(0.5);
    IntG := IntG*(State.Heap[J,4]-State.Heap[J,3])*Double(0.5);
    IntA := IntA*(State.Heap[J,4]-State.Heap[J,3])*Double(0.5);
    State.Heap[J,0] := AbsReal(IntG-IntK);
    State.Heap[J,1] := IntK;
    State.Heap[J,2] := IntA;
    State.SumErr := State.SumErr+State.Heap[J,0];
    State.SumAbs := State.SumAbs+State.Heap[J,2];
    J := J+1;
    goto lbl_16;
lbl_18:
    MHeapPush(State.Heap, State.HeapUsed-1, State.HeapWidth);
    MHeapPush(State.Heap, State.HeapUsed, State.HeapWidth);
    State.HeapUsed := State.HeapUsed+1;
    goto lbl_14;
lbl_15:
    Result := False;
    Exit;
    
    //
    // Saving state
    //
lbl_rcomm:
    Result := True;
    State.RState.IA[0] := I;
    State.RState.IA[1] := J;
    State.RState.IA[2] := NS;
    State.RState.IA[3] := Info;
    State.RState.RA[0] := C1;
    State.RState.RA[1] := C2;
    State.RState.RA[2] := IntG;
    State.RState.RA[3] := IntK;
    State.RState.RA[4] := IntA;
    State.RState.RA[5] := V;
    State.RState.RA[6] := TA;
    State.RState.RA[7] := TB;
    State.RState.RA[8] := QEps;
end;


procedure MHeapPop(var Heap : TReal2DArray;
     HeapSize : AlglibInteger;
     HeapWidth : AlglibInteger);
var
    I : AlglibInteger;
    P : AlglibInteger;
    T : Double;
    MaxCP : AlglibInteger;
begin
    if HeapSize=1 then
    begin
        Exit;
    end;
    I:=0;
    while I<=HeapWidth-1 do
    begin
        T := Heap[HeapSize-1,I];
        Heap[HeapSize-1,I] := Heap[0,I];
        Heap[0,I] := T;
        Inc(I);
    end;
    P := 0;
    while 2*P+1<HeapSize-1 do
    begin
        MaxCP := 2*P+1;
        if 2*P+2<HeapSize-1 then
        begin
            if AP_FP_Greater(Heap[2*P+2,0],Heap[2*P+1,0]) then
            begin
                MaxCP := 2*P+2;
            end;
        end;
        if AP_FP_Less(Heap[P,0],Heap[MaxCP,0]) then
        begin
            I:=0;
            while I<=HeapWidth-1 do
            begin
                T := Heap[P,I];
                Heap[P,I] := Heap[MaxCP,I];
                Heap[MaxCP,I] := T;
                Inc(I);
            end;
            P := MaxCP;
        end
        else
        begin
            Break;
        end;
    end;
end;


procedure MHeapPush(var Heap : TReal2DArray;
     HeapSize : AlglibInteger;
     HeapWidth : AlglibInteger);
var
    I : AlglibInteger;
    P : AlglibInteger;
    T : Double;
    Parent : AlglibInteger;
begin
    if HeapSize=0 then
    begin
        Exit;
    end;
    P := HeapSize;
    while P<>0 do
    begin
        Parent := (P-1) div 2;
        if AP_FP_Greater(Heap[P,0],Heap[Parent,0]) then
        begin
            I:=0;
            while I<=HeapWidth-1 do
            begin
                T := Heap[P,I];
                Heap[P,I] := Heap[Parent,I];
                Heap[Parent,I] := T;
                Inc(I);
            end;
            P := Parent;
        end
        else
        begin
            Break;
        end;
    end;
end;


procedure MHeapResize(var Heap : TReal2DArray;
     var HeapSize : AlglibInteger;
     NewHeapSize : AlglibInteger;
     HeapWidth : AlglibInteger);
var
    Tmp : TReal2DArray;
    I : AlglibInteger;
begin
    SetLength(Tmp, HeapSize, HeapWidth);
    I:=0;
    while I<=HeapSize-1 do
    begin
        APVMove(@Tmp[I][0], 0, HeapWidth-1, @Heap[I][0], 0, HeapWidth-1);
        Inc(I);
    end;
    SetLength(Heap, NewHeapSize, HeapWidth);
    I:=0;
    while I<=HeapSize-1 do
    begin
        APVMove(@Heap[I][0], 0, HeapWidth-1, @Tmp[I][0], 0, HeapWidth-1);
        Inc(I);
    end;
    HeapSize := NewHeapSize;
end;


end.