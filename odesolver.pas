{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Copyright 2009 by Sergey Bochkanov (ALGLIB project).

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
unit odesolver;
interface
uses Math, Sysutils, Ap;

type
ODESolverState = record
    N : AlglibInteger;
    M : AlglibInteger;
    XScale : Double;
    H : Double;
    Eps : Double;
    FracEps : Boolean;
    YC : TReal1DArray;
    EScale : TReal1DArray;
    XG : TReal1DArray;
    SolverType : AlglibInteger;
    X : Double;
    Y : TReal1DArray;
    DY : TReal1DArray;
    YTbl : TReal2DArray;
    RepTerminationType : AlglibInteger;
    RepNFEV : AlglibInteger;
    YN : TReal1DArray;
    YNS : TReal1DArray;
    RKA : TReal1DArray;
    RKC : TReal1DArray;
    RKCS : TReal1DArray;
    RKB : TReal2DArray;
    RKK : TReal2DArray;
    RState : RCommState;
end;


ODESolverReport = record
    NFEV : AlglibInteger;
    TerminationType : AlglibInteger;
end;



procedure ODESolverRKCK(const Y : TReal1DArray;
     N : AlglibInteger;
     const X : TReal1DArray;
     M : AlglibInteger;
     Eps : Double;
     H : Double;
     var State : ODESolverState);
function ODESolverIteration(var State : ODESolverState):Boolean;
procedure ODESolverResults(const State : ODESolverState;
     var M : AlglibInteger;
     var XTbl : TReal1DArray;
     var YTbl : TReal2DArray;
     var Rep : ODESolverReport);

implementation

const
    ODESolverMaxGrow = Double(3.0);
    ODESolverMaxShrink = Double(10.0);

procedure ODESolverInit(SolverType : AlglibInteger;
     const Y : TReal1DArray;
     N : AlglibInteger;
     const X : TReal1DArray;
     M : AlglibInteger;
     Eps : Double;
     H : Double;
     var State : ODESolverState);forward;


(*************************************************************************
Cash-Karp adaptive ODE solver.

This subroutine solves ODE  Y'=f(Y,x)  with  initial  conditions  Y(xs)=Ys
(here Y may be single variable or vector of N variables).

INPUT PARAMETERS:
    Y       -   initial conditions, array[0..N-1].
                contains values of Y[] at X[0]
    N       -   system size
    X       -   points at which Y should be tabulated, array[0..M-1]
                integrations starts at X[0], ends at X[M-1],  intermediate
                values at X[i] are returned too.
                SHOULD BE ORDERED BY ASCENDING OR BY DESCENDING!!!!
    M       -   number of intermediate points + first point + last point:
                * M>2 means that you need both Y(X[M-1]) and M-2 values at
                  intermediate points
                * M=2 means that you want just to integrate from  X[0]  to
                  X[1] and don't interested in intermediate values.
                * M=1 means that you don't want to integrate :)
                  it is degenerate case, but it will be handled correctly.
                * M<1 means error
    Eps     -   tolerance (absolute/relative error on each  step  will  be
                less than Eps). When passing:
                * Eps>0, it means desired ABSOLUTE error
                * Eps<0, it means desired RELATIVE error.  Relative errors
                  are calculated with respect to maximum values of  Y seen
                  so far. Be careful to use this criterion  when  starting
                  from Y[] that are close to zero.
    H       -   initial  step  lenth,  it  will  be adjusted automatically
                after the first  step.  If  H=0,  step  will  be  selected
                automatically  (usualy  it  will  be  equal  to  0.001  of
                min(x[i]-x[j])).

OUTPUT PARAMETERS
    State   -   structure which stores algorithm state between  subsequent
                calls of OdeSolverIteration. Used for reverse communication.
                This structure should be passed  to the OdeSolverIteration
                subroutine.

SEE ALSO
    AutoGKSmoothW, AutoGKSingular, AutoGKIteration, AutoGKResults.


  -- ALGLIB --
     Copyright 01.09.2009 by Bochkanov Sergey
*************************************************************************)
procedure ODESolverRKCK(const Y : TReal1DArray;
     N : AlglibInteger;
     const X : TReal1DArray;
     M : AlglibInteger;
     Eps : Double;
     H : Double;
     var State : ODESolverState);
begin
    ODESolverInit(0, Y, N, X, M, Eps, H, State);
end;


(*************************************************************************
One iteration of ODE solver.

Called after inialization of State structure with OdeSolverXXX subroutine.
See HTML docs for examples.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between subsequent
                calls and which is used for reverse communication. Must be
                initialized with OdeSolverXXX() call first.

If subroutine returned False, algorithm have finished its work.
If subroutine returned True, then user should:
* calculate F(State.X, State.Y)
* store it in State.DY
Here State.X is real, State.Y and State.DY are arrays[0..N-1] of reals.

  -- ALGLIB --
     Copyright 01.09.2009 by Bochkanov Sergey
*************************************************************************)
function ODESolverIteration(var State : ODESolverState):Boolean;
var
    N : AlglibInteger;
    M : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    XC : Double;
    V : Double;
    H : Double;
    H2 : Double;
    GridPoint : Boolean;
    Err : Double;
    MaxGrowPow : Double;
    KLimit : AlglibInteger;
label
lbl_3, lbl_6, lbl_8, lbl_0, lbl_10, lbl_7, lbl_5, lbl_1, lbl_rcomm;
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
        N := State.RState.IA[0];
        M := State.RState.IA[1];
        I := State.RState.IA[2];
        J := State.RState.IA[3];
        K := State.RState.IA[4];
        KLimit := State.RState.IA[5];
        GridPoint := State.RState.BA[0];
        XC := State.RState.RA[0];
        V := State.RState.RA[1];
        H := State.RState.RA[2];
        H2 := State.RState.RA[3];
        Err := State.RState.RA[4];
        MaxGrowPow := State.RState.RA[5];
    end
    else
    begin
        N := -983;
        M := -989;
        I := -834;
        J := 900;
        K := -287;
        KLimit := 364;
        GridPoint := False;
        XC := -338;
        V := -686;
        H := 912;
        H2 := 585;
        Err := 497;
        MaxGrowPow := -271;
    end;
    if State.RState.Stage=0 then
    begin
        goto lbl_0;
    end;
    
    //
    // Routine body
    //
    
    //
    // prepare
    //
    if State.RepTerminationType<>0 then
    begin
        Result := False;
        Exit;
    end;
    N := State.N;
    M := State.M;
    H := State.H;
    SetLength(State.Y, N);
    SetLength(State.DY, N);
    MaxGrowPow := Power(ODESolverMaxGrow, 5);
    State.RepNFEV := 0;
    
    //
    // some preliminary checks for internal errors
    // after this we assume that H>0 and M>1
    //
    Assert(AP_FP_Greater(State.H,0), 'ODESolver: internal error');
    Assert(M>1, 'ODESolverIteration: internal error');
    
    //
    // choose solver
    //
    if State.SolverType<>0 then
    begin
        goto lbl_1;
    end;
    
    //
    // Cask-Karp solver
    // Prepare coefficients table.
    // Check it for errors
    //
    SetLength(State.RKA, 6);
    State.RKA[0] := 0;
    State.RKA[1] := AP_Double(1)/5;
    State.RKA[2] := AP_Double(3)/10;
    State.RKA[3] := AP_Double(3)/5;
    State.RKA[4] := 1;
    State.RKA[5] := AP_Double(7)/8;
    SetLength(State.RKB, 6, 5);
    State.RKB[1,0] := AP_Double(1)/5;
    State.RKB[2,0] := AP_Double(3)/40;
    State.RKB[2,1] := AP_Double(9)/40;
    State.RKB[3,0] := AP_Double(3)/10;
    State.RKB[3,1] := -AP_Double(9)/10;
    State.RKB[3,2] := AP_Double(6)/5;
    State.RKB[4,0] := -AP_Double(11)/54;
    State.RKB[4,1] := AP_Double(5)/2;
    State.RKB[4,2] := -AP_Double(70)/27;
    State.RKB[4,3] := AP_Double(35)/27;
    State.RKB[5,0] := AP_Double(1631)/55296;
    State.RKB[5,1] := AP_Double(175)/512;
    State.RKB[5,2] := AP_Double(575)/13824;
    State.RKB[5,3] := AP_Double(44275)/110592;
    State.RKB[5,4] := AP_Double(253)/4096;
    SetLength(State.RKC, 6);
    State.RKC[0] := AP_Double(37)/378;
    State.RKC[1] := 0;
    State.RKC[2] := AP_Double(250)/621;
    State.RKC[3] := AP_Double(125)/594;
    State.RKC[4] := 0;
    State.RKC[5] := AP_Double(512)/1771;
    SetLength(State.RKCS, 6);
    State.RKCS[0] := AP_Double(2825)/27648;
    State.RKCS[1] := 0;
    State.RKCS[2] := AP_Double(18575)/48384;
    State.RKCS[3] := AP_Double(13525)/55296;
    State.RKCS[4] := AP_Double(277)/14336;
    State.RKCS[5] := AP_Double(1)/4;
    SetLength(State.RKK, 6, N);
    
    //
    // Main cycle consists of two iterations:
    // * outer where we travel from X[i-1] to X[i]
    // * inner where we travel inside [X[i-1],X[i]]
    //
    SetLength(State.YTbl, M, N);
    SetLength(State.EScale, N);
    SetLength(State.YN, N);
    SetLength(State.YNS, N);
    XC := State.XG[0];
    APVMove(@State.YTbl[0][0], 0, N-1, @State.YC[0], 0, N-1);
    J:=0;
    while J<=N-1 do
    begin
        State.EScale[J] := 0;
        Inc(J);
    end;
    I := 1;
lbl_3:
    if I>M-1 then
    begin
        goto lbl_5;
    end;
    
    //
    // begin inner iteration
    //
lbl_6:
    if False then
    begin
        goto lbl_7;
    end;
    
    //
    // truncate step if needed (beyond right boundary).
    // determine should we store X or not
    //
    if AP_FP_Greater_Eq(XC+H,State.XG[I]) then
    begin
        H := State.XG[I]-XC;
        GridPoint := True;
    end
    else
    begin
        GridPoint := False;
    end;
    
    //
    // Update error scale maximums
    //
    // These maximums are initialized by zeros,
    // then updated every iterations.
    //
    J:=0;
    while J<=N-1 do
    begin
        State.EScale[J] := Max(State.EScale[J], AbsReal(State.YC[J]));
        Inc(J);
    end;
    
    //
    // make one step:
    // 1. calculate all info needed to do step
    // 2. update errors scale maximums using values/derivatives
    //    obtained during (1)
    //
    // Take into account that we use scaling of X to reduce task
    // to the form where x[0] < x[1] < ... < x[n-1]. So X is
    // replaced by x=xscale*t, and dy/dx=f(y,x) is replaced
    // by dy/dt=xscale*f(y,xscale*t).
    //
    APVMove(@State.YN[0], 0, N-1, @State.YC[0], 0, N-1);
    APVMove(@State.YNS[0], 0, N-1, @State.YC[0], 0, N-1);
    K := 0;
lbl_8:
    if K>5 then
    begin
        goto lbl_10;
    end;
    
    //
    // prepare data for the next update of YN/YNS
    //
    State.X := State.XScale*(XC+State.RKA[K]*H);
    APVMove(@State.Y[0], 0, N-1, @State.YC[0], 0, N-1);
    J:=0;
    while J<=K-1 do
    begin
        V := State.RKB[K,J];
        APVAdd(@State.Y[0], 0, N-1, @State.RKK[J][0], 0, N-1, V);
        Inc(J);
    end;
    State.RState.Stage := 0;
    goto lbl_rcomm;
lbl_0:
    State.RepNFEV := State.RepNFEV+1;
    V := H*State.XScale;
    APVMove(@State.RKK[K][0], 0, N-1, @State.DY[0], 0, N-1, V);
    
    //
    // update YN/YNS
    //
    V := State.RKC[K];
    APVAdd(@State.YN[0], 0, N-1, @State.RKK[K][0], 0, N-1, V);
    V := State.RKCS[K];
    APVAdd(@State.YNS[0], 0, N-1, @State.RKK[K][0], 0, N-1, V);
    K := K+1;
    goto lbl_8;
lbl_10:
    
    //
    // estimate error
    //
    Err := 0;
    J:=0;
    while J<=N-1 do
    begin
        if  not State.FracEps then
        begin
            
            //
            // absolute error is estimated
            //
            Err := Max(Err, AbsReal(State.YN[J]-State.YNS[J]));
        end
        else
        begin
            
            //
            // Relative error is estimated
            //
            V := State.EScale[J];
            if AP_FP_Eq(V,0) then
            begin
                V := 1;
            end;
            Err := Max(Err, AbsReal(State.YN[J]-State.YNS[J])/V);
        end;
        Inc(J);
    end;
    
    //
    // calculate new step, restart if necessary
    //
    if AP_FP_Less_Eq(MaxGrowPow*Err,State.Eps) then
    begin
        H2 := ODESolverMaxGrow*H;
    end
    else
    begin
        H2 := H*Power(State.Eps/Err, Double(0.2));
    end;
    if AP_FP_Less(H2,H/ODESolverMaxShrink) then
    begin
        H2 := H/ODESolverMaxShrink;
    end;
    if AP_FP_Greater(Err,State.Eps) then
    begin
        H := H2;
        goto lbl_6;
    end;
    
    //
    // advance position
    //
    XC := XC+H;
    APVMove(@State.YC[0], 0, N-1, @State.YN[0], 0, N-1);
    
    //
    // update H
    //
    H := H2;
    
    //
    // break on grid point
    //
    if GridPoint then
    begin
        goto lbl_7;
    end;
    goto lbl_6;
lbl_7:
    
    //
    // save result
    //
    APVMove(@State.YTbl[I][0], 0, N-1, @State.YC[0], 0, N-1);
    I := I+1;
    goto lbl_3;
lbl_5:
    State.RepTerminationType := 1;
    Result := False;
    Exit;
lbl_1:
    Result := False;
    Exit;
    
    //
    // Saving state
    //
lbl_rcomm:
    Result := True;
    State.RState.IA[0] := N;
    State.RState.IA[1] := M;
    State.RState.IA[2] := I;
    State.RState.IA[3] := J;
    State.RState.IA[4] := K;
    State.RState.IA[5] := KLimit;
    State.RState.BA[0] := GridPoint;
    State.RState.RA[0] := XC;
    State.RState.RA[1] := V;
    State.RState.RA[2] := H;
    State.RState.RA[3] := H2;
    State.RState.RA[4] := Err;
    State.RState.RA[5] := MaxGrowPow;
end;


(*************************************************************************
ODE solver results

Called after OdeSolverIteration returned False.

INPUT PARAMETERS:
    State   -   algorithm state (used by OdeSolverIteration).

OUTPUT PARAMETERS:
    M       -   number of tabulated values, M>=1
    XTbl    -   array[0..M-1], values of X
    YTbl    -   array[0..M-1,0..N-1], values of Y in X[i]
    Rep     -   solver report:
                * Rep.TerminationType completetion code:
                    * -2    X is not ordered  by  ascending/descending  or
                            there are non-distinct X[],  i.e.  X[i]=X[i+1]
                    * -1    incorrect parameters were specified
                    *  1    task has been solved
                * Rep.NFEV contains number of function calculations

  -- ALGLIB --
     Copyright 01.09.2009 by Bochkanov Sergey
*************************************************************************)
procedure ODESolverResults(const State : ODESolverState;
     var M : AlglibInteger;
     var XTbl : TReal1DArray;
     var YTbl : TReal2DArray;
     var Rep : ODESolverReport);
var
    V : Double;
    I : AlglibInteger;
begin
    Rep.TerminationType := State.RepTerminationType;
    if Rep.TerminationType>0 then
    begin
        M := State.M;
        Rep.NFEV := State.RepNFEV;
        SetLength(XTbl, State.M);
        V := State.XScale;
        APVMove(@XTbl[0], 0, State.M-1, @State.XG[0], 0, State.M-1, V);
        SetLength(YTbl, State.M, State.N);
        I:=0;
        while I<=State.M-1 do
        begin
            APVMove(@YTbl[I][0], 0, State.N-1, @State.YTbl[I][0], 0, State.N-1);
            Inc(I);
        end;
    end
    else
    begin
        Rep.NFEV := 0;
    end;
end;


(*************************************************************************
Internal initialization subroutine
*************************************************************************)
procedure ODESolverInit(SolverType : AlglibInteger;
     const Y : TReal1DArray;
     N : AlglibInteger;
     const X : TReal1DArray;
     M : AlglibInteger;
     Eps : Double;
     H : Double;
     var State : ODESolverState);
var
    I : AlglibInteger;
    V : Double;
begin
    
    //
    // Prepare RComm
    //
    SetLength(State.RState.IA, 5+1);
    SetLength(State.RState.BA, 0+1);
    SetLength(State.RState.RA, 5+1);
    State.RState.Stage := -1;
    
    //
    // check parameters.
    //
    if (N<=0) or (M<1) or AP_FP_Eq(Eps,0) then
    begin
        State.RepTerminationType := -1;
        Exit;
    end;
    if AP_FP_Less(H,0) then
    begin
        H := -H;
    end;
    
    //
    // quick exit if necessary.
    // after this block we assume that M>1
    //
    if M=1 then
    begin
        State.RepNFEV := 0;
        State.RepTerminationType := 1;
        SetLength(State.YTbl, 1, N);
        APVMove(@State.YTbl[0][0], 0, N-1, @Y[0], 0, N-1);
        SetLength(State.XG, M);
        APVMove(@State.XG[0], 0, M-1, @X[0], 0, M-1);
        Exit;
    end;
    
    //
    // check again: correct order of X[]
    //
    if AP_FP_Eq(X[1],X[0]) then
    begin
        State.RepTerminationType := -2;
        Exit;
    end;
    I:=1;
    while I<=M-1 do
    begin
        if AP_FP_Greater(X[1],X[0]) and AP_FP_Less_Eq(X[I],X[I-1]) or AP_FP_Less(X[1],X[0]) and AP_FP_Greater_Eq(X[I],X[I-1]) then
        begin
            State.RepTerminationType := -2;
            Exit;
        end;
        Inc(I);
    end;
    
    //
    // auto-select H if necessary
    //
    if AP_FP_Eq(H,0) then
    begin
        V := AbsReal(X[1]-X[0]);
        I:=2;
        while I<=M-1 do
        begin
            V := Min(V, AbsReal(X[I]-X[I-1]));
            Inc(I);
        end;
        H := Double(0.001)*V;
    end;
    
    //
    // store parameters
    //
    State.N := N;
    State.M := M;
    State.H := H;
    State.Eps := AbsReal(Eps);
    State.FracEps := AP_FP_Less(Eps,0);
    SetLength(State.XG, M);
    APVMove(@State.XG[0], 0, M-1, @X[0], 0, M-1);
    if AP_FP_Greater(X[1],X[0]) then
    begin
        State.XScale := 1;
    end
    else
    begin
        State.XScale := -1;
        APVMul(@State.XG[0], 0, M-1, -1);
    end;
    SetLength(State.YC, N);
    APVMove(@State.YC[0], 0, N-1, @Y[0], 0, N-1);
    State.SolverType := SolverType;
    State.RepTerminationType := 0;
end;


end.