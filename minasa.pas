{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Copyright (c) 2010, Sergey Bochkanov (ALGLIB project).

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
unit minasa;
interface
uses Math, Sysutils, Ap, linmin;

type
MinASAState = record
    N : AlglibInteger;
    EpsG : Double;
    EpsF : Double;
    EpsX : Double;
    MaxIts : AlglibInteger;
    XRep : Boolean;
    StpMax : Double;
    CGType : AlglibInteger;
    K : AlglibInteger;
    NFEV : AlglibInteger;
    MCStage : AlglibInteger;
    BndL : TReal1DArray;
    BndU : TReal1DArray;
    CurAlgo : AlglibInteger;
    ACount : AlglibInteger;
    Mu : Double;
    FInit : Double;
    DGInit : Double;
    AK : TReal1DArray;
    XK : TReal1DArray;
    DK : TReal1DArray;
    AN : TReal1DArray;
    XN : TReal1DArray;
    DN : TReal1DArray;
    D : TReal1DArray;
    FOld : Double;
    Stp : Double;
    WORK : TReal1DArray;
    YK : TReal1DArray;
    GC : TReal1DArray;
    X : TReal1DArray;
    F : Double;
    G : TReal1DArray;
    NeedFG : Boolean;
    XUpdated : Boolean;
    RState : RCommState;
    RepIterationsCount : AlglibInteger;
    RepNFEV : AlglibInteger;
    RepTerminationType : AlglibInteger;
    DebugRestartsCount : AlglibInteger;
    LState : LINMINState;
    BetaHS : Double;
    BetaDY : Double;
end;


MinASAReport = record
    IterationsCount : AlglibInteger;
    NFEV : AlglibInteger;
    TerminationType : AlglibInteger;
    ActiveConstraints : AlglibInteger;
end;



procedure MinASACreate(N : AlglibInteger;
     const X : TReal1DArray;
     const BndL : TReal1DArray;
     const BndU : TReal1DArray;
     var State : MinASAState);
procedure MinASASetCond(var State : MinASAState;
     EpsG : Double;
     EpsF : Double;
     EpsX : Double;
     MaxIts : AlglibInteger);
procedure MinASASetXRep(var State : MinASAState; NeedXRep : Boolean);
procedure MinASASetAlgorithm(var State : MinASAState;
     AlgoType : AlglibInteger);
procedure MinASASetStpMax(var State : MinASAState; StpMax : Double);
function MinASAIteration(var State : MinASAState):Boolean;
procedure MinASAResults(const State : MinASAState;
     var X : TReal1DArray;
     var Rep : MinASAReport);

implementation

const
    N1 = 2;
    N2 = 2;
    STPMIN = Double(1.0E-300);
    GPAFTol = Double(0.0001);
    GPADecay = Double(0.5);
    ASARho = Double(0.5);

function ASABoundVal(X : Double; B1 : Double; B2 : Double):Double;forward;
function ASABoundedAntiGradNorm(const State : MinASAState):Double;forward;
function ASAGINorm(const State : MinASAState):Double;forward;
function ASAD1Norm(const State : MinASAState):Double;forward;
function ASAUIsEmpty(const State : MinASAState):Boolean;forward;
function ASAWantToUnstick(const State : MinASAState):Boolean;forward;
procedure ClearRequestFields(var State : MinASAState);forward;


(*************************************************************************
              NONLINEAR BOUND CONSTRAINED OPTIMIZATION USING
                               MODIFIED
                   WILLIAM W. HAGER AND HONGCHAO ZHANG
                         ACTIVE SET ALGORITHM

The  subroutine  minimizes  function  F(x)  of  N  arguments  with   bound
constraints: BndL[i] <= x[i] <= BndU[i]

This method is  globally  convergent  as  long  as  grad(f)  is  Lipschitz
continuous on a level set: L = { x : f(x)<=f(x0) }.

INPUT PARAMETERS:
    N       -   problem dimension. N>0
    X       -   initial solution approximation, array[0..N-1].
    BndL    -   lower bounds, array[0..N-1].
                all elements MUST be specified,  i.e.  all  variables  are
                bounded. However, if some (all) variables  are  unbounded,
                you may specify very small number as bound: -1000,  -1.0E6
                or -1.0E300, or something like that.
    BndU    -   upper bounds, array[0..N-1].
                all elements MUST be specified,  i.e.  all  variables  are
                bounded. However, if some (all) variables  are  unbounded,
                you may specify very large number as bound: +1000,  +1.0E6
                or +1.0E300, or something like that.
    EpsG    -   positive number which  defines  a  precision  of  search.  The
                subroutine finishes its work if the condition ||G|| < EpsG  is
                satisfied, where ||.|| means Euclidian norm, G - gradient, X -
                current approximation.
    EpsF    -   positive number which  defines  a  precision  of  search.  The
                subroutine finishes its work if on iteration  number  k+1  the
                condition |F(k+1)-F(k)| <= EpsF*max{|F(k)|, |F(k+1)|, 1}    is
                satisfied.
    EpsX    -   positive number which  defines  a  precision  of  search.  The
                subroutine finishes its work if on iteration number k+1    the
                condition |X(k+1)-X(k)| <= EpsX is fulfilled.
    MaxIts  -   maximum number of iterations. If MaxIts=0, the number of
                iterations is unlimited.

OUTPUT PARAMETERS:
    State - structure used for reverse communication.

This function  initializes  State   structure  with  default  optimization
parameters (stopping conditions, step size, etc.).  Use  MinASASet??????()
functions to tune optimization parameters.

After   all   optimization   parameters   are   tuned,   you   should  use
MinASAIteration() function to advance algorithm iterations.

NOTES:

1. you may tune stopping conditions with MinASASetCond() function
2. if target function contains exp() or other fast growing functions,  and
   optimization algorithm makes too large steps which leads  to  overflow,
   use MinASASetStpMax() function to bound algorithm's steps.

  -- ALGLIB --
     Copyright 25.03.2010 by Bochkanov Sergey
*************************************************************************)
procedure MinASACreate(N : AlglibInteger;
     const X : TReal1DArray;
     const BndL : TReal1DArray;
     const BndU : TReal1DArray;
     var State : MinASAState);
var
    I : AlglibInteger;
begin
    Assert(N>=1, 'MinASA: N too small!');
    I:=0;
    while I<=N-1 do
    begin
        Assert(AP_FP_Less_Eq(BndL[I],BndU[I]), 'MinASA: inconsistent bounds!');
        Assert(AP_FP_Less_Eq(BndL[I],X[I]), 'MinASA: infeasible X!');
        Assert(AP_FP_Less_Eq(X[I],BndU[I]), 'MinASA: infeasible X!');
        Inc(I);
    end;
    
    //
    // Initialize
    //
    State.N := N;
    MinASASetCond(State, 0, 0, 0, 0);
    MinASASetXRep(State, False);
    MinASASetStpMax(State, 0);
    MinASASetAlgorithm(State, -1);
    SetLength(State.BndL, N);
    SetLength(State.BndU, N);
    SetLength(State.AK, N);
    SetLength(State.XK, N);
    SetLength(State.DK, N);
    SetLength(State.AN, N);
    SetLength(State.XN, N);
    SetLength(State.DN, N);
    SetLength(State.X, N);
    SetLength(State.D, N);
    SetLength(State.G, N);
    SetLength(State.GC, N);
    SetLength(State.WORK, N);
    SetLength(State.YK, N);
    APVMove(@State.BndL[0], 0, N-1, @BndL[0], 0, N-1);
    APVMove(@State.BndU[0], 0, N-1, @BndU[0], 0, N-1);
    
    //
    // Prepare first run
    //
    APVMove(@State.X[0], 0, N-1, @X[0], 0, N-1);
    SetLength(State.RState.IA, 3+1);
    SetLength(State.RState.BA, 1+1);
    SetLength(State.RState.RA, 2+1);
    State.RState.Stage := -1;
end;


(*************************************************************************
This function sets stopping conditions for the ASA optimization algorithm.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be initialized
                with MinASACreate()
    EpsG    -   >=0
                The  subroutine  finishes  its  work   if   the  condition
                ||G||<EpsG is satisfied, where ||.|| means Euclidian norm,
                G - gradient.
    EpsF    -   >=0
                The  subroutine  finishes  its work if on k+1-th iteration
                the  condition  |F(k+1)-F(k)|<=EpsF*max{|F(k)|,|F(k+1)|,1}
                is satisfied.
    EpsX    -   >=0
                The subroutine finishes its work if  on  k+1-th  iteration
                the condition |X(k+1)-X(k)| <= EpsX is fulfilled.
    MaxIts  -   maximum number of iterations. If MaxIts=0, the  number  of
                iterations is unlimited.

Passing EpsG=0, EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to
automatic stopping criterion selection (small EpsX).

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************)
procedure MinASASetCond(var State : MinASAState;
     EpsG : Double;
     EpsF : Double;
     EpsX : Double;
     MaxIts : AlglibInteger);
begin
    Assert(AP_FP_Greater_Eq(EpsG,0), 'MinASASetCond: negative EpsG!');
    Assert(AP_FP_Greater_Eq(EpsF,0), 'MinASASetCond: negative EpsF!');
    Assert(AP_FP_Greater_Eq(EpsX,0), 'MinASASetCond: negative EpsX!');
    Assert(MaxIts>=0, 'MinASASetCond: negative MaxIts!');
    if AP_FP_Eq(EpsG,0) and AP_FP_Eq(EpsF,0) and AP_FP_Eq(EpsX,0) and (MaxIts=0) then
    begin
        EpsX := Double(1.0E-6);
    end;
    State.EpsG := EpsG;
    State.EpsF := EpsF;
    State.EpsX := EpsX;
    State.MaxIts := MaxIts;
end;


(*************************************************************************
This function turns on/off reporting.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be
                initialized with MinASACreate()
    NeedXRep-   whether iteration reports are needed or not

Usually  algorithm  returns from  MinASAIteration()  only  when  it  needs
function/gradient. However, with this function we can let  it  stop  after
each  iteration  (one  iteration  may  include   more  than  one  function
evaluation), which is indicated by XUpdated field.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************)
procedure MinASASetXRep(var State : MinASAState; NeedXRep : Boolean);
begin
    State.XRep := NeedXRep;
end;


(*************************************************************************
This function sets optimization algorithm.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be
                initialized with MinASACreate()
    UAType  -   algorithm type:
                * -1    automatic selection of the best algorithm
                * 0     DY (Dai and Yuan) algorithm
                * 1     Hybrid DY-HS algorithm

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************)
procedure MinASASetAlgorithm(var State : MinASAState;
     AlgoType : AlglibInteger);
begin
    Assert((AlgoType>=-1) and (AlgoType<=1), 'MinASASetAlgorithm: incorrect AlgoType!');
    if AlgoType=-1 then
    begin
        AlgoType := 1;
    end;
    State.CGType := AlgoType;
end;


(*************************************************************************
This function sets maximum step length

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be
                initialized with MinCGCreate()
    StpMax  -   maximum step length, >=0. Set StpMax to 0.0,  if you don't
                want to limit step length.

Use this subroutine when you optimize target function which contains exp()
or  other  fast  growing  functions,  and optimization algorithm makes too
large  steps  which  leads  to overflow. This function allows us to reject
steps  that  are  too  large  (and  therefore  expose  us  to the possible
overflow) without actually calculating function value at the x+stp*d.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************)
procedure MinASASetStpMax(var State : MinASAState; StpMax : Double);
begin
    Assert(AP_FP_Greater_Eq(StpMax,0), 'MinASASetStpMax: StpMax<0!');
    State.StpMax := StpMax;
end;


(*************************************************************************
One ASA iteration

Called after initialization with MinASACreate.
See HTML documentation for examples.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be initialized
                with MinASACreate.
RESULT:
* if function returned False, iterative proces has converged.
  Use MinLBFGSResults() to obtain optimization results.
* if subroutine returned True, then, depending on structure fields, we
  have one of the following situations


=== FUNC/GRAD REQUEST ===
State.NeedFG is True => function value/gradient are needed.
Caller should calculate function value State.F and gradient
State.G[0..N-1] at State.X[0..N-1] and call MinLBFGSIteration() again.

=== NEW INTERATION IS REPORTED ===
State.XUpdated is True => one more iteration was made.
State.X contains current position, State.F contains function value at X.
You can read info from these fields, but never modify  them  because  they
contain the only copy of optimization algorithm state.

One and only one of these fields (NeedFG, XUpdated) is true on return. New
iterations are reported only when reports  are  explicitly  turned  on  by
MinLBFGSSetXRep() function, so if you never called it, you can expect that
NeedFG is always True.


  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************)
function MinASAIteration(var State : MinASAState):Boolean;
var
    N : AlglibInteger;
    I : AlglibInteger;
    BetaK : Double;
    V : Double;
    VV : Double;
    MCINFO : AlglibInteger;
    B : Boolean;
    StepFound : Boolean;
    DiffCnt : AlglibInteger;
label
lbl_0, lbl_1, lbl_15, lbl_17, lbl_21, lbl_2, lbl_23, lbl_24, lbl_25, lbl_27, lbl_3, lbl_28, lbl_26, lbl_4, lbl_29, lbl_5, lbl_33, lbl_31, lbl_6, lbl_37, lbl_35, lbl_7, lbl_41, lbl_39, lbl_8, lbl_45, lbl_43, lbl_22, lbl_19, lbl_49, lbl_51, lbl_9, lbl_52, lbl_10, lbl_53, lbl_11, lbl_57, lbl_55, lbl_12, lbl_61, lbl_59, lbl_13, lbl_67, lbl_65, lbl_14, lbl_71, lbl_69, lbl_63, lbl_50, lbl_47, lbl_18, lbl_rcomm;
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
        I := State.RState.IA[1];
        MCINFO := State.RState.IA[2];
        DiffCnt := State.RState.IA[3];
        B := State.RState.BA[0];
        StepFound := State.RState.BA[1];
        BetaK := State.RState.RA[0];
        V := State.RState.RA[1];
        VV := State.RState.RA[2];
    end
    else
    begin
        N := -983;
        I := -989;
        MCINFO := -834;
        DiffCnt := 900;
        B := True;
        StepFound := False;
        BetaK := 214;
        V := -338;
        VV := -686;
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
    if State.RState.Stage=3 then
    begin
        goto lbl_3;
    end;
    if State.RState.Stage=4 then
    begin
        goto lbl_4;
    end;
    if State.RState.Stage=5 then
    begin
        goto lbl_5;
    end;
    if State.RState.Stage=6 then
    begin
        goto lbl_6;
    end;
    if State.RState.Stage=7 then
    begin
        goto lbl_7;
    end;
    if State.RState.Stage=8 then
    begin
        goto lbl_8;
    end;
    if State.RState.Stage=9 then
    begin
        goto lbl_9;
    end;
    if State.RState.Stage=10 then
    begin
        goto lbl_10;
    end;
    if State.RState.Stage=11 then
    begin
        goto lbl_11;
    end;
    if State.RState.Stage=12 then
    begin
        goto lbl_12;
    end;
    if State.RState.Stage=13 then
    begin
        goto lbl_13;
    end;
    if State.RState.Stage=14 then
    begin
        goto lbl_14;
    end;
    
    //
    // Routine body
    //
    
    //
    // Prepare
    //
    N := State.N;
    State.RepTerminationType := 0;
    State.RepIterationsCount := 0;
    State.RepNFEV := 0;
    State.DebugRestartsCount := 0;
    State.CGType := 1;
    APVMove(@State.XK[0], 0, N-1, @State.X[0], 0, N-1);
    I:=0;
    while I<=N-1 do
    begin
        if AP_FP_Eq(State.XK[I],State.BndL[I]) or AP_FP_Eq(State.XK[I],State.BndU[I]) then
        begin
            State.AK[I] := 0;
        end
        else
        begin
            State.AK[I] := 1;
        end;
        Inc(I);
    end;
    State.Mu := Double(0.1);
    State.CurAlgo := 0;
    
    //
    // Calculate F/G, initialize algorithm
    //
    ClearRequestFields(State);
    State.NeedFG := True;
    State.RState.Stage := 0;
    goto lbl_rcomm;
lbl_0:
    if not State.XRep then
    begin
        goto lbl_15;
    end;
    
    //
    // progress report
    //
    ClearRequestFields(State);
    State.XUpdated := True;
    State.RState.Stage := 1;
    goto lbl_rcomm;
lbl_1:
lbl_15:
    if AP_FP_Less_Eq(ASABoundedAntiGradNorm(State),State.EpsG) then
    begin
        State.RepTerminationType := 4;
        Result := False;
        Exit;
    end;
    State.RepNFEV := State.RepNFEV+1;
    
    //
    // Main cycle
    //
    // At the beginning of new iteration:
    // * CurAlgo stores current algorithm selector
    // * State.XK, State.F and State.G store current X/F/G
    // * State.AK stores current set of active constraints
    //
lbl_17:
    if False then
    begin
        goto lbl_18;
    end;
    
    //
    // GPA algorithm
    //
    if State.CurAlgo<>0 then
    begin
        goto lbl_19;
    end;
    State.K := 0;
    State.ACount := 0;
lbl_21:
    if False then
    begin
        goto lbl_22;
    end;
    
    //
    // Determine Dk = proj(xk - gk)-xk
    //
    I:=0;
    while I<=N-1 do
    begin
        State.D[I] := ASABoundVal(State.XK[I]-State.G[I], State.BndL[I], State.BndU[I])-State.XK[I];
        Inc(I);
    end;
    
    //
    // Armijo line search.
    // * exact search with alpha=1 is tried first,
    //   'exact' means that we evaluate f() EXACTLY at
    //   bound(x-g,bndl,bndu), without intermediate floating
    //   point operations.
    // * alpha<1 are tried if explicit search wasn't successful
    // Result is placed into XN.
    //
    // Two types of search are needed because we can't
    // just use second type with alpha=1 because in finite
    // precision arithmetics (x1-x0)+x0 may differ from x1.
    // So while x1 is correctly bounded (it lie EXACTLY on
    // boundary, if it is active), (x1-x0)+x0 may be
    // not bounded.
    //
    V := APVDotProduct(@State.D[0], 0, N-1, @State.G[0], 0, N-1);
    State.DGInit := V;
    State.FInit := State.F;
    if not(AP_FP_Less_Eq(ASAD1Norm(State),State.StpMax) or AP_FP_Eq(State.StpMax,0)) then
    begin
        goto lbl_23;
    end;
    
    //
    // Try alpha=1 step first
    //
    I:=0;
    while I<=N-1 do
    begin
        State.X[I] := ASABoundVal(State.XK[I]-State.G[I], State.BndL[I], State.BndU[I]);
        Inc(I);
    end;
    ClearRequestFields(State);
    State.NeedFG := True;
    State.RState.Stage := 2;
    goto lbl_rcomm;
lbl_2:
    State.RepNFEV := State.RepNFEV+1;
    StepFound := AP_FP_Less_Eq(State.F,State.FInit+GPAFTol*State.DGInit);
    goto lbl_24;
lbl_23:
    StepFound := False;
lbl_24:
    if not StepFound then
    begin
        goto lbl_25;
    end;
    
    //
    // we are at the boundary(ies)
    //
    APVMove(@State.XN[0], 0, N-1, @State.X[0], 0, N-1);
    State.Stp := 1;
    goto lbl_26;
lbl_25:
    
    //
    // alpha=1 is too large, try smaller values
    //
    State.Stp := 1;
    LinMinNormalizeD(State.D, State.Stp, N);
    State.DGInit := State.DGInit/State.Stp;
    State.Stp := GPADecay*State.Stp;
    if AP_FP_Greater(State.StpMax,0) then
    begin
        State.Stp := Min(State.Stp, State.StpMax);
    end;
lbl_27:
    if False then
    begin
        goto lbl_28;
    end;
    V := State.Stp;
    APVMove(@State.X[0], 0, N-1, @State.XK[0], 0, N-1);
    APVAdd(@State.X[0], 0, N-1, @State.D[0], 0, N-1, V);
    ClearRequestFields(State);
    State.NeedFG := True;
    State.RState.Stage := 3;
    goto lbl_rcomm;
lbl_3:
    State.RepNFEV := State.RepNFEV+1;
    if AP_FP_Less_Eq(State.Stp,STPMIN) then
    begin
        goto lbl_28;
    end;
    if AP_FP_Less_Eq(State.F,State.FInit+State.Stp*GPAFTol*State.DGInit) then
    begin
        goto lbl_28;
    end;
    State.Stp := State.Stp*GPADecay;
    goto lbl_27;
lbl_28:
    APVMove(@State.XN[0], 0, N-1, @State.X[0], 0, N-1);
lbl_26:
    State.RepIterationsCount := State.RepIterationsCount+1;
    if not State.XRep then
    begin
        goto lbl_29;
    end;
    
    //
    // progress report
    //
    ClearRequestFields(State);
    State.XUpdated := True;
    State.RState.Stage := 4;
    goto lbl_rcomm;
lbl_4:
lbl_29:
    
    //
    // Calculate new set of active constraints.
    // Reset counter if active set was changed.
    // Prepare for the new iteration
    //
    I:=0;
    while I<=N-1 do
    begin
        if AP_FP_Eq(State.XN[I],State.BndL[I]) or AP_FP_Eq(State.XN[I],State.BndU[I]) then
        begin
            State.AN[I] := 0;
        end
        else
        begin
            State.AN[I] := 1;
        end;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        if AP_FP_Neq(State.AK[I],State.AN[I]) then
        begin
            State.ACount := -1;
            Break;
        end;
        Inc(I);
    end;
    State.ACount := State.ACount+1;
    APVMove(@State.XK[0], 0, N-1, @State.XN[0], 0, N-1);
    APVMove(@State.AK[0], 0, N-1, @State.AN[0], 0, N-1);
    
    //
    // Stopping conditions
    //
    if not((State.RepIterationsCount>=State.MaxIts) and (State.MaxIts>0)) then
    begin
        goto lbl_31;
    end;
    
    //
    // Too many iterations
    //
    State.RepTerminationType := 5;
    if not State.XRep then
    begin
        goto lbl_33;
    end;
    ClearRequestFields(State);
    State.XUpdated := True;
    State.RState.Stage := 5;
    goto lbl_rcomm;
lbl_5:
lbl_33:
    Result := False;
    Exit;
lbl_31:
    if AP_FP_Greater(ASABoundedAntiGradNorm(State),State.EpsG) then
    begin
        goto lbl_35;
    end;
    
    //
    // Gradient is small enough
    //
    State.RepTerminationType := 4;
    if not State.XRep then
    begin
        goto lbl_37;
    end;
    ClearRequestFields(State);
    State.XUpdated := True;
    State.RState.Stage := 6;
    goto lbl_rcomm;
lbl_6:
lbl_37:
    Result := False;
    Exit;
lbl_35:
    V := APVDotProduct(@State.D[0], 0, N-1, @State.D[0], 0, N-1);
    if AP_FP_Greater(Sqrt(V)*State.Stp,State.EpsX) then
    begin
        goto lbl_39;
    end;
    
    //
    // Step size is too small, no further improvement is
    // possible
    //
    State.RepTerminationType := 2;
    if not State.XRep then
    begin
        goto lbl_41;
    end;
    ClearRequestFields(State);
    State.XUpdated := True;
    State.RState.Stage := 7;
    goto lbl_rcomm;
lbl_7:
lbl_41:
    Result := False;
    Exit;
lbl_39:
    if AP_FP_Greater(State.FInit-State.F,State.EpsF*Max(AbsReal(State.FInit), Max(AbsReal(State.F), Double(1.0)))) then
    begin
        goto lbl_43;
    end;
    
    //
    // F(k+1)-F(k) is small enough
    //
    State.RepTerminationType := 1;
    if not State.XRep then
    begin
        goto lbl_45;
    end;
    ClearRequestFields(State);
    State.XUpdated := True;
    State.RState.Stage := 8;
    goto lbl_rcomm;
lbl_8:
lbl_45:
    Result := False;
    Exit;
lbl_43:
    
    //
    // Decide - should we switch algorithm or not
    //
    if ASAUIsEmpty(State) then
    begin
        if AP_FP_Greater_Eq(ASAGINorm(State),State.Mu*ASAD1Norm(State)) then
        begin
            State.CurAlgo := 1;
            goto lbl_22;
        end
        else
        begin
            State.Mu := State.Mu*ASARho;
        end;
    end
    else
    begin
        if State.ACount=N1 then
        begin
            if AP_FP_Greater_Eq(ASAGINorm(State),State.Mu*ASAD1Norm(State)) then
            begin
                State.CurAlgo := 1;
                goto lbl_22;
            end;
        end;
    end;
    
    //
    // Next iteration
    //
    State.K := State.K+1;
    goto lbl_21;
lbl_22:
lbl_19:
    
    //
    // CG algorithm
    //
    if State.CurAlgo<>1 then
    begin
        goto lbl_47;
    end;
    
    //
    // first, check that there are non-active constraints.
    // move to GPA algorithm, if all constraints are active
    //
    B := True;
    I:=0;
    while I<=N-1 do
    begin
        if AP_FP_Neq(State.AK[I],0) then
        begin
            B := False;
            Break;
        end;
        Inc(I);
    end;
    if B then
    begin
        State.CurAlgo := 0;
        goto lbl_17;
    end;
    
    //
    // CG iterations
    //
    State.FOld := State.F;
    APVMove(@State.XK[0], 0, N-1, @State.X[0], 0, N-1);
    I:=0;
    while I<=N-1 do
    begin
        State.DK[I] := -State.G[I]*State.AK[I];
        State.GC[I] := State.G[I]*State.AK[I];
        Inc(I);
    end;
lbl_49:
    if False then
    begin
        goto lbl_50;
    end;
    
    //
    // Store G[k] for later calculation of Y[k]
    //
    I:=0;
    while I<=N-1 do
    begin
        State.YK[I] := -State.GC[I];
        Inc(I);
    end;
    
    //
    // Make a CG step in direction given by DK[]:
    // * calculate step. Step projection into feasible set
    //   is used. It has several benefits: a) step may be
    //   found with usual line search, b) multiple constraints
    //   may be activated with one step, c) activated constraints
    //   are detected in a natural way - just compare x[i] with
    //   bounds
    // * update active set, set B to True, if there
    //   were changes in the set.
    //
    APVMove(@State.D[0], 0, N-1, @State.DK[0], 0, N-1);
    APVMove(@State.XN[0], 0, N-1, @State.XK[0], 0, N-1);
    State.MCStage := 0;
    State.Stp := 1;
    LinMinNormalizeD(State.D, State.Stp, N);
    MCSRCH(N, State.XN, State.F, State.GC, State.D, State.Stp, State.StpMax, MCINFO, State.NFEV, State.WORK, State.LState, State.MCStage);
lbl_51:
    if State.MCStage=0 then
    begin
        goto lbl_52;
    end;
    
    //
    // preprocess data: bound State.XN so it belongs to the
    // feasible set and store it in the State.X
    //
    I:=0;
    while I<=N-1 do
    begin
        State.X[I] := ASABoundVal(State.XN[I], State.BndL[I], State.BndU[I]);
        Inc(I);
    end;
    
    //
    // RComm
    //
    ClearRequestFields(State);
    State.NeedFG := True;
    State.RState.Stage := 9;
    goto lbl_rcomm;
lbl_9:
    
    //
    // postprocess data: zero components of G corresponding to
    // the active constraints
    //
    I:=0;
    while I<=N-1 do
    begin
        if AP_FP_Eq(State.X[I],State.BndL[I]) or AP_FP_Eq(State.X[I],State.BndU[I]) then
        begin
            State.GC[I] := 0;
        end
        else
        begin
            State.GC[I] := State.G[I];
        end;
        Inc(I);
    end;
    MCSRCH(N, State.XN, State.F, State.GC, State.D, State.Stp, State.StpMax, MCINFO, State.NFEV, State.WORK, State.LState, State.MCStage);
    goto lbl_51;
lbl_52:
    DiffCnt := 0;
    I:=0;
    while I<=N-1 do
    begin
        
        //
        // XN contains unprojected result, project it,
        // save copy to X (will be used for progress reporting)
        //
        State.XN[I] := ASABoundVal(State.XN[I], State.BndL[I], State.BndU[I]);
        
        //
        // update active set
        //
        if AP_FP_Eq(State.XN[I],State.BndL[I]) or AP_FP_Eq(State.XN[I],State.BndU[I]) then
        begin
            State.AN[I] := 0;
        end
        else
        begin
            State.AN[I] := 1;
        end;
        if AP_FP_Neq(State.AN[I],State.AK[I]) then
        begin
            DiffCnt := DiffCnt+1;
        end;
        State.AK[I] := State.AN[I];
        Inc(I);
    end;
    APVMove(@State.XK[0], 0, N-1, @State.XN[0], 0, N-1);
    State.RepNFEV := State.RepNFEV+State.NFEV;
    State.RepIterationsCount := State.RepIterationsCount+1;
    if not State.XRep then
    begin
        goto lbl_53;
    end;
    
    //
    // progress report
    //
    ClearRequestFields(State);
    State.XUpdated := True;
    State.RState.Stage := 10;
    goto lbl_rcomm;
lbl_10:
lbl_53:
    
    //
    // Check stopping conditions.
    //
    if AP_FP_Greater(ASABoundedAntiGradNorm(State),State.EpsG) then
    begin
        goto lbl_55;
    end;
    
    //
    // Gradient is small enough
    //
    State.RepTerminationType := 4;
    if not State.XRep then
    begin
        goto lbl_57;
    end;
    ClearRequestFields(State);
    State.XUpdated := True;
    State.RState.Stage := 11;
    goto lbl_rcomm;
lbl_11:
lbl_57:
    Result := False;
    Exit;
lbl_55:
    if not((State.RepIterationsCount>=State.MaxIts) and (State.MaxIts>0)) then
    begin
        goto lbl_59;
    end;
    
    //
    // Too many iterations
    //
    State.RepTerminationType := 5;
    if not State.XRep then
    begin
        goto lbl_61;
    end;
    ClearRequestFields(State);
    State.XUpdated := True;
    State.RState.Stage := 12;
    goto lbl_rcomm;
lbl_12:
lbl_61:
    Result := False;
    Exit;
lbl_59:
    if not(AP_FP_Greater_Eq(ASAGINorm(State),State.Mu*ASAD1Norm(State)) and (DiffCnt=0)) then
    begin
        goto lbl_63;
    end;
    
    //
    // These conditions are explicitly or implicitly
    // related to the current step size and influenced
    // by changes in the active constraints.
    //
    // For these reasons they are checked only when we don't
    // want to 'unstick' at the end of the iteration and there
    // were no changes in the active set.
    //
    // NOTE: consition |G|>=Mu*|D1| must be exactly opposite
    // to the condition used to switch back to GPA. At least
    // one inequality must be strict, otherwise infinite cycle
    // may occur when |G|=Mu*|D1| (we DON'T test stopping
    // conditions and we DON'T switch to GPA, so we cycle
    // indefinitely).
    //
    if AP_FP_Greater(State.FOld-State.F,State.EpsF*Max(AbsReal(State.FOld), Max(AbsReal(State.F), Double(1.0)))) then
    begin
        goto lbl_65;
    end;
    
    //
    // F(k+1)-F(k) is small enough
    //
    State.RepTerminationType := 1;
    if not State.XRep then
    begin
        goto lbl_67;
    end;
    ClearRequestFields(State);
    State.XUpdated := True;
    State.RState.Stage := 13;
    goto lbl_rcomm;
lbl_13:
lbl_67:
    Result := False;
    Exit;
lbl_65:
    V := APVDotProduct(@State.D[0], 0, N-1, @State.D[0], 0, N-1);
    if AP_FP_Greater(Sqrt(V)*State.Stp,State.EpsX) then
    begin
        goto lbl_69;
    end;
    
    //
    // X(k+1)-X(k) is small enough
    //
    State.RepTerminationType := 2;
    if not State.XRep then
    begin
        goto lbl_71;
    end;
    ClearRequestFields(State);
    State.XUpdated := True;
    State.RState.Stage := 14;
    goto lbl_rcomm;
lbl_14:
lbl_71:
    Result := False;
    Exit;
lbl_69:
lbl_63:
    
    //
    // Check conditions for switching
    //
    if AP_FP_Less(ASAGINorm(State),State.Mu*ASAD1Norm(State)) then
    begin
        State.CurAlgo := 0;
        goto lbl_50;
    end;
    if DiffCnt>0 then
    begin
        if ASAUIsEmpty(State) or (DiffCnt>=N2) then
        begin
            State.CurAlgo := 1;
        end
        else
        begin
            State.CurAlgo := 0;
        end;
        goto lbl_50;
    end;
    
    //
    // Calculate D(k+1)
    //
    // Line search may result in:
    // * maximum feasible step being taken (already processed)
    // * point satisfying Wolfe conditions
    // * some kind of error (CG is restarted by assigning 0.0 to Beta)
    //
    if MCINFO=1 then
    begin
        
        //
        // Standard Wolfe conditions are satisfied:
        // * calculate Y[K] and BetaK
        //
        APVAdd(@State.YK[0], 0, N-1, @State.GC[0], 0, N-1);
        VV := APVDotProduct(@State.YK[0], 0, N-1, @State.DK[0], 0, N-1);
        V := APVDotProduct(@State.GC[0], 0, N-1, @State.GC[0], 0, N-1);
        State.BetaDY := V/VV;
        V := APVDotProduct(@State.GC[0], 0, N-1, @State.YK[0], 0, N-1);
        State.BetaHS := V/VV;
        if State.CGType=0 then
        begin
            BetaK := State.BetaDY;
        end;
        if State.CGType=1 then
        begin
            BetaK := Max(0, Min(State.BetaDY, State.BetaHS));
        end;
    end
    else
    begin
        
        //
        // Something is wrong (may be function is too wild or too flat).
        //
        // We'll set BetaK=0, which will restart CG algorithm.
        // We can stop later (during normal checks) if stopping conditions are met.
        //
        BetaK := 0;
        State.DebugRestartsCount := State.DebugRestartsCount+1;
    end;
    APVMoveNeg(@State.DN[0], 0, N-1, @State.GC[0], 0, N-1);
    APVAdd(@State.DN[0], 0, N-1, @State.DK[0], 0, N-1, BetaK);
    APVMove(@State.DK[0], 0, N-1, @State.DN[0], 0, N-1);
    
    //
    // update other information
    //
    State.FOld := State.F;
    State.K := State.K+1;
    goto lbl_49;
lbl_50:
lbl_47:
    goto lbl_17;
lbl_18:
    Result := False;
    Exit;
    
    //
    // Saving state
    //
lbl_rcomm:
    Result := True;
    State.RState.IA[0] := N;
    State.RState.IA[1] := I;
    State.RState.IA[2] := MCINFO;
    State.RState.IA[3] := DiffCnt;
    State.RState.BA[0] := B;
    State.RState.BA[1] := StepFound;
    State.RState.RA[0] := BetaK;
    State.RState.RA[1] := V;
    State.RState.RA[2] := VV;
end;


(*************************************************************************
Conjugate gradient results

Called after MinASA returned False.

INPUT PARAMETERS:
    State   -   algorithm state (used by MinASAIteration).

OUTPUT PARAMETERS:
    X       -   array[0..N-1], solution
    Rep     -   optimization report:
                * Rep.TerminationType completetion code:
                    * -2    rounding errors prevent further improvement.
                            X contains best point found.
                    * -1    incorrect parameters were specified
                    *  1    relative function improvement is no more than
                            EpsF.
                    *  2    relative step is no more than EpsX.
                    *  4    gradient norm is no more than EpsG
                    *  5    MaxIts steps was taken
                    *  7    stopping conditions are too stringent,
                            further improvement is impossible
                * Rep.IterationsCount contains iterations count
                * NFEV countains number of function calculations
                * ActiveConstraints contains number of active constraints

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************)
procedure MinASAResults(const State : MinASAState;
     var X : TReal1DArray;
     var Rep : MinASAReport);
var
    I : AlglibInteger;
begin
    SetLength(X, State.N-1+1);
    APVMove(@X[0], 0, State.N-1, @State.X[0], 0, State.N-1);
    Rep.IterationsCount := State.RepIterationsCount;
    Rep.NFEV := State.RepNFEV;
    Rep.TerminationType := State.RepTerminationType;
    Rep.ActiveConstraints := 0;
    I:=0;
    while I<=State.N-1 do
    begin
        if AP_FP_Eq(State.AK[I],0) then
        begin
            Rep.ActiveConstraints := Rep.ActiveConstraints+1;
        end;
        Inc(I);
    end;
end;


(*************************************************************************
'bound' value: map X to [B1,B2]

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************)
function ASABoundVal(X : Double; B1 : Double; B2 : Double):Double;
begin
    if AP_FP_Less_Eq(X,B1) then
    begin
        Result := B1;
        Exit;
    end;
    if AP_FP_Greater_Eq(X,B2) then
    begin
        Result := B2;
        Exit;
    end;
    Result := X;
end;


(*************************************************************************
Returns norm of bounded anti-gradient.

Bounded antigradient is a vector obtained from  anti-gradient  by  zeroing
components which point outwards:
    result = norm(v)
    v[i]=0     if ((-g[i]<0)and(x[i]=bndl[i])) or
                  ((-g[i]>0)and(x[i]=bndu[i]))
    v[i]=-g[i] otherwise

This function may be used to check a stopping criterion.

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************)
function ASABoundedAntiGradNorm(const State : MinASAState):Double;
var
    I : AlglibInteger;
    V : Double;
begin
    Result := 0;
    I:=0;
    while I<=State.N-1 do
    begin
        V := -State.G[I];
        if AP_FP_Eq(State.X[I],State.BndL[I]) and AP_FP_Less(-State.G[I],0) then
        begin
            V := 0;
        end;
        if AP_FP_Eq(State.X[I],State.BndU[I]) and AP_FP_Greater(-State.G[I],0) then
        begin
            V := 0;
        end;
        Result := Result+AP_Sqr(V);
        Inc(I);
    end;
    Result := Sqrt(Result);
end;


(*************************************************************************
Returns norm of GI(x).

GI(x) is  a  gradient  vector  whose  components  associated  with  active
constraints are zeroed. It  differs  from  bounded  anti-gradient  because
components  of   GI(x)   are   zeroed  independently  of  sign(g[i]),  and
anti-gradient's components are zeroed with respect to both constraint  and
sign.

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************)
function ASAGINorm(const State : MinASAState):Double;
var
    I : AlglibInteger;
    V : Double;
begin
    Result := 0;
    I:=0;
    while I<=State.N-1 do
    begin
        if AP_FP_Neq(State.X[I],State.BndL[I]) and AP_FP_Neq(State.X[I],State.BndU[I]) then
        begin
            Result := Result+AP_Sqr(State.G[I]);
        end;
        Inc(I);
    end;
    Result := Sqrt(Result);
end;


(*************************************************************************
Returns norm(D1(State.X))

For a meaning of D1 see 'NEW ACTIVE SET ALGORITHM FOR BOX CONSTRAINED
OPTIMIZATION' by WILLIAM W. HAGER AND HONGCHAO ZHANG.

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************)
function ASAD1Norm(const State : MinASAState):Double;
var
    I : AlglibInteger;
begin
    Result := 0;
    I:=0;
    while I<=State.N-1 do
    begin
        Result := Result+AP_Sqr(ASABoundVal(State.X[I]-State.G[I], State.BndL[I], State.BndU[I])-State.X[I]);
        Inc(I);
    end;
    Result := Sqrt(Result);
end;


(*************************************************************************
Returns True, if U set is empty.

* State.X is used as point,
* State.G - as gradient,
* D is calculated within function (because State.D may have different
  meaning depending on current optimization algorithm)

For a meaning of U see 'NEW ACTIVE SET ALGORITHM FOR BOX CONSTRAINED
OPTIMIZATION' by WILLIAM W. HAGER AND HONGCHAO ZHANG.

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************)
function ASAUIsEmpty(const State : MinASAState):Boolean;
var
    I : AlglibInteger;
    D : Double;
    D2 : Double;
    D32 : Double;
begin
    D := ASAD1Norm(State);
    D2 := Sqrt(D);
    D32 := D*D2;
    Result := True;
    I:=0;
    while I<=State.N-1 do
    begin
        if AP_FP_Greater_Eq(AbsReal(State.G[I]),D2) and AP_FP_Greater_Eq(Min(State.X[I]-State.BndL[I], State.BndU[I]-State.X[I]),D32) then
        begin
            Result := False;
            Exit;
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Returns True, if optimizer "want  to  unstick"  from  one  of  the  active
constraints, i.e. there is such active constraint with index I that either
lower bound is active and g[i]<0, or upper bound is active and g[i]>0.

State.X is used as current point, State.X - as gradient.
  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************)
function ASAWantToUnstick(const State : MinASAState):Boolean;
var
    I : AlglibInteger;
begin
    Result := False;
    I:=0;
    while I<=State.N-1 do
    begin
        if AP_FP_Eq(State.X[I],State.BndL[I]) and AP_FP_Less(State.G[I],0) then
        begin
            Result := True;
        end;
        if AP_FP_Eq(State.X[I],State.BndU[I]) and AP_FP_Greater(State.G[I],0) then
        begin
            Result := True;
        end;
        if Result then
        begin
            Exit;
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Clears request fileds (to be sure that we don't forgot to clear something)
*************************************************************************)
procedure ClearRequestFields(var State : MinASAState);
begin
    State.NeedFG := False;
    State.XUpdated := False;
end;


end.