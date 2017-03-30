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
unit mincg;
interface
uses Math, Sysutils, Ap, linmin;

type
MinCGState = record
    N : AlglibInteger;
    EpsG : Double;
    EpsF : Double;
    EpsX : Double;
    MaxIts : AlglibInteger;
    StpMax : Double;
    XRep : Boolean;
    CGType : AlglibInteger;
    NFEV : AlglibInteger;
    MCStage : AlglibInteger;
    K : AlglibInteger;
    XK : TReal1DArray;
    DK : TReal1DArray;
    XN : TReal1DArray;
    DN : TReal1DArray;
    D : TReal1DArray;
    FOld : Double;
    Stp : Double;
    WORK : TReal1DArray;
    YK : TReal1DArray;
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


MinCGReport = record
    IterationsCount : AlglibInteger;
    NFEV : AlglibInteger;
    TerminationType : AlglibInteger;
end;



procedure MinCGCreate(N : AlglibInteger;
     const X : TReal1DArray;
     var State : MinCGState);
procedure MinCGSetCond(var State : MinCGState;
     EpsG : Double;
     EpsF : Double;
     EpsX : Double;
     MaxIts : AlglibInteger);
procedure MinCGSetXRep(var State : MinCGState; NeedXRep : Boolean);
procedure MinCGSetCGType(var State : MinCGState; CGType : AlglibInteger);
procedure MinCGSetStpMax(var State : MinCGState; StpMax : Double);
function MinCGIteration(var State : MinCGState):Boolean;
procedure MinCGResults(const State : MinCGState;
     var X : TReal1DArray;
     var Rep : MinCGReport);

implementation

procedure ClearRequestFields(var State : MinCGState);forward;


(*************************************************************************
        NONLINEAR CONJUGATE GRADIENT METHOD

The subroutine minimizes function F(x) of N arguments by using one of  the
nonlinear conjugate gradient methods.

These CG methods are globally convergent (even on non-convex functions) as
long as grad(f) is Lipschitz continuous in  a  some  neighborhood  of  the
L = { x : f(x)<=f(x0) }.

INPUT PARAMETERS:
    N       -   problem dimension. N>0
    X       -   initial solution approximation, array[0..N-1].
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

See also MinCGIteration, MinCGResults

NOTE:

Passing EpsG=0, EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to
automatic stopping criterion selection (small EpsX).

  -- ALGLIB --
     Copyright 25.03.2010 by Bochkanov Sergey
*************************************************************************)
procedure MinCGCreate(N : AlglibInteger;
     const X : TReal1DArray;
     var State : MinCGState);
begin
    Assert(N>=1, 'MinCGCreate: N too small!');
    
    //
    // Initialize
    //
    State.N := N;
    MinCGSetCond(State, 0, 0, 0, 0);
    MinCGSetXRep(State, False);
    MinCGSetStpMax(State, 0);
    MinCGSetCGType(State, -1);
    SetLength(State.XK, N);
    SetLength(State.DK, N);
    SetLength(State.XN, N);
    SetLength(State.DN, N);
    SetLength(State.X, N);
    SetLength(State.D, N);
    SetLength(State.G, N);
    SetLength(State.WORK, N);
    SetLength(State.YK, N);
    
    //
    // Prepare first run
    //
    APVMove(@State.X[0], 0, N-1, @X[0], 0, N-1);
    SetLength(State.RState.IA, 2+1);
    SetLength(State.RState.RA, 2+1);
    State.RState.Stage := -1;
end;


(*************************************************************************
This function sets stopping conditions for CG optimization algorithm.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be initialized
                with MinCGCreate()
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
procedure MinCGSetCond(var State : MinCGState;
     EpsG : Double;
     EpsF : Double;
     EpsX : Double;
     MaxIts : AlglibInteger);
begin
    Assert(AP_FP_Greater_Eq(EpsG,0), 'MinCGSetCond: negative EpsG!');
    Assert(AP_FP_Greater_Eq(EpsF,0), 'MinCGSetCond: negative EpsF!');
    Assert(AP_FP_Greater_Eq(EpsX,0), 'MinCGSetCond: negative EpsX!');
    Assert(MaxIts>=0, 'MinCGSetCond: negative MaxIts!');
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
                initialized with MinCGCreate()
    NeedXRep-   whether iteration reports are needed or not

Usually  algorithm  returns  from  MinCGIteration()  only  when  it  needs
function/gradient. However, with this function we can let  it  stop  after
each  iteration  (one  iteration  may  include   more  than  one  function
evaluation), which is indicated by XUpdated field.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************)
procedure MinCGSetXRep(var State : MinCGState; NeedXRep : Boolean);
begin
    State.XRep := NeedXRep;
end;


(*************************************************************************
This function sets CG algorithm.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be
                initialized with MinCGCreate()
    CGType  -   algorithm type:
                * -1    automatic selection of the best algorithm
                * 0     DY (Dai and Yuan) algorithm
                * 1     Hybrid DY-HS algorithm

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************)
procedure MinCGSetCGType(var State : MinCGState; CGType : AlglibInteger);
begin
    Assert((CGType>=-1) and (CGType<=1), 'MinCGSetCGType: incorrect CGType!');
    if CGType=-1 then
    begin
        CGType := 1;
    end;
    State.CGType := CGType;
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
procedure MinCGSetStpMax(var State : MinCGState; StpMax : Double);
begin
    Assert(AP_FP_Greater_Eq(StpMax,0), 'MinCGSetStpMax: StpMax<0!');
    State.StpMax := StpMax;
end;


(*************************************************************************
One conjugate gradient iteration

Called after initialization with MinCG.
See HTML documentation for examples.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be initialized
                with MinCG.

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
     Copyright 20.04.2009 by Bochkanov Sergey
*************************************************************************)
function MinCGIteration(var State : MinCGState):Boolean;
var
    N : AlglibInteger;
    I : AlglibInteger;
    BetaK : Double;
    V : Double;
    VV : Double;
    MCINFO : AlglibInteger;
label
lbl_0, lbl_1, lbl_4, lbl_6, lbl_8, lbl_2, lbl_9, lbl_3, lbl_10, lbl_7, lbl_rcomm;
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
        BetaK := State.RState.RA[0];
        V := State.RState.RA[1];
        VV := State.RState.RA[2];
    end
    else
    begin
        N := -983;
        I := -989;
        MCINFO := -834;
        BetaK := 900;
        V := -287;
        VV := 364;
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
        goto lbl_4;
    end;
    ClearRequestFields(State);
    State.XUpdated := True;
    State.RState.Stage := 1;
    goto lbl_rcomm;
lbl_1:
lbl_4:
    V := APVDotProduct(@State.G[0], 0, N-1, @State.G[0], 0, N-1);
    V := Sqrt(V);
    if AP_FP_Eq(V,0) then
    begin
        State.RepTerminationType := 4;
        Result := False;
        Exit;
    end;
    State.RepNFEV := 1;
    State.K := 0;
    State.FOld := State.F;
    APVMove(@State.XK[0], 0, N-1, @State.X[0], 0, N-1);
    APVMoveNeg(@State.DK[0], 0, N-1, @State.G[0], 0, N-1);
    
    //
    // Main cycle
    //
lbl_6:
    if False then
    begin
        goto lbl_7;
    end;
    
    //
    // Store G[k] for later calculation of Y[k]
    //
    APVMoveNeg(@State.YK[0], 0, N-1, @State.G[0], 0, N-1);
    
    //
    // Calculate X(k+1): minimize F(x+alpha*d)
    //
    APVMove(@State.D[0], 0, N-1, @State.DK[0], 0, N-1);
    APVMove(@State.X[0], 0, N-1, @State.XK[0], 0, N-1);
    State.MCStage := 0;
    State.Stp := Double(1.0);
    LinMinNormalizeD(State.D, State.Stp, N);
    MCSRCH(N, State.X, State.F, State.G, State.D, State.Stp, State.StpMax, MCINFO, State.NFEV, State.WORK, State.LState, State.MCStage);
lbl_8:
    if State.MCStage=0 then
    begin
        goto lbl_9;
    end;
    ClearRequestFields(State);
    State.NeedFG := True;
    State.RState.Stage := 2;
    goto lbl_rcomm;
lbl_2:
    MCSRCH(N, State.X, State.F, State.G, State.D, State.Stp, State.StpMax, MCINFO, State.NFEV, State.WORK, State.LState, State.MCStage);
    goto lbl_8;
lbl_9:
    if not State.XRep then
    begin
        goto lbl_10;
    end;
    ClearRequestFields(State);
    State.XUpdated := True;
    State.RState.Stage := 3;
    goto lbl_rcomm;
lbl_3:
lbl_10:
    APVMove(@State.XN[0], 0, N-1, @State.X[0], 0, N-1);
    if MCINFO=1 then
    begin
        
        //
        // Standard Wolfe conditions hold
        // Calculate Y[K] and BetaK
        //
        APVAdd(@State.YK[0], 0, N-1, @State.G[0], 0, N-1);
        VV := APVDotProduct(@State.YK[0], 0, N-1, @State.DK[0], 0, N-1);
        V := APVDotProduct(@State.G[0], 0, N-1, @State.G[0], 0, N-1);
        State.BetaDY := V/VV;
        V := APVDotProduct(@State.G[0], 0, N-1, @State.YK[0], 0, N-1);
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
    
    //
    // Calculate D(k+1)
    //
    APVMoveNeg(@State.DN[0], 0, N-1, @State.G[0], 0, N-1);
    APVAdd(@State.DN[0], 0, N-1, @State.DK[0], 0, N-1, BetaK);
    
    //
    // Update information and Hessian.
    // Check stopping conditions.
    //
    State.RepNFEV := State.RepNFEV+State.NFEV;
    State.RepIterationsCount := State.RepIterationsCount+1;
    if (State.RepIterationsCount>=State.MaxIts) and (State.MaxIts>0) then
    begin
        
        //
        // Too many iterations
        //
        State.RepTerminationType := 5;
        Result := False;
        Exit;
    end;
    V := APVDotProduct(@State.G[0], 0, N-1, @State.G[0], 0, N-1);
    if AP_FP_Less_Eq(Sqrt(V),State.EpsG) then
    begin
        
        //
        // Gradient is small enough
        //
        State.RepTerminationType := 4;
        Result := False;
        Exit;
    end;
    if AP_FP_Less_Eq(State.FOld-State.F,State.EpsF*Max(AbsReal(State.FOld), Max(AbsReal(State.F), Double(1.0)))) then
    begin
        
        //
        // F(k+1)-F(k) is small enough
        //
        State.RepTerminationType := 1;
        Result := False;
        Exit;
    end;
    V := APVDotProduct(@State.D[0], 0, N-1, @State.D[0], 0, N-1);
    if AP_FP_Less_Eq(Sqrt(V)*State.Stp,State.EpsX) then
    begin
        
        //
        // X(k+1)-X(k) is small enough
        //
        State.RepTerminationType := 2;
        Result := False;
        Exit;
    end;
    
    //
    // Shift Xk/Dk, update other information
    //
    APVMove(@State.XK[0], 0, N-1, @State.XN[0], 0, N-1);
    APVMove(@State.DK[0], 0, N-1, @State.DN[0], 0, N-1);
    State.FOld := State.F;
    State.K := State.K+1;
    goto lbl_6;
lbl_7:
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
    State.RState.RA[0] := BetaK;
    State.RState.RA[1] := V;
    State.RState.RA[2] := VV;
end;


(*************************************************************************
Conjugate gradient results

Called after MinCG returned False.

INPUT PARAMETERS:
    State   -   algorithm state (used by MinCGIteration).

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

  -- ALGLIB --
     Copyright 20.04.2009 by Bochkanov Sergey
*************************************************************************)
procedure MinCGResults(const State : MinCGState;
     var X : TReal1DArray;
     var Rep : MinCGReport);
begin
    SetLength(X, State.N-1+1);
    APVMove(@X[0], 0, State.N-1, @State.XN[0], 0, State.N-1);
    Rep.IterationsCount := State.RepIterationsCount;
    Rep.NFEV := State.RepNFEV;
    Rep.TerminationType := State.RepTerminationType;
end;


(*************************************************************************
Clears request fileds (to be sure that we don't forgot to clear something)
*************************************************************************)
procedure ClearRequestFields(var State : MinCGState);
begin
    State.NeedFG := False;
    State.XUpdated := False;
end;


end.