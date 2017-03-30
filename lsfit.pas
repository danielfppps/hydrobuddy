{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Copyright (c) 2006-2009, Sergey Bochkanov (ALGLIB project).

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
unit lsfit;
interface
uses Math, Sysutils, Ap, blas, reflections, creflections, hqrnd, matgen, ablasf, ablas, trfac, trlinsolve, safesolve, rcond, matinv, hblas, sblas, ortfac, rotations, bdsvd, svd, xblas, densesolver, linmin, minlbfgs, minlm;

type
(*************************************************************************
Least squares fitting report:
    TaskRCond       reciprocal of task's condition number
    RMSError        RMS error
    AvgError        average error
    AvgRelError     average relative error (for non-zero Y[I])
    MaxError        maximum error
*************************************************************************)
LSFitReport = record
    TaskRCond : Double;
    RMSError : Double;
    AvgError : Double;
    AvgRelError : Double;
    MaxError : Double;
end;


LSFitState = record
    N : AlglibInteger;
    M : AlglibInteger;
    K : AlglibInteger;
    EpsF : Double;
    EpsX : Double;
    MaxIts : AlglibInteger;
    StpMax : Double;
    TaskX : TReal2DArray;
    TaskY : TReal1DArray;
    W : TReal1DArray;
    CheapFG : Boolean;
    HaveHess : Boolean;
    NeedF : Boolean;
    NeedFG : Boolean;
    NeedFGH : Boolean;
    PointIndex : AlglibInteger;
    X : TReal1DArray;
    C : TReal1DArray;
    F : Double;
    G : TReal1DArray;
    H : TReal2DArray;
    RepTerminationType : AlglibInteger;
    RepRMSError : Double;
    RepAvgError : Double;
    RepAvgRelError : Double;
    RepMaxError : Double;
    OptState : MinLMState;
    OptRep : MinLMReport;
    RState : RCommState;
end;



procedure LSFitLinearW(const Y : TReal1DArray;
     const W : TReal1DArray;
     const FMatrix : TReal2DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var C : TReal1DArray;
     var Rep : LSFitReport);
procedure LSFitLinearWC(Y : TReal1DArray;
     const W : TReal1DArray;
     const FMatrix : TReal2DArray;
     CMatrix : TReal2DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     K : AlglibInteger;
     var Info : AlglibInteger;
     var C : TReal1DArray;
     var Rep : LSFitReport);
procedure LSFitLinear(const Y : TReal1DArray;
     const FMatrix : TReal2DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var C : TReal1DArray;
     var Rep : LSFitReport);
procedure LSFitLinearC(Y : TReal1DArray;
     const FMatrix : TReal2DArray;
     const CMatrix : TReal2DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     K : AlglibInteger;
     var Info : AlglibInteger;
     var C : TReal1DArray;
     var Rep : LSFitReport);
procedure LSFitNonlinearWFG(const X : TReal2DArray;
     const Y : TReal1DArray;
     const W : TReal1DArray;
     const C : TReal1DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     K : AlglibInteger;
     CheapFG : Boolean;
     var State : LSFitState);
procedure LSFitNonlinearFG(const X : TReal2DArray;
     const Y : TReal1DArray;
     const C : TReal1DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     K : AlglibInteger;
     CheapFG : Boolean;
     var State : LSFitState);
procedure LSFitNonlinearWFGH(const X : TReal2DArray;
     const Y : TReal1DArray;
     const W : TReal1DArray;
     const C : TReal1DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     K : AlglibInteger;
     var State : LSFitState);
procedure LSFitNonlinearFGH(const X : TReal2DArray;
     const Y : TReal1DArray;
     const C : TReal1DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     K : AlglibInteger;
     var State : LSFitState);
procedure LSFitNonlinearSetCond(var State : LSFitState;
     EpsF : Double;
     EpsX : Double;
     MaxIts : AlglibInteger);
procedure LSFitNonlinearSetStpMax(var State : LSFitState; StpMax : Double);
function LSFitNonlinearIteration(var State : LSFitState):Boolean;
procedure LSFitNonlinearResults(const State : LSFitState;
     var Info : AlglibInteger;
     var C : TReal1DArray;
     var Rep : LSFitReport);
procedure LSFitScaleXY(var X : TReal1DArray;
     var Y : TReal1DArray;
     N : AlglibInteger;
     var XC : TReal1DArray;
     var YC : TReal1DArray;
     const DC : TInteger1DArray;
     K : AlglibInteger;
     var XA : Double;
     var XB : Double;
     var SA : Double;
     var SB : Double;
     var XOriginal : TReal1DArray;
     var YOriginal : TReal1DArray);

implementation

procedure LSFitLinearInternal(const Y : TReal1DArray;
     const W : TReal1DArray;
     const FMatrix : TReal2DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var C : TReal1DArray;
     var Rep : LSFitReport);forward;
procedure LSFitClearRequestFields(var State : LSFitState);forward;


(*************************************************************************
Weighted linear least squares fitting.

QR decomposition is used to reduce task to MxM, then triangular solver  or
SVD-based solver is used depending on condition number of the  system.  It
allows to maximize speed and retain decent accuracy.

INPUT PARAMETERS:
    Y       -   array[0..N-1] Function values in  N  points.
    W       -   array[0..N-1]  Weights  corresponding to function  values.
                Each summand in square  sum  of  approximation  deviations
                from  given  values  is  multiplied  by  the   square   of
                corresponding weight.
    FMatrix -   a table of basis functions values, array[0..N-1, 0..M-1].
                FMatrix[I, J] - value of J-th basis function in I-th point.
    N       -   number of points used. N>=1.
    M       -   number of basis functions, M>=1.

OUTPUT PARAMETERS:
    Info    -   error code:
                * -4    internal SVD decomposition subroutine failed (very
                        rare and for degenerate systems only)
                * -1    incorrect N/M were specified
                *  1    task is solved
    C       -   decomposition coefficients, array[0..M-1]
    Rep     -   fitting report. Following fields are set:
                * Rep.TaskRCond     reciprocal of condition number
                * RMSError          rms error on the (X,Y).
                * AvgError          average error on the (X,Y).
                * AvgRelError       average relative error on the non-zero Y
                * MaxError          maximum error
                                    NON-WEIGHTED ERRORS ARE CALCULATED

SEE ALSO
    LSFitLinear
    LSFitLinearC
    LSFitLinearWC

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure LSFitLinearW(const Y : TReal1DArray;
     const W : TReal1DArray;
     const FMatrix : TReal2DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var C : TReal1DArray;
     var Rep : LSFitReport);
begin
    LSFitLinearInternal(Y, W, FMatrix, N, M, Info, C, Rep);
end;


(*************************************************************************
Weighted constained linear least squares fitting.

This  is  variation  of LSFitLinearW(), which searchs for min|A*x=b| given
that  K  additional  constaints  C*x=bc are satisfied. It reduces original
task to modified one: min|B*y-d| WITHOUT constraints,  then LSFitLinearW()
is called.

INPUT PARAMETERS:
    Y       -   array[0..N-1] Function values in  N  points.
    W       -   array[0..N-1]  Weights  corresponding to function  values.
                Each summand in square  sum  of  approximation  deviations
                from  given  values  is  multiplied  by  the   square   of
                corresponding weight.
    FMatrix -   a table of basis functions values, array[0..N-1, 0..M-1].
                FMatrix[I,J] - value of J-th basis function in I-th point.
    CMatrix -   a table of constaints, array[0..K-1,0..M].
                I-th row of CMatrix corresponds to I-th linear constraint:
                CMatrix[I,0]*C[0] + ... + CMatrix[I,M-1]*C[M-1] = CMatrix[I,M]
    N       -   number of points used. N>=1.
    M       -   number of basis functions, M>=1.
    K       -   number of constraints, 0 <= K < M
                K=0 corresponds to absence of constraints.

OUTPUT PARAMETERS:
    Info    -   error code:
                * -4    internal SVD decomposition subroutine failed (very
                        rare and for degenerate systems only)
                * -3    either   too   many  constraints  (M   or   more),
                        degenerate  constraints   (some   constraints  are
                        repetead twice) or inconsistent  constraints  were
                        specified.
                * -1    incorrect N/M/K were specified
                *  1    task is solved
    C       -   decomposition coefficients, array[0..M-1]
    Rep     -   fitting report. Following fields are set:
                * RMSError          rms error on the (X,Y).
                * AvgError          average error on the (X,Y).
                * AvgRelError       average relative error on the non-zero Y
                * MaxError          maximum error
                                    NON-WEIGHTED ERRORS ARE CALCULATED

IMPORTANT:
    this subroitine doesn't calculate task's condition number for K<>0.

SEE ALSO
    LSFitLinear
    LSFitLinearC
    LSFitLinearWC

  -- ALGLIB --
     Copyright 07.09.2009 by Bochkanov Sergey
*************************************************************************)
procedure LSFitLinearWC(Y : TReal1DArray;
     const W : TReal1DArray;
     const FMatrix : TReal2DArray;
     CMatrix : TReal2DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     K : AlglibInteger;
     var Info : AlglibInteger;
     var C : TReal1DArray;
     var Rep : LSFitReport);
var
    I : AlglibInteger;
    J : AlglibInteger;
    Tau : TReal1DArray;
    Q : TReal2DArray;
    F2 : TReal2DArray;
    Tmp : TReal1DArray;
    C0 : TReal1DArray;
    V : Double;
begin
    Y := DynamicArrayCopy(Y);
    CMatrix := DynamicArrayCopy(CMatrix);
    if (N<1) or (M<1) or (K<0) then
    begin
        Info := -1;
        Exit;
    end;
    if K>=M then
    begin
        Info := -3;
        Exit;
    end;
    
    //
    // Solve
    //
    if K=0 then
    begin
        
        //
        // no constraints
        //
        LSFitLinearInternal(Y, W, FMatrix, N, M, Info, C, Rep);
    end
    else
    begin
        
        //
        // First, find general form solution of constraints system:
        // * factorize C = L*Q
        // * unpack Q
        // * fill upper part of C with zeros (for RCond)
        //
        // We got C=C0+Q2'*y where Q2 is lower M-K rows of Q.
        //
        RMatrixLQ(CMatrix, K, M, Tau);
        RMatrixLQUnpackQ(CMatrix, K, M, Tau, M, Q);
        I:=0;
        while I<=K-1 do
        begin
            J:=I+1;
            while J<=M-1 do
            begin
                CMatrix[I,J] := Double(0.0);
                Inc(J);
            end;
            Inc(I);
        end;
        if AP_FP_Less(RMatrixLURCondInf(CMatrix, K),1000*MachineEpsilon) then
        begin
            Info := -3;
            Exit;
        end;
        SetLength(Tmp, K);
        I:=0;
        while I<=K-1 do
        begin
            if I>0 then
            begin
                V := APVDotProduct(@CMatrix[I][0], 0, I-1, @Tmp[0], 0, I-1);
            end
            else
            begin
                V := 0;
            end;
            Tmp[I] := (CMatrix[I,M]-V)/CMatrix[I,I];
            Inc(I);
        end;
        SetLength(C0, M);
        I:=0;
        while I<=M-1 do
        begin
            C0[I] := 0;
            Inc(I);
        end;
        I:=0;
        while I<=K-1 do
        begin
            V := Tmp[I];
            APVAdd(@C0[0], 0, M-1, @Q[I][0], 0, M-1, V);
            Inc(I);
        end;
        
        //
        // Second, prepare modified matrix F2 = F*Q2' and solve modified task
        //
        SetLength(Tmp, Max(N, M)+1);
        SetLength(F2, N, M-K);
        MatrixVectorMultiply(FMatrix, 0, N-1, 0, M-1, False, C0, 0, M-1, -Double(1.0), Y, 0, N-1, Double(1.0));
        MatrixMatrixMultiply(FMatrix, 0, N-1, 0, M-1, False, Q, K, M-1, 0, M-1, True, Double(1.0), F2, 0, N-1, 0, M-K-1, Double(0.0), Tmp);
        LSFitLinearInternal(Y, W, F2, N, M-K, Info, Tmp, Rep);
        Rep.TaskRCond := -1;
        if Info<=0 then
        begin
            Exit;
        end;
        
        //
        // then, convert back to original answer: C = C0 + Q2'*Y0
        //
        SetLength(C, M);
        APVMove(@C[0], 0, M-1, @C0[0], 0, M-1);
        MatrixVectorMultiply(Q, K, M-1, 0, M-1, True, Tmp, 0, M-K-1, Double(1.0), C, 0, M-1, Double(1.0));
    end;
end;


(*************************************************************************
Linear least squares fitting, without weights.

See LSFitLinearW for more information.

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure LSFitLinear(const Y : TReal1DArray;
     const FMatrix : TReal2DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var C : TReal1DArray;
     var Rep : LSFitReport);
var
    W : TReal1DArray;
    I : AlglibInteger;
begin
    if N<1 then
    begin
        Info := -1;
        Exit;
    end;
    SetLength(W, N);
    I:=0;
    while I<=N-1 do
    begin
        W[I] := 1;
        Inc(I);
    end;
    LSFitLinearInternal(Y, W, FMatrix, N, M, Info, C, Rep);
end;


(*************************************************************************
Constained linear least squares fitting, without weights.

See LSFitLinearWC() for more information.

  -- ALGLIB --
     Copyright 07.09.2009 by Bochkanov Sergey
*************************************************************************)
procedure LSFitLinearC(Y : TReal1DArray;
     const FMatrix : TReal2DArray;
     const CMatrix : TReal2DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     K : AlglibInteger;
     var Info : AlglibInteger;
     var C : TReal1DArray;
     var Rep : LSFitReport);
var
    W : TReal1DArray;
    I : AlglibInteger;
begin
    Y := DynamicArrayCopy(Y);
    if N<1 then
    begin
        Info := -1;
        Exit;
    end;
    SetLength(W, N);
    I:=0;
    while I<=N-1 do
    begin
        W[I] := 1;
        Inc(I);
    end;
    LSFitLinearWC(Y, W, FMatrix, CMatrix, N, M, K, Info, C, Rep);
end;


(*************************************************************************
Weighted nonlinear least squares fitting using gradient and Hessian.

Nonlinear task min(F(c)) is solved, where

    F(c) = (w[0]*(f(x[0],c)-y[0]))^2 + ... + (w[n-1]*(f(x[n-1],c)-y[n-1]))^2,
    
    * N is a number of points,
    * M is a dimension of a space points belong to,
    * K is a dimension of a space of parameters being fitted,
    * w is an N-dimensional vector of weight coefficients,
    * x is a set of N points, each of them is an M-dimensional vector,
    * c is a K-dimensional vector of parameters being fitted
    
This subroutine uses only f(x[i],c) and its gradient.
    
INPUT PARAMETERS:
    X       -   array[0..N-1,0..M-1], points (one row = one point)
    Y       -   array[0..N-1], function values.
    W       -   weights, array[0..N-1]
    C       -   array[0..K-1], initial approximation to the solution,
    N       -   number of points, N>1
    M       -   dimension of space
    K       -   number of parameters being fitted
    CheapFG -   boolean flag, which is:
                * True  if both function and gradient calculation complexity
                        are less than O(M^2).  An improved  algorithm  can
                        be  used  which corresponds  to  FGJ  scheme  from
                        MINLM unit.
                * False otherwise.
                        Standard Jacibian-bases  Levenberg-Marquardt  algo
                        will be used (FJ scheme).

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state between subsequent
                calls  of   LSFitNonlinearIteration.   Used  for  reverse
                communication.  This  structure   should   be  passed  to
                LSFitNonlinearIteration subroutine.

See also:
    LSFitNonlinearIteration
    LSFitNonlinearResults
    LSFitNonlinearFG (fitting without weights)
    LSFitNonlinearWFGH (fitting using Hessian)
    LSFitNonlinearFGH (fitting using Hessian, without weights)


  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure LSFitNonlinearWFG(const X : TReal2DArray;
     const Y : TReal1DArray;
     const W : TReal1DArray;
     const C : TReal1DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     K : AlglibInteger;
     CheapFG : Boolean;
     var State : LSFitState);
var
    I : AlglibInteger;
begin
    State.N := N;
    State.M := M;
    State.K := K;
    LSFitNonLinearSetCond(State, Double(0.0), Double(0.0), 0);
    LSFitNonLinearSetStpMax(State, Double(0.0));
    State.CheapFG := CheapFG;
    State.HaveHess := False;
    if (N>=1) and (M>=1) and (K>=1) then
    begin
        SetLength(State.TaskX, N, M);
        SetLength(State.TaskY, N);
        SetLength(State.W, N);
        SetLength(State.C, K);
        APVMove(@State.C[0], 0, K-1, @C[0], 0, K-1);
        APVMove(@State.W[0], 0, N-1, @W[0], 0, N-1);
        I:=0;
        while I<=N-1 do
        begin
            APVMove(@State.TaskX[I][0], 0, M-1, @X[I][0], 0, M-1);
            State.TaskY[I] := Y[I];
            Inc(I);
        end;
    end;
    SetLength(State.RState.IA, 4+1);
    SetLength(State.RState.RA, 1+1);
    State.RState.Stage := -1;
end;


(*************************************************************************
Nonlinear least squares fitting, no individual weights.
See LSFitNonlinearWFG for more information.

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure LSFitNonlinearFG(const X : TReal2DArray;
     const Y : TReal1DArray;
     const C : TReal1DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     K : AlglibInteger;
     CheapFG : Boolean;
     var State : LSFitState);
var
    I : AlglibInteger;
begin
    State.N := N;
    State.M := M;
    State.K := K;
    LSFitNonLinearSetCond(State, Double(0.0), Double(0.0), 0);
    LSFitNonLinearSetStpMax(State, Double(0.0));
    State.CheapFG := CheapFG;
    State.HaveHess := False;
    if (N>=1) and (M>=1) and (K>=1) then
    begin
        SetLength(State.TaskX, N, M);
        SetLength(State.TaskY, N);
        SetLength(State.W, N);
        SetLength(State.C, K);
        APVMove(@State.C[0], 0, K-1, @C[0], 0, K-1);
        I:=0;
        while I<=N-1 do
        begin
            APVMove(@State.TaskX[I][0], 0, M-1, @X[I][0], 0, M-1);
            State.TaskY[I] := Y[I];
            State.W[I] := 1;
            Inc(I);
        end;
    end;
    SetLength(State.RState.IA, 4+1);
    SetLength(State.RState.RA, 1+1);
    State.RState.Stage := -1;
end;


(*************************************************************************
Weighted nonlinear least squares fitting using gradient/Hessian.

Nonlinear task min(F(c)) is solved, where

    F(c) = (w[0]*(f(x[0],c)-y[0]))^2 + ... + (w[n-1]*(f(x[n-1],c)-y[n-1]))^2,

    * N is a number of points,
    * M is a dimension of a space points belong to,
    * K is a dimension of a space of parameters being fitted,
    * w is an N-dimensional vector of weight coefficients,
    * x is a set of N points, each of them is an M-dimensional vector,
    * c is a K-dimensional vector of parameters being fitted

This subroutine uses f(x[i],c), its gradient and its Hessian.

See LSFitNonlinearWFG() subroutine for information about function
parameters.

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure LSFitNonlinearWFGH(const X : TReal2DArray;
     const Y : TReal1DArray;
     const W : TReal1DArray;
     const C : TReal1DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     K : AlglibInteger;
     var State : LSFitState);
var
    I : AlglibInteger;
begin
    State.N := N;
    State.M := M;
    State.K := K;
    LSFitNonLinearSetCond(State, Double(0.0), Double(0.0), 0);
    LSFitNonLinearSetStpMax(State, Double(0.0));
    State.CheapFG := True;
    State.HaveHess := True;
    if (N>=1) and (M>=1) and (K>=1) then
    begin
        SetLength(State.TaskX, N, M);
        SetLength(State.TaskY, N);
        SetLength(State.W, N);
        SetLength(State.C, K);
        APVMove(@State.C[0], 0, K-1, @C[0], 0, K-1);
        APVMove(@State.W[0], 0, N-1, @W[0], 0, N-1);
        I:=0;
        while I<=N-1 do
        begin
            APVMove(@State.TaskX[I][0], 0, M-1, @X[I][0], 0, M-1);
            State.TaskY[I] := Y[I];
            Inc(I);
        end;
    end;
    SetLength(State.RState.IA, 4+1);
    SetLength(State.RState.RA, 1+1);
    State.RState.Stage := -1;
end;


(*************************************************************************
Nonlinear least squares fitting using gradient/Hessian without  individual
weights. See LSFitNonlinearWFGH() for more information.


  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure LSFitNonlinearFGH(const X : TReal2DArray;
     const Y : TReal1DArray;
     const C : TReal1DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     K : AlglibInteger;
     var State : LSFitState);
var
    I : AlglibInteger;
begin
    State.N := N;
    State.M := M;
    State.K := K;
    LSFitNonLinearSetCond(State, Double(0.0), Double(0.0), 0);
    LSFitNonLinearSetStpMax(State, Double(0.0));
    State.CheapFG := True;
    State.HaveHess := True;
    if (N>=1) and (M>=1) and (K>=1) then
    begin
        SetLength(State.TaskX, N, M);
        SetLength(State.TaskY, N);
        SetLength(State.W, N);
        SetLength(State.C, K);
        APVMove(@State.C[0], 0, K-1, @C[0], 0, K-1);
        I:=0;
        while I<=N-1 do
        begin
            APVMove(@State.TaskX[I][0], 0, M-1, @X[I][0], 0, M-1);
            State.TaskY[I] := Y[I];
            State.W[I] := 1;
            Inc(I);
        end;
    end;
    SetLength(State.RState.IA, 4+1);
    SetLength(State.RState.RA, 1+1);
    State.RState.Stage := -1;
end;


(*************************************************************************
Stopping conditions for nonlinear least squares fitting.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be initialized
                with LSFitNonLinearCreate???()
    EpsF    -   stopping criterion. Algorithm stops if
                |F(k+1)-F(k)| <= EpsF*max{|F(k)|, |F(k+1)|, 1}
    EpsX    -   stopping criterion. Algorithm stops if
                |X(k+1)-X(k)| <= EpsX*(1+|X(k)|)
    MaxIts  -   stopping criterion. Algorithm stops after MaxIts iterations.
                MaxIts=0 means no stopping criterion.

NOTE

Passing EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to automatic
stopping criterion selection (according to the scheme used by MINLM unit).


  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure LSFitNonlinearSetCond(var State : LSFitState;
     EpsF : Double;
     EpsX : Double;
     MaxIts : AlglibInteger);
begin
    Assert(AP_FP_Greater_Eq(EpsF,0), 'LSFitNonlinearSetCond: negative EpsF!');
    Assert(AP_FP_Greater_Eq(EpsX,0), 'LSFitNonlinearSetCond: negative EpsX!');
    Assert(MaxIts>=0, 'LSFitNonlinearSetCond: negative MaxIts!');
    State.EpsF := EpsF;
    State.EpsX := EpsX;
    State.MaxIts := MaxIts;
end;


(*************************************************************************
This function sets maximum step length

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be
                initialized with LSFitNonLinearCreate???()
    StpMax  -   maximum step length, >=0. Set StpMax to 0.0,  if you don't
                want to limit step length.

Use this subroutine when you optimize target function which contains exp()
or  other  fast  growing  functions,  and optimization algorithm makes too
large  steps  which  leads  to overflow. This function allows us to reject
steps  that  are  too  large  (and  therefore  expose  us  to the possible
overflow) without actually calculating function value at the x+stp*d.

NOTE: non-zero StpMax leads to moderate  performance  degradation  because
intermediate  step  of  preconditioned L-BFGS optimization is incompatible
with limits on step size.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************)
procedure LSFitNonlinearSetStpMax(var State : LSFitState; StpMax : Double);
begin
    Assert(AP_FP_Greater_Eq(StpMax,0), 'LSFitNonlinearSetStpMax: StpMax<0!');
    State.StpMax := StpMax;
end;


(*************************************************************************
Nonlinear least squares fitting. Algorithm iteration.

Called after inialization of the State structure with  LSFitNonlinearXXX()
subroutine. See HTML docs for examples.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state between  subsequent
                calls and which is used for reverse communication. Must be
                initialized with LSFitNonlinearXXX() call first.

RESULT
1. If subroutine returned False, iterative algorithm has converged.
2. If subroutine returned True, then if:
* if State.NeedF=True,      function value F(X,C) is required
* if State.NeedFG=True,     function value F(X,C) and gradient  dF/dC(X,C)
                            are required
* if State.NeedFGH=True     function value F(X,C), gradient dF/dC(X,C) and
                            Hessian are required

One and only one of this fields can be set at time.

Function, its gradient and Hessian are calculated at  (X,C),  where  X  is
stored in State.X[0..M-1] and C is stored in State.C[0..K-1].

Results are stored:
* function value            -   in State.F
* gradient                  -   in State.G[0..K-1]
* Hessian                   -   in State.H[0..K-1,0..K-1]

  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************)
function LSFitNonlinearIteration(var State : LSFitState):Boolean;
var
    N : AlglibInteger;
    M : AlglibInteger;
    K : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    RelCnt : Double;
label
lbl_5, lbl_9, lbl_0, lbl_11, lbl_7, lbl_14, lbl_1, lbl_16, lbl_12, lbl_19, lbl_2, lbl_21, lbl_17, lbl_24, lbl_3, lbl_26, lbl_22, lbl_6, lbl_29, lbl_4, lbl_31, lbl_27, lbl_rcomm;
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
        K := State.RState.IA[2];
        I := State.RState.IA[3];
        J := State.RState.IA[4];
        V := State.RState.RA[0];
        RelCnt := State.RState.RA[1];
    end
    else
    begin
        N := -983;
        M := -989;
        K := -834;
        I := 900;
        J := -287;
        V := 364;
        RelCnt := 214;
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
    
    //
    // Routine body
    //
    
    //
    // check params
    //
    if (State.N<1) or (State.M<1) or (State.K<1) or AP_FP_Less(State.EpsF,0) or AP_FP_Less(State.EpsX,0) or (State.MaxIts<0) then
    begin
        State.RepTerminationType := -1;
        Result := False;
        Exit;
    end;
    
    //
    // init
    //
    N := State.N;
    M := State.M;
    K := State.K;
    SetLength(State.X, M);
    SetLength(State.G, K);
    if State.HaveHess then
    begin
        SetLength(State.H, K, K);
    end;
    
    //
    // initialize LM optimizer
    //
    if State.HaveHess then
    begin
        
        //
        // use Hessian.
        // transform stopping conditions.
        //
        MinLMCreateFGH(K, State.C, State.OptState);
    end
    else
    begin
        
        //
        // use one of gradient-based schemes (depending on gradient cost).
        // transform stopping conditions.
        //
        if State.CheapFG then
        begin
            MinLMCreateFGJ(K, N, State.C, State.OptState);
        end
        else
        begin
            MinLMCreateFJ(K, N, State.C, State.OptState);
        end;
    end;
    MinLMSetCond(State.OptState, Double(0.0), State.EpsF, State.EpsX, State.MaxIts);
    MinLMSetStpMax(State.OptState, State.StpMax);
    
    //
    // Optimize
    //
lbl_5:
    if not MinLMIteration(State.OptState) then
    begin
        goto lbl_6;
    end;
    if not State.OptState.NeedF then
    begin
        goto lbl_7;
    end;
    
    //
    // calculate F = sum (wi*(f(xi,c)-yi))^2
    //
    State.OptState.F := 0;
    I := 0;
lbl_9:
    if I>N-1 then
    begin
        goto lbl_11;
    end;
    APVMove(@State.C[0], 0, K-1, @State.OptState.X[0], 0, K-1);
    APVMove(@State.X[0], 0, M-1, @State.TaskX[I][0], 0, M-1);
    State.PointIndex := I;
    LSFitClearRequestFields(State);
    State.NeedF := True;
    State.RState.Stage := 0;
    goto lbl_rcomm;
lbl_0:
    State.OptState.F := State.OptState.F+AP_Sqr(State.W[I]*(State.F-State.TaskY[I]));
    I := I+1;
    goto lbl_9;
lbl_11:
    goto lbl_5;
lbl_7:
    if not State.OptState.NeedFG then
    begin
        goto lbl_12;
    end;
    
    //
    // calculate F/gradF
    //
    State.OptState.F := 0;
    I:=0;
    while I<=K-1 do
    begin
        State.OptState.G[I] := 0;
        Inc(I);
    end;
    I := 0;
lbl_14:
    if I>N-1 then
    begin
        goto lbl_16;
    end;
    APVMove(@State.C[0], 0, K-1, @State.OptState.X[0], 0, K-1);
    APVMove(@State.X[0], 0, M-1, @State.TaskX[I][0], 0, M-1);
    State.PointIndex := I;
    LSFitClearRequestFields(State);
    State.NeedFG := True;
    State.RState.Stage := 1;
    goto lbl_rcomm;
lbl_1:
    State.OptState.F := State.OptState.F+AP_Sqr(State.W[I]*(State.F-State.TaskY[I]));
    V := AP_Sqr(State.W[I])*2*(State.F-State.TaskY[I]);
    APVAdd(@State.OptState.G[0], 0, K-1, @State.G[0], 0, K-1, V);
    I := I+1;
    goto lbl_14;
lbl_16:
    goto lbl_5;
lbl_12:
    if not State.OptState.NeedFiJ then
    begin
        goto lbl_17;
    end;
    
    //
    // calculate Fi/jac(Fi)
    //
    I := 0;
lbl_19:
    if I>N-1 then
    begin
        goto lbl_21;
    end;
    APVMove(@State.C[0], 0, K-1, @State.OptState.X[0], 0, K-1);
    APVMove(@State.X[0], 0, M-1, @State.TaskX[I][0], 0, M-1);
    State.PointIndex := I;
    LSFitClearRequestFields(State);
    State.NeedFG := True;
    State.RState.Stage := 2;
    goto lbl_rcomm;
lbl_2:
    State.OptState.FI[I] := State.W[I]*(State.F-State.TaskY[I]);
    V := State.W[I];
    APVMove(@State.OptState.J[I][0], 0, K-1, @State.G[0], 0, K-1, V);
    I := I+1;
    goto lbl_19;
lbl_21:
    goto lbl_5;
lbl_17:
    if not State.OptState.NeedFGH then
    begin
        goto lbl_22;
    end;
    
    //
    // calculate F/grad(F)/hess(F)
    //
    State.OptState.F := 0;
    I:=0;
    while I<=K-1 do
    begin
        State.OptState.G[I] := 0;
        Inc(I);
    end;
    I:=0;
    while I<=K-1 do
    begin
        J:=0;
        while J<=K-1 do
        begin
            State.OptState.H[I,J] := 0;
            Inc(J);
        end;
        Inc(I);
    end;
    I := 0;
lbl_24:
    if I>N-1 then
    begin
        goto lbl_26;
    end;
    APVMove(@State.C[0], 0, K-1, @State.OptState.X[0], 0, K-1);
    APVMove(@State.X[0], 0, M-1, @State.TaskX[I][0], 0, M-1);
    State.PointIndex := I;
    LSFitClearRequestFields(State);
    State.NeedFGH := True;
    State.RState.Stage := 3;
    goto lbl_rcomm;
lbl_3:
    State.OptState.F := State.OptState.F+AP_Sqr(State.W[I]*(State.F-State.TaskY[I]));
    V := AP_Sqr(State.W[I])*2*(State.F-State.TaskY[I]);
    APVAdd(@State.OptState.G[0], 0, K-1, @State.G[0], 0, K-1, V);
    J:=0;
    while J<=K-1 do
    begin
        V := 2*AP_Sqr(State.W[I])*State.G[J];
        APVAdd(@State.OptState.H[J][0], 0, K-1, @State.G[0], 0, K-1, V);
        V := 2*AP_Sqr(State.W[I])*(State.F-State.TaskY[I]);
        APVAdd(@State.OptState.H[J][0], 0, K-1, @State.H[J][0], 0, K-1, V);
        Inc(J);
    end;
    I := I+1;
    goto lbl_24;
lbl_26:
    goto lbl_5;
lbl_22:
    goto lbl_5;
lbl_6:
    MinLMResults(State.OptState, State.C, State.OptRep);
    State.RepTerminationType := State.OptRep.TerminationType;
    
    //
    // calculate errors
    //
    if State.RepTerminationType<=0 then
    begin
        goto lbl_27;
    end;
    State.RepRMSError := 0;
    State.RepAvgError := 0;
    State.RepAvgRelError := 0;
    State.RepMaxError := 0;
    RelCnt := 0;
    I := 0;
lbl_29:
    if I>N-1 then
    begin
        goto lbl_31;
    end;
    APVMove(@State.C[0], 0, K-1, @State.C[0], 0, K-1);
    APVMove(@State.X[0], 0, M-1, @State.TaskX[I][0], 0, M-1);
    State.PointIndex := I;
    LSFitClearRequestFields(State);
    State.NeedF := True;
    State.RState.Stage := 4;
    goto lbl_rcomm;
lbl_4:
    V := State.F;
    State.RepRMSError := State.RepRMSError+AP_Sqr(V-State.TaskY[I]);
    State.RepAvgError := State.RepAvgError+AbsReal(V-State.TaskY[I]);
    if AP_FP_Neq(State.TaskY[I],0) then
    begin
        State.RepAvgRelError := State.RepAvgRelError+AbsReal(V-State.TaskY[I])/AbsReal(State.TaskY[I]);
        RelCnt := RelCnt+1;
    end;
    State.RepMaxError := Max(State.RepMaxError, AbsReal(V-State.TaskY[I]));
    I := I+1;
    goto lbl_29;
lbl_31:
    State.RepRMSError := Sqrt(State.RepRMSError/N);
    State.RepAvgError := State.RepAvgError/N;
    if AP_FP_Neq(RelCnt,0) then
    begin
        State.RepAvgRelError := State.RepAvgRelError/RelCnt;
    end;
lbl_27:
    Result := False;
    Exit;
    
    //
    // Saving state
    //
lbl_rcomm:
    Result := True;
    State.RState.IA[0] := N;
    State.RState.IA[1] := M;
    State.RState.IA[2] := K;
    State.RState.IA[3] := I;
    State.RState.IA[4] := J;
    State.RState.RA[0] := V;
    State.RState.RA[1] := RelCnt;
end;


(*************************************************************************
Nonlinear least squares fitting results.

Called after LSFitNonlinearIteration() returned False.

INPUT PARAMETERS:
    State   -   algorithm state (used by LSFitNonlinearIteration).

OUTPUT PARAMETERS:
    Info    -   completetion code:
                    * -1    incorrect parameters were specified
                    *  1    relative function improvement is no more than
                            EpsF.
                    *  2    relative step is no more than EpsX.
                    *  4    gradient norm is no more than EpsG
                    *  5    MaxIts steps was taken
    C       -   array[0..K-1], solution
    Rep     -   optimization report. Following fields are set:
                * Rep.TerminationType completetion code:
                * RMSError          rms error on the (X,Y).
                * AvgError          average error on the (X,Y).
                * AvgRelError       average relative error on the non-zero Y
                * MaxError          maximum error
                                    NON-WEIGHTED ERRORS ARE CALCULATED


  -- ALGLIB --
     Copyright 17.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure LSFitNonlinearResults(const State : LSFitState;
     var Info : AlglibInteger;
     var C : TReal1DArray;
     var Rep : LSFitReport);
begin
    Info := State.RepTerminationType;
    if Info>0 then
    begin
        SetLength(C, State.K);
        APVMove(@C[0], 0, State.K-1, @State.C[0], 0, State.K-1);
        Rep.RMSError := State.RepRMSError;
        Rep.AvgError := State.RepAvgError;
        Rep.AvgRelError := State.RepAvgRelError;
        Rep.MaxError := State.RepMaxError;
    end;
end;


procedure LSFitScaleXY(var X : TReal1DArray;
     var Y : TReal1DArray;
     N : AlglibInteger;
     var XC : TReal1DArray;
     var YC : TReal1DArray;
     const DC : TInteger1DArray;
     K : AlglibInteger;
     var XA : Double;
     var XB : Double;
     var SA : Double;
     var SB : Double;
     var XOriginal : TReal1DArray;
     var YOriginal : TReal1DArray);
var
    XMin : Double;
    XMax : Double;
    I : AlglibInteger;
begin
    Assert(N>=1, 'LSFitScaleXY: incorrect N');
    Assert(K>=0, 'LSFitScaleXY: incorrect K');
    
    //
    // Calculate xmin/xmax.
    // Force xmin<>xmax.
    //
    XMin := X[0];
    XMax := X[0];
    I:=1;
    while I<=N-1 do
    begin
        XMin := Min(XMin, X[I]);
        XMax := Max(XMax, X[I]);
        Inc(I);
    end;
    I:=0;
    while I<=K-1 do
    begin
        XMin := Min(XMin, XC[I]);
        XMax := Max(XMax, XC[I]);
        Inc(I);
    end;
    if AP_FP_Eq(XMin,XMax) then
    begin
        if AP_FP_Eq(XMin,0) then
        begin
            XMin := -1;
            XMax := +1;
        end
        else
        begin
            XMin := Double(0.5)*XMin;
        end;
    end;
    
    //
    // Transform abscissas: map [XA,XB] to [0,1]
    //
    // Store old X[] in XOriginal[] (it will be used
    // to calculate relative error).
    //
    SetLength(XOriginal, N);
    APVMove(@XOriginal[0], 0, N-1, @X[0], 0, N-1);
    XA := XMin;
    XB := XMax;
    I:=0;
    while I<=N-1 do
    begin
        X[I] := 2*(X[I]-Double(0.5)*(XA+XB))/(XB-XA);
        Inc(I);
    end;
    I:=0;
    while I<=K-1 do
    begin
        Assert(DC[I]>=0, 'LSFitScaleXY: internal error!');
        XC[I] := 2*(XC[I]-Double(0.5)*(XA+XB))/(XB-XA);
        YC[I] := YC[I]*Power(Double(0.5)*(XB-XA), DC[I]);
        Inc(I);
    end;
    
    //
    // Transform function values: map [SA,SB] to [0,1]
    // SA = mean(Y),
    // SB = SA+stddev(Y).
    //
    // Store old Y[] in YOriginal[] (it will be used
    // to calculate relative error).
    //
    SetLength(YOriginal, N);
    APVMove(@YOriginal[0], 0, N-1, @Y[0], 0, N-1);
    SA := 0;
    I:=0;
    while I<=N-1 do
    begin
        SA := SA+Y[I];
        Inc(I);
    end;
    SA := SA/N;
    SB := 0;
    I:=0;
    while I<=N-1 do
    begin
        SB := SB+AP_Sqr(Y[I]-SA);
        Inc(I);
    end;
    SB := Sqrt(SB/N)+SA;
    if AP_FP_Eq(SB,SA) then
    begin
        SB := 2*SA;
    end;
    if AP_FP_Eq(SB,SA) then
    begin
        SB := SA+1;
    end;
    I:=0;
    while I<=N-1 do
    begin
        Y[I] := (Y[I]-SA)/(SB-SA);
        Inc(I);
    end;
    I:=0;
    while I<=K-1 do
    begin
        if DC[I]=0 then
        begin
            YC[I] := (YC[I]-SA)/(SB-SA);
        end
        else
        begin
            YC[I] := YC[I]/(SB-SA);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Internal fitting subroutine
*************************************************************************)
procedure LSFitLinearInternal(const Y : TReal1DArray;
     const W : TReal1DArray;
     const FMatrix : TReal2DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var C : TReal1DArray;
     var Rep : LSFitReport);
var
    Threshold : Double;
    FT : TReal2DArray;
    Q : TReal2DArray;
    L : TReal2DArray;
    R : TReal2DArray;
    B : TReal1DArray;
    WMod : TReal1DArray;
    Tau : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    SV : TReal1DArray;
    U : TReal2DArray;
    VT : TReal2DArray;
    Tmp : TReal1DArray;
    UTB : TReal1DArray;
    SUTB : TReal1DArray;
    RelCnt : AlglibInteger;
begin
    if (N<1) or (M<1) then
    begin
        Info := -1;
        Exit;
    end;
    Info := 1;
    Threshold := Sqrt(MachineEpsilon);
    
    //
    // Degenerate case, needs special handling
    //
    if N<M then
    begin
        
        //
        // Create design matrix.
        //
        SetLength(FT, N, M);
        SetLength(B, N);
        SetLength(WMod, N);
        J:=0;
        while J<=N-1 do
        begin
            V := W[J];
            APVMove(@FT[J][0], 0, M-1, @FMatrix[J][0], 0, M-1, V);
            B[J] := W[J]*Y[J];
            WMod[J] := 1;
            Inc(J);
        end;
        
        //
        // LQ decomposition and reduction to M=N
        //
        SetLength(C, M);
        I:=0;
        while I<=M-1 do
        begin
            C[I] := 0;
            Inc(I);
        end;
        Rep.TaskRCond := 0;
        RMatrixLQ(FT, N, M, Tau);
        RMatrixLQUnpackQ(FT, N, M, Tau, N, Q);
        RMatrixLQUnpackL(FT, N, M, L);
        LSFitLinearInternal(B, WMod, L, N, N, Info, Tmp, Rep);
        if Info<=0 then
        begin
            Exit;
        end;
        I:=0;
        while I<=N-1 do
        begin
            V := Tmp[I];
            APVAdd(@C[0], 0, M-1, @Q[I][0], 0, M-1, V);
            Inc(I);
        end;
        Exit;
    end;
    
    //
    // N>=M. Generate design matrix and reduce to N=M using
    // QR decomposition.
    //
    SetLength(FT, N, M);
    SetLength(B, N);
    J:=0;
    while J<=N-1 do
    begin
        V := W[J];
        APVMove(@FT[J][0], 0, M-1, @FMatrix[J][0], 0, M-1, V);
        B[J] := W[J]*Y[J];
        Inc(J);
    end;
    RMatrixQR(FT, N, M, Tau);
    RMatrixQRUnpackQ(FT, N, M, Tau, M, Q);
    RMatrixQRUnpackR(FT, N, M, R);
    SetLength(Tmp, M);
    I:=0;
    while I<=M-1 do
    begin
        Tmp[I] := 0;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        V := B[I];
        APVAdd(@Tmp[0], 0, M-1, @Q[I][0], 0, M-1, V);
        Inc(I);
    end;
    SetLength(B, M);
    APVMove(@B[0], 0, M-1, @Tmp[0], 0, M-1);
    
    //
    // R contains reduced MxM design upper triangular matrix,
    // B contains reduced Mx1 right part.
    //
    // Determine system condition number and decide
    // should we use triangular solver (faster) or
    // SVD-based solver (more stable).
    //
    // We can use LU-based RCond estimator for this task.
    //
    Rep.TaskRCond := RMatrixLURCondInf(R, M);
    if AP_FP_Greater(Rep.TaskRCond,Threshold) then
    begin
        
        //
        // use QR-based solver
        //
        SetLength(C, M);
        C[M-1] := B[M-1]/R[M-1,M-1];
        I:=M-2;
        while I>=0 do
        begin
            V := APVDotProduct(@R[I][0], I+1, M-1, @C[0], I+1, M-1);
            C[I] := (B[I]-V)/R[I,I];
            Dec(I);
        end;
    end
    else
    begin
        
        //
        // use SVD-based solver
        //
        if  not RMatrixSVD(R, M, M, 1, 1, 2, SV, U, VT) then
        begin
            Info := -4;
            Exit;
        end;
        SetLength(UTB, M);
        SetLength(SUTB, M);
        I:=0;
        while I<=M-1 do
        begin
            UTB[I] := 0;
            Inc(I);
        end;
        I:=0;
        while I<=M-1 do
        begin
            V := B[I];
            APVAdd(@UTB[0], 0, M-1, @U[I][0], 0, M-1, V);
            Inc(I);
        end;
        if AP_FP_Greater(SV[0],0) then
        begin
            Rep.TaskRCond := SV[M-1]/SV[0];
            I:=0;
            while I<=M-1 do
            begin
                if AP_FP_Greater(SV[I],Threshold*SV[0]) then
                begin
                    SUTB[I] := UTB[I]/SV[I];
                end
                else
                begin
                    SUTB[I] := 0;
                end;
                Inc(I);
            end;
        end
        else
        begin
            Rep.TaskRCond := 0;
            I:=0;
            while I<=M-1 do
            begin
                SUTB[I] := 0;
                Inc(I);
            end;
        end;
        SetLength(C, M);
        I:=0;
        while I<=M-1 do
        begin
            C[I] := 0;
            Inc(I);
        end;
        I:=0;
        while I<=M-1 do
        begin
            V := SUTB[I];
            APVAdd(@C[0], 0, M-1, @VT[I][0], 0, M-1, V);
            Inc(I);
        end;
    end;
    
    //
    // calculate errors
    //
    Rep.RMSError := 0;
    Rep.AvgError := 0;
    Rep.AvgRelError := 0;
    Rep.MaxError := 0;
    RelCnt := 0;
    I:=0;
    while I<=N-1 do
    begin
        V := APVDotProduct(@FMatrix[I][0], 0, M-1, @C[0], 0, M-1);
        Rep.RMSError := Rep.RMSError+AP_Sqr(V-Y[I]);
        Rep.AvgError := Rep.AvgError+AbsReal(V-Y[I]);
        if AP_FP_Neq(Y[I],0) then
        begin
            Rep.AvgRelError := Rep.AvgRelError+AbsReal(V-Y[I])/AbsReal(Y[I]);
            RelCnt := RelCnt+1;
        end;
        Rep.MaxError := Max(Rep.MaxError, AbsReal(V-Y[I]));
        Inc(I);
    end;
    Rep.RMSError := Sqrt(Rep.RMSError/N);
    Rep.AvgError := Rep.AvgError/N;
    if RelCnt<>0 then
    begin
        Rep.AvgRelError := Rep.AvgRelError/RelCnt;
    end;
end;


(*************************************************************************
Internal subroutine
*************************************************************************)
procedure LSFitClearRequestFields(var State : LSFitState);
begin
    State.NeedF := False;
    State.NeedFG := False;
    State.NeedFGH := False;
end;


end.