{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Copyright (c) 2007-2008, Sergey Bochkanov (ALGLIB project).

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
unit densesolver;
interface
uses Math, Sysutils, Ap, hblas, reflections, creflections, sblas, ablasf, ablas, ortfac, blas, rotations, bdsvd, svd, hqrnd, matgen, trfac, trlinsolve, safesolve, rcond, xblas;

type
DenseSolverReport = record
    R1 : Double;
    RInf : Double;
end;


DenseSolverLSReport = record
    R2 : Double;
    CX : TReal2DArray;
    N : AlglibInteger;
    K : AlglibInteger;
end;



procedure RMatrixSolve(const A : TReal2DArray;
     N : AlglibInteger;
     const B : TReal1DArray;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TReal1DArray);
procedure RMatrixSolveM(const A : TReal2DArray;
     N : AlglibInteger;
     const B : TReal2DArray;
     M : AlglibInteger;
     RFS : Boolean;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TReal2DArray);
procedure RMatrixLUSolve(const LUA : TReal2DArray;
     const P : TInteger1DArray;
     N : AlglibInteger;
     const B : TReal1DArray;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TReal1DArray);
procedure RMatrixLUSolveM(const LUA : TReal2DArray;
     const P : TInteger1DArray;
     N : AlglibInteger;
     const B : TReal2DArray;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TReal2DArray);
procedure RMatrixMixedSolve(const A : TReal2DArray;
     const LUA : TReal2DArray;
     const P : TInteger1DArray;
     N : AlglibInteger;
     const B : TReal1DArray;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TReal1DArray);
procedure RMatrixMixedSolveM(const A : TReal2DArray;
     const LUA : TReal2DArray;
     const P : TInteger1DArray;
     N : AlglibInteger;
     const B : TReal2DArray;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TReal2DArray);
procedure CMatrixSolveM(const A : TComplex2DArray;
     N : AlglibInteger;
     const B : TComplex2DArray;
     M : AlglibInteger;
     RFS : Boolean;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TComplex2DArray);
procedure CMatrixSolve(const A : TComplex2DArray;
     N : AlglibInteger;
     const B : TComplex1DArray;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TComplex1DArray);
procedure CMatrixLUSolveM(const LUA : TComplex2DArray;
     const P : TInteger1DArray;
     N : AlglibInteger;
     const B : TComplex2DArray;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TComplex2DArray);
procedure CMatrixLUSolve(const LUA : TComplex2DArray;
     const P : TInteger1DArray;
     N : AlglibInteger;
     const B : TComplex1DArray;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TComplex1DArray);
procedure CMatrixMixedSolveM(const A : TComplex2DArray;
     const LUA : TComplex2DArray;
     const P : TInteger1DArray;
     N : AlglibInteger;
     const B : TComplex2DArray;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TComplex2DArray);
procedure CMatrixMixedSolve(const A : TComplex2DArray;
     const LUA : TComplex2DArray;
     const P : TInteger1DArray;
     N : AlglibInteger;
     const B : TComplex1DArray;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TComplex1DArray);
procedure SPDMatrixSolveM(const A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     const B : TReal2DArray;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TReal2DArray);
procedure SPDMatrixSolve(const A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     const B : TReal1DArray;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TReal1DArray);
procedure SPDMatrixCholeskySolveM(const CHA : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     const B : TReal2DArray;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TReal2DArray);
procedure SPDMatrixCholeskySolve(const CHA : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     const B : TReal1DArray;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TReal1DArray);
procedure HPDMatrixSolveM(const A : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     const B : TComplex2DArray;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TComplex2DArray);
procedure HPDMatrixSolve(const A : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     const B : TComplex1DArray;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TComplex1DArray);
procedure HPDMatrixCholeskySolveM(const CHA : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     const B : TComplex2DArray;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TComplex2DArray);
procedure HPDMatrixCholeskySolve(const CHA : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     const B : TComplex1DArray;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TComplex1DArray);
procedure RMatrixSolveLS(const A : TReal2DArray;
     NRows : AlglibInteger;
     NCols : AlglibInteger;
     const B : TReal1DArray;
     Threshold : Double;
     var Info : AlglibInteger;
     var Rep : DenseSolverLSReport;
     var X : TReal1DArray);

implementation

procedure RMatrixLUSolveInternal(const LUA : TReal2DArray;
     const P : TInteger1DArray;
     const ScaleA : Double;
     N : AlglibInteger;
     const A : TReal2DArray;
     HaveA : Boolean;
     const B : TReal2DArray;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TReal2DArray);forward;
procedure SPDMatrixCholeskySolveInternal(const CHA : TReal2DArray;
     const SqrtScaleA : Double;
     N : AlglibInteger;
     IsUpper : Boolean;
     const A : TReal2DArray;
     HaveA : Boolean;
     const B : TReal2DArray;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TReal2DArray);forward;
procedure CMatrixLUSolveInternal(const LUA : TComplex2DArray;
     const P : TInteger1DArray;
     const ScaleA : Double;
     N : AlglibInteger;
     const A : TComplex2DArray;
     HaveA : Boolean;
     const B : TComplex2DArray;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TComplex2DArray);forward;
procedure HPDMatrixCholeskySolveInternal(const CHA : TComplex2DArray;
     const SqrtScaleA : Double;
     N : AlglibInteger;
     IsUpper : Boolean;
     const A : TComplex2DArray;
     HaveA : Boolean;
     const B : TComplex2DArray;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TComplex2DArray);forward;
function DenseSolverRFSMax(N : AlglibInteger;
     R1 : Double;
     RInf : Double):AlglibInteger;forward;
function DenseSolverRFSMaxV2(N : AlglibInteger;
     R2 : Double):AlglibInteger;forward;
procedure RBasicLUSolve(const LUA : TReal2DArray;
     const P : TInteger1DArray;
     ScaleA : Double;
     N : AlglibInteger;
     var XB : TReal1DArray;
     var Tmp : TReal1DArray);forward;
procedure SPDBasicCholeskySolve(const CHA : TReal2DArray;
     SqrtScaleA : Double;
     N : AlglibInteger;
     IsUpper : Boolean;
     var XB : TReal1DArray;
     var Tmp : TReal1DArray);forward;
procedure CBasicLUSolve(const LUA : TComplex2DArray;
     const P : TInteger1DArray;
     ScaleA : Double;
     N : AlglibInteger;
     var XB : TComplex1DArray;
     var Tmp : TComplex1DArray);forward;
procedure HPDBasicCholeskySolve(const CHA : TComplex2DArray;
     SqrtScaleA : Double;
     N : AlglibInteger;
     IsUpper : Boolean;
     var XB : TComplex1DArray;
     var Tmp : TComplex1DArray);forward;


(*************************************************************************
Dense solver.

This  subroutine  solves  a  system  A*x=b,  where A is NxN non-denegerate
real matrix, x and b are vectors.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(N^3) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   return code:
                * -3    A is singular, or VERY close to singular.
                        X is filled by zeros in such cases.
                * -1    N<=0 was passed
                *  1    task is solved (but matrix A may be ill-conditioned,
                        check R1/RInf parameters for condition numbers).
    Rep     -   solver report, see below for more info
    X       -   array[0..N-1], it contains:
                * solution of A*x=b if A is non-singular (well-conditioned
                  or ill-conditioned, but not very close to singular)
                * zeros,  if  A  is  singular  or  VERY  close to singular
                  (in this case Info=-3).

SOLVER REPORT

Subroutine sets following fields of the Rep structure:
* R1        reciprocal of condition number: 1/cond(A), 1-norm.
* RInf      reciprocal of condition number: 1/cond(A), inf-norm.

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
procedure RMatrixSolve(const A : TReal2DArray;
     N : AlglibInteger;
     const B : TReal1DArray;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TReal1DArray);
var
    BM : TReal2DArray;
    XM : TReal2DArray;
    i_ : AlglibInteger;
begin
    if N<=0 then
    begin
        Info := -1;
        Exit;
    end;
    SetLength(BM, N, 1);
    for i_ := 0 to N-1 do
    begin
        BM[i_,0] := B[i_];
    end;
    RMatrixSolveM(A, N, BM, 1, True, Info, Rep, XM);
    SetLength(X, N);
    for i_ := 0 to N-1 do
    begin
        X[i_] := XM[i_,0];
    end;
end;


(*************************************************************************
Dense solver.

Similar to RMatrixSolve() but solves task with multiple right parts (where
b and x are NxM matrices).

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* optional iterative refinement
* O(N^3+M*N^2) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size
    RFS     -   iterative refinement switch:
                * True - refinement is used.
                  Less performance, more precision.
                * False - refinement is not used.
                  More performance, less precision.

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
procedure RMatrixSolveM(const A : TReal2DArray;
     N : AlglibInteger;
     const B : TReal2DArray;
     M : AlglibInteger;
     RFS : Boolean;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TReal2DArray);
var
    DA : TReal2DArray;
    EmptyA : TReal2DArray;
    P : TInteger1DArray;
    ScaleA : Double;
    I : AlglibInteger;
    J : AlglibInteger;
begin
    
    //
    // prepare: check inputs, allocate space...
    //
    if (N<=0) or (M<=0) then
    begin
        Info := -1;
        Exit;
    end;
    SetLength(DA, N, N);
    
    //
    // 1. scale matrix, max(|A[i,j]|)
    // 2. factorize scaled matrix
    // 3. solve
    //
    ScaleA := 0;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            ScaleA := Max(ScaleA, AbsReal(A[I,J]));
            Inc(J);
        end;
        Inc(I);
    end;
    if AP_FP_Eq(ScaleA,0) then
    begin
        ScaleA := 1;
    end;
    ScaleA := 1/ScaleA;
    I:=0;
    while I<=N-1 do
    begin
        APVMove(@DA[I][0], 0, N-1, @A[I][0], 0, N-1);
        Inc(I);
    end;
    RMatrixLU(DA, N, N, P);
    if RFS then
    begin
        RMatrixLUSolveInternal(DA, P, ScaleA, N, A, True, B, M, Info, Rep, X);
    end
    else
    begin
        RMatrixLUSolveInternal(DA, P, ScaleA, N, EmptyA, False, B, M, Info, Rep, X);
    end;
end;


(*************************************************************************
Dense solver.

This  subroutine  solves  a  system  A*X=B,  where A is NxN non-denegerate
real matrix given by its LU decomposition, X and B are NxM real matrices.

Algorithm features:
* automatic detection of degenerate cases
* O(N^2) complexity
* condition number estimation

No iterative refinement  is provided because exact form of original matrix
is not known to subroutine. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
    P       -   array[0..N-1], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve
    
  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
procedure RMatrixLUSolve(const LUA : TReal2DArray;
     const P : TInteger1DArray;
     N : AlglibInteger;
     const B : TReal1DArray;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TReal1DArray);
var
    BM : TReal2DArray;
    XM : TReal2DArray;
    i_ : AlglibInteger;
begin
    if N<=0 then
    begin
        Info := -1;
        Exit;
    end;
    SetLength(BM, N, 1);
    for i_ := 0 to N-1 do
    begin
        BM[i_,0] := B[i_];
    end;
    RMatrixLUSolveM(LUA, P, N, BM, 1, Info, Rep, XM);
    SetLength(X, N);
    for i_ := 0 to N-1 do
    begin
        X[i_] := XM[i_,0];
    end;
end;


(*************************************************************************
Dense solver.

Similar to RMatrixLUSolve() but solves task with multiple right parts
(where b and x are NxM matrices).

Algorithm features:
* automatic detection of degenerate cases
* O(M*N^2) complexity
* condition number estimation

No iterative refinement  is provided because exact form of original matrix
is not known to subroutine. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
    P       -   array[0..N-1], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
procedure RMatrixLUSolveM(const LUA : TReal2DArray;
     const P : TInteger1DArray;
     N : AlglibInteger;
     const B : TReal2DArray;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TReal2DArray);
var
    EmptyA : TReal2DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    ScaleA : Double;
begin
    
    //
    // prepare: check inputs, allocate space...
    //
    if (N<=0) or (M<=0) then
    begin
        Info := -1;
        Exit;
    end;
    
    //
    // 1. scale matrix, max(|U[i,j]|)
    //    we assume that LU is in its normal form, i.e. |L[i,j]|<=1
    // 2. solve
    //
    ScaleA := 0;
    I:=0;
    while I<=N-1 do
    begin
        J:=I;
        while J<=N-1 do
        begin
            ScaleA := Max(ScaleA, AbsReal(LUA[I,J]));
            Inc(J);
        end;
        Inc(I);
    end;
    if AP_FP_Eq(ScaleA,0) then
    begin
        ScaleA := 1;
    end;
    ScaleA := 1/ScaleA;
    RMatrixLUSolveInternal(LUA, P, ScaleA, N, EmptyA, False, B, M, Info, Rep, X);
end;


(*************************************************************************
Dense solver.

This  subroutine  solves  a  system  A*x=b,  where BOTH ORIGINAL A AND ITS
LU DECOMPOSITION ARE KNOWN. You can use it if for some  reasons  you  have
both A and its LU decomposition.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(N^2) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
    P       -   array[0..N-1], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolveM
    Rep     -   same as in RMatrixSolveM
    X       -   same as in RMatrixSolveM

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
procedure RMatrixMixedSolve(const A : TReal2DArray;
     const LUA : TReal2DArray;
     const P : TInteger1DArray;
     N : AlglibInteger;
     const B : TReal1DArray;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TReal1DArray);
var
    BM : TReal2DArray;
    XM : TReal2DArray;
    i_ : AlglibInteger;
begin
    if N<=0 then
    begin
        Info := -1;
        Exit;
    end;
    SetLength(BM, N, 1);
    for i_ := 0 to N-1 do
    begin
        BM[i_,0] := B[i_];
    end;
    RMatrixMixedSolveM(A, LUA, P, N, BM, 1, Info, Rep, XM);
    SetLength(X, N);
    for i_ := 0 to N-1 do
    begin
        X[i_] := XM[i_,0];
    end;
end;


(*************************************************************************
Dense solver.

Similar to RMatrixMixedSolve() but  solves task with multiple right  parts
(where b and x are NxM matrices).

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(M*N^2) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
    P       -   array[0..N-1], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolveM
    Rep     -   same as in RMatrixSolveM
    X       -   same as in RMatrixSolveM

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
procedure RMatrixMixedSolveM(const A : TReal2DArray;
     const LUA : TReal2DArray;
     const P : TInteger1DArray;
     N : AlglibInteger;
     const B : TReal2DArray;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TReal2DArray);
var
    ScaleA : Double;
    I : AlglibInteger;
    J : AlglibInteger;
begin
    
    //
    // prepare: check inputs, allocate space...
    //
    if (N<=0) or (M<=0) then
    begin
        Info := -1;
        Exit;
    end;
    
    //
    // 1. scale matrix, max(|A[i,j]|)
    // 2. factorize scaled matrix
    // 3. solve
    //
    ScaleA := 0;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            ScaleA := Max(ScaleA, AbsReal(A[I,J]));
            Inc(J);
        end;
        Inc(I);
    end;
    if AP_FP_Eq(ScaleA,0) then
    begin
        ScaleA := 1;
    end;
    ScaleA := 1/ScaleA;
    RMatrixLUSolveInternal(LUA, P, ScaleA, N, A, True, B, M, Info, Rep, X);
end;


(*************************************************************************
Dense solver. Same as RMatrixSolveM(), but for complex matrices.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(N^3+M*N^2) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size
    RFS     -   iterative refinement switch:
                * True - refinement is used.
                  Less performance, more precision.
                * False - refinement is not used.
                  More performance, less precision.

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
procedure CMatrixSolveM(const A : TComplex2DArray;
     N : AlglibInteger;
     const B : TComplex2DArray;
     M : AlglibInteger;
     RFS : Boolean;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TComplex2DArray);
var
    DA : TComplex2DArray;
    EmptyA : TComplex2DArray;
    P : TInteger1DArray;
    ScaleA : Double;
    I : AlglibInteger;
    J : AlglibInteger;
    i_ : AlglibInteger;
begin
    
    //
    // prepare: check inputs, allocate space...
    //
    if (N<=0) or (M<=0) then
    begin
        Info := -1;
        Exit;
    end;
    SetLength(DA, N, N);
    
    //
    // 1. scale matrix, max(|A[i,j]|)
    // 2. factorize scaled matrix
    // 3. solve
    //
    ScaleA := 0;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            ScaleA := Max(ScaleA, AbsComplex(A[I,J]));
            Inc(J);
        end;
        Inc(I);
    end;
    if AP_FP_Eq(ScaleA,0) then
    begin
        ScaleA := 1;
    end;
    ScaleA := 1/ScaleA;
    I:=0;
    while I<=N-1 do
    begin
        for i_ := 0 to N-1 do
        begin
            DA[I,i_] := A[I,i_];
        end;
        Inc(I);
    end;
    CMatrixLU(DA, N, N, P);
    if RFS then
    begin
        CMatrixLUSolveInternal(DA, P, ScaleA, N, A, True, B, M, Info, Rep, X);
    end
    else
    begin
        CMatrixLUSolveInternal(DA, P, ScaleA, N, EmptyA, False, B, M, Info, Rep, X);
    end;
end;


(*************************************************************************
Dense solver. Same as RMatrixSolve(), but for complex matrices.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(N^3) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
procedure CMatrixSolve(const A : TComplex2DArray;
     N : AlglibInteger;
     const B : TComplex1DArray;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TComplex1DArray);
var
    BM : TComplex2DArray;
    XM : TComplex2DArray;
    i_ : AlglibInteger;
begin
    if N<=0 then
    begin
        Info := -1;
        Exit;
    end;
    SetLength(BM, N, 1);
    for i_ := 0 to N-1 do
    begin
        BM[i_,0] := B[i_];
    end;
    CMatrixSolveM(A, N, BM, 1, True, Info, Rep, XM);
    SetLength(X, N);
    for i_ := 0 to N-1 do
    begin
        X[i_] := XM[i_,0];
    end;
end;


(*************************************************************************
Dense solver. Same as RMatrixLUSolveM(), but for complex matrices.

Algorithm features:
* automatic detection of degenerate cases
* O(M*N^2) complexity
* condition number estimation

No iterative refinement  is provided because exact form of original matrix
is not known to subroutine. Use CMatrixSolve or CMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
    P       -   array[0..N-1], pivots array, RMatrixLU result
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
procedure CMatrixLUSolveM(const LUA : TComplex2DArray;
     const P : TInteger1DArray;
     N : AlglibInteger;
     const B : TComplex2DArray;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TComplex2DArray);
var
    EmptyA : TComplex2DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    ScaleA : Double;
begin
    
    //
    // prepare: check inputs, allocate space...
    //
    if (N<=0) or (M<=0) then
    begin
        Info := -1;
        Exit;
    end;
    
    //
    // 1. scale matrix, max(|U[i,j]|)
    //    we assume that LU is in its normal form, i.e. |L[i,j]|<=1
    // 2. solve
    //
    ScaleA := 0;
    I:=0;
    while I<=N-1 do
    begin
        J:=I;
        while J<=N-1 do
        begin
            ScaleA := Max(ScaleA, AbsComplex(LUA[I,J]));
            Inc(J);
        end;
        Inc(I);
    end;
    if AP_FP_Eq(ScaleA,0) then
    begin
        ScaleA := 1;
    end;
    ScaleA := 1/ScaleA;
    CMatrixLUSolveInternal(LUA, P, ScaleA, N, EmptyA, False, B, M, Info, Rep, X);
end;


(*************************************************************************
Dense solver. Same as RMatrixLUSolve(), but for complex matrices.

Algorithm features:
* automatic detection of degenerate cases
* O(N^2) complexity
* condition number estimation

No iterative refinement is provided because exact form of original matrix
is not known to subroutine. Use CMatrixSolve or CMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    LUA     -   array[0..N-1,0..N-1], LU decomposition, CMatrixLU result
    P       -   array[0..N-1], pivots array, CMatrixLU result
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
procedure CMatrixLUSolve(const LUA : TComplex2DArray;
     const P : TInteger1DArray;
     N : AlglibInteger;
     const B : TComplex1DArray;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TComplex1DArray);
var
    BM : TComplex2DArray;
    XM : TComplex2DArray;
    i_ : AlglibInteger;
begin
    if N<=0 then
    begin
        Info := -1;
        Exit;
    end;
    SetLength(BM, N, 1);
    for i_ := 0 to N-1 do
    begin
        BM[i_,0] := B[i_];
    end;
    CMatrixLUSolveM(LUA, P, N, BM, 1, Info, Rep, XM);
    SetLength(X, N);
    for i_ := 0 to N-1 do
    begin
        X[i_] := XM[i_,0];
    end;
end;


(*************************************************************************
Dense solver. Same as RMatrixMixedSolveM(), but for complex matrices.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(M*N^2) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    LUA     -   array[0..N-1,0..N-1], LU decomposition, CMatrixLU result
    P       -   array[0..N-1], pivots array, CMatrixLU result
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolveM
    Rep     -   same as in RMatrixSolveM
    X       -   same as in RMatrixSolveM

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
procedure CMatrixMixedSolveM(const A : TComplex2DArray;
     const LUA : TComplex2DArray;
     const P : TInteger1DArray;
     N : AlglibInteger;
     const B : TComplex2DArray;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TComplex2DArray);
var
    ScaleA : Double;
    I : AlglibInteger;
    J : AlglibInteger;
begin
    
    //
    // prepare: check inputs, allocate space...
    //
    if (N<=0) or (M<=0) then
    begin
        Info := -1;
        Exit;
    end;
    
    //
    // 1. scale matrix, max(|A[i,j]|)
    // 2. factorize scaled matrix
    // 3. solve
    //
    ScaleA := 0;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            ScaleA := Max(ScaleA, AbsComplex(A[I,J]));
            Inc(J);
        end;
        Inc(I);
    end;
    if AP_FP_Eq(ScaleA,0) then
    begin
        ScaleA := 1;
    end;
    ScaleA := 1/ScaleA;
    CMatrixLUSolveInternal(LUA, P, ScaleA, N, A, True, B, M, Info, Rep, X);
end;


(*************************************************************************
Dense solver. Same as RMatrixMixedSolve(), but for complex matrices.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* iterative refinement
* O(N^2) complexity

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    LUA     -   array[0..N-1,0..N-1], LU decomposition, CMatrixLU result
    P       -   array[0..N-1], pivots array, CMatrixLU result
    N       -   size of A
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolveM
    Rep     -   same as in RMatrixSolveM
    X       -   same as in RMatrixSolveM

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
procedure CMatrixMixedSolve(const A : TComplex2DArray;
     const LUA : TComplex2DArray;
     const P : TInteger1DArray;
     N : AlglibInteger;
     const B : TComplex1DArray;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TComplex1DArray);
var
    BM : TComplex2DArray;
    XM : TComplex2DArray;
    i_ : AlglibInteger;
begin
    if N<=0 then
    begin
        Info := -1;
        Exit;
    end;
    SetLength(BM, N, 1);
    for i_ := 0 to N-1 do
    begin
        BM[i_,0] := B[i_];
    end;
    CMatrixMixedSolveM(A, LUA, P, N, BM, 1, Info, Rep, XM);
    SetLength(X, N);
    for i_ := 0 to N-1 do
    begin
        X[i_] := XM[i_,0];
    end;
end;


(*************************************************************************
Dense solver. Same as RMatrixSolveM(), but for symmetric positive definite
matrices.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* O(N^3+M*N^2) complexity
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    IsUpper -   what half of A is provided
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve.
                Returns -3 for non-SPD matrices.
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
procedure SPDMatrixSolveM(const A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     const B : TReal2DArray;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TReal2DArray);
var
    DA : TReal2DArray;
    SqrtScaleA : Double;
    I : AlglibInteger;
    J : AlglibInteger;
    J1 : AlglibInteger;
    J2 : AlglibInteger;
begin
    
    //
    // prepare: check inputs, allocate space...
    //
    if (N<=0) or (M<=0) then
    begin
        Info := -1;
        Exit;
    end;
    SetLength(DA, N, N);
    
    //
    // 1. scale matrix, max(|A[i,j]|)
    // 2. factorize scaled matrix
    // 3. solve
    //
    SqrtScaleA := 0;
    I:=0;
    while I<=N-1 do
    begin
        if IsUpper then
        begin
            J1 := I;
            J2 := N-1;
        end
        else
        begin
            J1 := 0;
            J2 := I;
        end;
        J:=J1;
        while J<=J2 do
        begin
            SqrtScaleA := Max(SqrtScaleA, AbsReal(A[I,J]));
            Inc(J);
        end;
        Inc(I);
    end;
    if AP_FP_Eq(SqrtScaleA,0) then
    begin
        SqrtScaleA := 1;
    end;
    SqrtScaleA := 1/SqrtScaleA;
    SqrtScaleA := Sqrt(SqrtScaleA);
    I:=0;
    while I<=N-1 do
    begin
        if IsUpper then
        begin
            J1 := I;
            J2 := N-1;
        end
        else
        begin
            J1 := 0;
            J2 := I;
        end;
        APVMove(@DA[I][0], J1, J2, @A[I][0], J1, J2);
        Inc(I);
    end;
    if  not SPDMatrixCholesky(DA, N, IsUpper) then
    begin
        SetLength(X, N, M);
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                X[I,J] := 0;
                Inc(J);
            end;
            Inc(I);
        end;
        Rep.R1 := 0;
        Rep.RInf := 0;
        Info := -3;
        Exit;
    end;
    Info := 1;
    SPDMatrixCholeskySolveInternal(DA, SqrtScaleA, N, IsUpper, A, True, B, M, Info, Rep, X);
end;


(*************************************************************************
Dense solver. Same as RMatrixSolve(), but for SPD matrices.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* O(N^3) complexity
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    IsUpper -   what half of A is provided
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
                Returns -3 for non-SPD matrices.
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
procedure SPDMatrixSolve(const A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     const B : TReal1DArray;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TReal1DArray);
var
    BM : TReal2DArray;
    XM : TReal2DArray;
    i_ : AlglibInteger;
begin
    if N<=0 then
    begin
        Info := -1;
        Exit;
    end;
    SetLength(BM, N, 1);
    for i_ := 0 to N-1 do
    begin
        BM[i_,0] := B[i_];
    end;
    SPDMatrixSolveM(A, N, IsUpper, BM, 1, Info, Rep, XM);
    SetLength(X, N);
    for i_ := 0 to N-1 do
    begin
        X[i_] := XM[i_,0];
    end;
end;


(*************************************************************************
Dense solver. Same as RMatrixLUSolveM(), but for SPD matrices  represented
by their Cholesky decomposition.

Algorithm features:
* automatic detection of degenerate cases
* O(M*N^2) complexity
* condition number estimation
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    CHA     -   array[0..N-1,0..N-1], Cholesky decomposition,
                SPDMatrixCholesky result
    N       -   size of CHA
    IsUpper -   what half of CHA is provided
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
procedure SPDMatrixCholeskySolveM(const CHA : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     const B : TReal2DArray;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TReal2DArray);
var
    EmptyA : TReal2DArray;
    SqrtScaleA : Double;
    I : AlglibInteger;
    J : AlglibInteger;
    J1 : AlglibInteger;
    J2 : AlglibInteger;
begin
    
    //
    // prepare: check inputs, allocate space...
    //
    if (N<=0) or (M<=0) then
    begin
        Info := -1;
        Exit;
    end;
    
    //
    // 1. scale matrix, max(|U[i,j]|)
    // 2. factorize scaled matrix
    // 3. solve
    //
    SqrtScaleA := 0;
    I:=0;
    while I<=N-1 do
    begin
        if IsUpper then
        begin
            J1 := I;
            J2 := N-1;
        end
        else
        begin
            J1 := 0;
            J2 := I;
        end;
        J:=J1;
        while J<=J2 do
        begin
            SqrtScaleA := Max(SqrtScaleA, AbsReal(CHA[I,J]));
            Inc(J);
        end;
        Inc(I);
    end;
    if AP_FP_Eq(SqrtScaleA,0) then
    begin
        SqrtScaleA := 1;
    end;
    SqrtScaleA := 1/SqrtScaleA;
    SPDMatrixCholeskySolveInternal(CHA, SqrtScaleA, N, IsUpper, EmptyA, False, B, M, Info, Rep, X);
end;


(*************************************************************************
Dense solver. Same as RMatrixLUSolve(), but for  SPD matrices  represented
by their Cholesky decomposition.

Algorithm features:
* automatic detection of degenerate cases
* O(N^2) complexity
* condition number estimation
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    CHA     -   array[0..N-1,0..N-1], Cholesky decomposition,
                SPDMatrixCholesky result
    N       -   size of A
    IsUpper -   what half of CHA is provided
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
procedure SPDMatrixCholeskySolve(const CHA : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     const B : TReal1DArray;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TReal1DArray);
var
    BM : TReal2DArray;
    XM : TReal2DArray;
    i_ : AlglibInteger;
begin
    if N<=0 then
    begin
        Info := -1;
        Exit;
    end;
    SetLength(BM, N, 1);
    for i_ := 0 to N-1 do
    begin
        BM[i_,0] := B[i_];
    end;
    SPDMatrixCholeskySolveM(CHA, N, IsUpper, BM, 1, Info, Rep, XM);
    SetLength(X, N);
    for i_ := 0 to N-1 do
    begin
        X[i_] := XM[i_,0];
    end;
end;


(*************************************************************************
Dense solver. Same as RMatrixSolveM(), but for Hermitian positive definite
matrices.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* O(N^3+M*N^2) complexity
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    IsUpper -   what half of A is provided
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve.
                Returns -3 for non-HPD matrices.
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
procedure HPDMatrixSolveM(const A : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     const B : TComplex2DArray;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TComplex2DArray);
var
    DA : TComplex2DArray;
    SqrtScaleA : Double;
    I : AlglibInteger;
    J : AlglibInteger;
    J1 : AlglibInteger;
    J2 : AlglibInteger;
    i_ : AlglibInteger;
begin
    
    //
    // prepare: check inputs, allocate space...
    //
    if (N<=0) or (M<=0) then
    begin
        Info := -1;
        Exit;
    end;
    SetLength(DA, N, N);
    
    //
    // 1. scale matrix, max(|A[i,j]|)
    // 2. factorize scaled matrix
    // 3. solve
    //
    SqrtScaleA := 0;
    I:=0;
    while I<=N-1 do
    begin
        if IsUpper then
        begin
            J1 := I;
            J2 := N-1;
        end
        else
        begin
            J1 := 0;
            J2 := I;
        end;
        J:=J1;
        while J<=J2 do
        begin
            SqrtScaleA := Max(SqrtScaleA, AbsComplex(A[I,J]));
            Inc(J);
        end;
        Inc(I);
    end;
    if AP_FP_Eq(SqrtScaleA,0) then
    begin
        SqrtScaleA := 1;
    end;
    SqrtScaleA := 1/SqrtScaleA;
    SqrtScaleA := Sqrt(SqrtScaleA);
    I:=0;
    while I<=N-1 do
    begin
        if IsUpper then
        begin
            J1 := I;
            J2 := N-1;
        end
        else
        begin
            J1 := 0;
            J2 := I;
        end;
        for i_ := J1 to J2 do
        begin
            DA[I,i_] := A[I,i_];
        end;
        Inc(I);
    end;
    if  not HPDMatrixCholesky(DA, N, IsUpper) then
    begin
        SetLength(X, N, M);
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                X[I,J] := C_Complex(0);
                Inc(J);
            end;
            Inc(I);
        end;
        Rep.R1 := 0;
        Rep.RInf := 0;
        Info := -3;
        Exit;
    end;
    Info := 1;
    HPDMatrixCholeskySolveInternal(DA, SqrtScaleA, N, IsUpper, A, True, B, M, Info, Rep, X);
end;


(*************************************************************************
Dense solver. Same as RMatrixSolve(),  but for Hermitian positive definite
matrices.

Algorithm features:
* automatic detection of degenerate cases
* condition number estimation
* O(N^3) complexity
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    IsUpper -   what half of A is provided
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
                Returns -3 for non-HPD matrices.
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
procedure HPDMatrixSolve(const A : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     const B : TComplex1DArray;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TComplex1DArray);
var
    BM : TComplex2DArray;
    XM : TComplex2DArray;
    i_ : AlglibInteger;
begin
    if N<=0 then
    begin
        Info := -1;
        Exit;
    end;
    SetLength(BM, N, 1);
    for i_ := 0 to N-1 do
    begin
        BM[i_,0] := B[i_];
    end;
    HPDMatrixSolveM(A, N, IsUpper, BM, 1, Info, Rep, XM);
    SetLength(X, N);
    for i_ := 0 to N-1 do
    begin
        X[i_] := XM[i_,0];
    end;
end;


(*************************************************************************
Dense solver. Same as RMatrixLUSolveM(), but for HPD matrices  represented
by their Cholesky decomposition.

Algorithm features:
* automatic detection of degenerate cases
* O(M*N^2) complexity
* condition number estimation
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    CHA     -   array[0..N-1,0..N-1], Cholesky decomposition,
                HPDMatrixCholesky result
    N       -   size of CHA
    IsUpper -   what half of CHA is provided
    B       -   array[0..N-1,0..M-1], right part
    M       -   right part size

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
procedure HPDMatrixCholeskySolveM(const CHA : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     const B : TComplex2DArray;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TComplex2DArray);
var
    EmptyA : TComplex2DArray;
    SqrtScaleA : Double;
    I : AlglibInteger;
    J : AlglibInteger;
    J1 : AlglibInteger;
    J2 : AlglibInteger;
begin
    
    //
    // prepare: check inputs, allocate space...
    //
    if (N<=0) or (M<=0) then
    begin
        Info := -1;
        Exit;
    end;
    
    //
    // 1. scale matrix, max(|U[i,j]|)
    // 2. factorize scaled matrix
    // 3. solve
    //
    SqrtScaleA := 0;
    I:=0;
    while I<=N-1 do
    begin
        if IsUpper then
        begin
            J1 := I;
            J2 := N-1;
        end
        else
        begin
            J1 := 0;
            J2 := I;
        end;
        J:=J1;
        while J<=J2 do
        begin
            SqrtScaleA := Max(SqrtScaleA, AbsComplex(CHA[I,J]));
            Inc(J);
        end;
        Inc(I);
    end;
    if AP_FP_Eq(SqrtScaleA,0) then
    begin
        SqrtScaleA := 1;
    end;
    SqrtScaleA := 1/SqrtScaleA;
    HPDMatrixCholeskySolveInternal(CHA, SqrtScaleA, N, IsUpper, EmptyA, False, B, M, Info, Rep, X);
end;


(*************************************************************************
Dense solver. Same as RMatrixLUSolve(), but for  HPD matrices  represented
by their Cholesky decomposition.

Algorithm features:
* automatic detection of degenerate cases
* O(N^2) complexity
* condition number estimation
* matrix is represented by its upper or lower triangle

No iterative refinement is provided because such partial representation of
matrix does not allow efficient calculation of extra-precise  matrix-vector
products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
need iterative refinement.

INPUT PARAMETERS
    CHA     -   array[0..N-1,0..N-1], Cholesky decomposition,
                SPDMatrixCholesky result
    N       -   size of A
    IsUpper -   what half of CHA is provided
    B       -   array[0..N-1], right part

OUTPUT PARAMETERS
    Info    -   same as in RMatrixSolve
    Rep     -   same as in RMatrixSolve
    X       -   same as in RMatrixSolve

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
procedure HPDMatrixCholeskySolve(const CHA : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     const B : TComplex1DArray;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TComplex1DArray);
var
    BM : TComplex2DArray;
    XM : TComplex2DArray;
    i_ : AlglibInteger;
begin
    if N<=0 then
    begin
        Info := -1;
        Exit;
    end;
    SetLength(BM, N, 1);
    for i_ := 0 to N-1 do
    begin
        BM[i_,0] := B[i_];
    end;
    HPDMatrixCholeskySolveM(CHA, N, IsUpper, BM, 1, Info, Rep, XM);
    SetLength(X, N);
    for i_ := 0 to N-1 do
    begin
        X[i_] := XM[i_,0];
    end;
end;


(*************************************************************************
Dense solver.

This subroutine finds solution of the linear system A*X=B with non-square,
possibly degenerate A.  System  is  solved in the least squares sense, and
general least squares solution  X = X0 + CX*y  which  minimizes |A*X-B| is
returned. If A is non-degenerate, solution in the  usual sense is returned

Algorithm features:
* automatic detection of degenerate cases
* iterative refinement
* O(N^3) complexity

INPUT PARAMETERS
    A       -   array[0..NRows-1,0..NCols-1], system matrix
    NRows   -   vertical size of A
    NCols   -   horizontal size of A
    B       -   array[0..NCols-1], right part
    Threshold-  a number in [0,1]. Singular values  beyond  Threshold  are
                considered  zero.  Set  it to 0.0, if you don't understand
                what it means, so the solver will choose good value on its
                own.
                
OUTPUT PARAMETERS
    Info    -   return code:
                * -4    SVD subroutine failed
                * -1    if NRows<=0 or NCols<=0 or Threshold<0 was passed
                *  1    if task is solved
    Rep     -   solver report, see below for more info
    X       -   array[0..N-1,0..M-1], it contains:
                * solution of A*X=B if A is non-singular (well-conditioned
                  or ill-conditioned, but not very close to singular)
                * zeros,  if  A  is  singular  or  VERY  close to singular
                  (in this case Info=-3).

SOLVER REPORT

Subroutine sets following fields of the Rep structure:
* R2        reciprocal of condition number: 1/cond(A), 2-norm.
* N         = NCols
* K         dim(Null(A))
* CX        array[0..N-1,0..K-1], kernel of A.
            Columns of CX store such vectors that A*CX[i]=0.

  -- ALGLIB --
     Copyright 24.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure RMatrixSolveLS(const A : TReal2DArray;
     NRows : AlglibInteger;
     NCols : AlglibInteger;
     const B : TReal1DArray;
     Threshold : Double;
     var Info : AlglibInteger;
     var Rep : DenseSolverLSReport;
     var X : TReal1DArray);
var
    SV : TReal1DArray;
    U : TReal2DArray;
    VT : TReal2DArray;
    RP : TReal1DArray;
    UTB : TReal1DArray;
    SUTB : TReal1DArray;
    Tmp : TReal1DArray;
    TA : TReal1DArray;
    TX : TReal1DArray;
    Buf : TReal1DArray;
    W : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    NSV : AlglibInteger;
    KernelIdx : AlglibInteger;
    V : Double;
    VErr : Double;
    SVDFailed : Boolean;
    ZeroA : Boolean;
    RFS : AlglibInteger;
    NRFS : AlglibInteger;
    TerminateNextTime : Boolean;
    SmallErr : Boolean;
    i_ : AlglibInteger;
begin
    if (NRows<=0) or (NCols<=0) or AP_FP_Less(Threshold,0) then
    begin
        Info := -1;
        Exit;
    end;
    if AP_FP_Eq(Threshold,0) then
    begin
        Threshold := 1000*MachineEpsilon;
    end;
    
    //
    // Factorize A first
    //
    SVDFailed :=  not RMatrixSVD(A, NRows, NCols, 1, 2, 2, SV, U, VT);
    ZeroA := AP_FP_Eq(SV[0],0);
    if SVDFailed or ZeroA then
    begin
        if SVDFailed then
        begin
            Info := -4;
        end
        else
        begin
            Info := 1;
        end;
        SetLength(X, NCols);
        I:=0;
        while I<=NCols-1 do
        begin
            X[I] := 0;
            Inc(I);
        end;
        Rep.N := NCols;
        Rep.K := NCols;
        SetLength(Rep.CX, NCols, NCols);
        I:=0;
        while I<=NCols-1 do
        begin
            J:=0;
            while J<=NCols-1 do
            begin
                if I=J then
                begin
                    Rep.CX[I,J] := 1;
                end
                else
                begin
                    Rep.CX[I,J] := 0;
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        Rep.R2 := 0;
        Exit;
    end;
    NSV := Min(NCols, NRows);
    if NSV=NCols then
    begin
        Rep.R2 := SV[NSV-1]/SV[0];
    end
    else
    begin
        Rep.R2 := 0;
    end;
    Rep.N := NCols;
    Info := 1;
    
    //
    // Iterative refinement of xc combined with solution:
    // 1. xc = 0
    // 2. calculate r = bc-A*xc using extra-precise dot product
    // 3. solve A*y = r
    // 4. update x:=x+r
    // 5. goto 2
    //
    // This cycle is executed until one of two things happens:
    // 1. maximum number of iterations reached
    // 2. last iteration decreased error to the lower limit
    //
    SetLength(UTB, NSV);
    SetLength(SUTB, NSV);
    SetLength(X, NCols);
    SetLength(Tmp, NCols);
    SetLength(TA, NCols+1);
    SetLength(TX, NCols+1);
    SetLength(Buf, NCols+1);
    I:=0;
    while I<=NCols-1 do
    begin
        X[I] := 0;
        Inc(I);
    end;
    KernelIdx := NSV;
    I:=0;
    while I<=NSV-1 do
    begin
        if AP_FP_Less_Eq(SV[I],Threshold*SV[0]) then
        begin
            KernelIdx := I;
            Break;
        end;
        Inc(I);
    end;
    Rep.K := NCols-KernelIdx;
    NRFS := DenseSolverRFSMaxV2(NCols, Rep.R2);
    TerminateNextTime := False;
    SetLength(RP, NRows);
    RFS:=0;
    while RFS<=NRFS do
    begin
        if TerminateNextTime then
        begin
            Break;
        end;
        
        //
        // calculate right part
        //
        if RFS=0 then
        begin
            APVMove(@RP[0], 0, NRows-1, @B[0], 0, NRows-1);
        end
        else
        begin
            SmallErr := True;
            I:=0;
            while I<=NRows-1 do
            begin
                APVMove(@TA[0], 0, NCols-1, @A[I][0], 0, NCols-1);
                TA[NCols] := -1;
                APVMove(@TX[0], 0, NCols-1, @X[0], 0, NCols-1);
                TX[NCols] := B[I];
                XDot(TA, TX, NCols+1, Buf, V, VErr);
                RP[I] := -V;
                SmallErr := SmallErr and AP_FP_Less(AbsReal(V),4*VErr);
                Inc(I);
            end;
            if SmallErr then
            begin
                TerminateNextTime := True;
            end;
        end;
        
        //
        // solve A*dx = rp
        //
        I:=0;
        while I<=NCols-1 do
        begin
            Tmp[I] := 0;
            Inc(I);
        end;
        I:=0;
        while I<=NSV-1 do
        begin
            UTB[I] := 0;
            Inc(I);
        end;
        I:=0;
        while I<=NRows-1 do
        begin
            V := RP[I];
            APVAdd(@UTB[0], 0, NSV-1, @U[I][0], 0, NSV-1, V);
            Inc(I);
        end;
        I:=0;
        while I<=NSV-1 do
        begin
            if I<KernelIdx then
            begin
                SUTB[I] := UTB[I]/SV[I];
            end
            else
            begin
                SUTB[I] := 0;
            end;
            Inc(I);
        end;
        I:=0;
        while I<=NSV-1 do
        begin
            V := SUTB[I];
            APVAdd(@Tmp[0], 0, NCols-1, @VT[I][0], 0, NCols-1, V);
            Inc(I);
        end;
        
        //
        // update x:  x:=x+dx
        //
        APVAdd(@X[0], 0, NCols-1, @Tmp[0], 0, NCols-1);
        Inc(RFS);
    end;
    
    //
    // fill CX
    //
    if Rep.K>0 then
    begin
        SetLength(Rep.CX, NCols, Rep.K);
        I:=0;
        while I<=Rep.K-1 do
        begin
            for i_ := 0 to NCols-1 do
            begin
                Rep.CX[i_,I] := VT[KernelIdx+I,i_];
            end;
            Inc(I);
        end;
    end;
end;


(*************************************************************************
Internal LU solver

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
procedure RMatrixLUSolveInternal(const LUA : TReal2DArray;
     const P : TInteger1DArray;
     const ScaleA : Double;
     N : AlglibInteger;
     const A : TReal2DArray;
     HaveA : Boolean;
     const B : TReal2DArray;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    RFS : AlglibInteger;
    NRFS : AlglibInteger;
    XC : TReal1DArray;
    Y : TReal1DArray;
    BC : TReal1DArray;
    XA : TReal1DArray;
    XB : TReal1DArray;
    TX : TReal1DArray;
    V : Double;
    VErr : Double;
    MXB : Double;
    ScaleRight : Double;
    SmallErr : Boolean;
    TerminateNextTime : Boolean;
    i_ : AlglibInteger;
begin
    Assert(AP_FP_Greater(ScaleA,0));
    
    //
    // prepare: check inputs, allocate space...
    //
    if (N<=0) or (M<=0) then
    begin
        Info := -1;
        Exit;
    end;
    I:=0;
    while I<=N-1 do
    begin
        if (P[I]>N-1) or (P[I]<I) then
        begin
            Info := -1;
            Exit;
        end;
        Inc(I);
    end;
    SetLength(X, N, M);
    SetLength(Y, N);
    SetLength(XC, N);
    SetLength(BC, N);
    SetLength(TX, N+1);
    SetLength(XA, N+1);
    SetLength(XB, N+1);
    
    //
    // estimate condition number, test for near singularity
    //
    Rep.R1 := RMatrixLURCond1(LUA, N);
    Rep.RInf := RMatrixLURCondInf(LUA, N);
    if AP_FP_Less(Rep.R1,RCondThreshold) or AP_FP_Less(Rep.RInf,RCondThreshold) then
    begin
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                X[I,J] := 0;
                Inc(J);
            end;
            Inc(I);
        end;
        Rep.R1 := 0;
        Rep.RInf := 0;
        Info := -3;
        Exit;
    end;
    Info := 1;
    
    //
    // solve
    //
    K:=0;
    while K<=M-1 do
    begin
        
        //
        // copy B to contiguous storage
        //
        for i_ := 0 to N-1 do
        begin
            BC[i_] := B[i_,K];
        end;
        
        //
        // Scale right part:
        // * MX stores max(|Bi|)
        // * ScaleRight stores actual scaling applied to B when solving systems
        //   it is chosen to make |scaleRight*b| close to 1.
        //
        MXB := 0;
        I:=0;
        while I<=N-1 do
        begin
            MXB := Max(MXB, AbsReal(BC[I]));
            Inc(I);
        end;
        if AP_FP_Eq(MXB,0) then
        begin
            MXB := 1;
        end;
        ScaleRight := 1/MXB;
        
        //
        // First, non-iterative part of solution process.
        // We use separate code for this task because
        // XDot is quite slow and we want to save time.
        //
        APVMove(@XC[0], 0, N-1, @BC[0], 0, N-1, ScaleRight);
        RBasicLUSolve(LUA, P, ScaleA, N, XC, TX);
        
        //
        // Iterative refinement of xc:
        // * calculate r = bc-A*xc using extra-precise dot product
        // * solve A*y = r
        // * update x:=x+r
        //
        // This cycle is executed until one of two things happens:
        // 1. maximum number of iterations reached
        // 2. last iteration decreased error to the lower limit
        //
        if HaveA then
        begin
            NRFS := DenseSolverRFSMax(N, Rep.R1, Rep.RInf);
            TerminateNextTime := False;
            RFS:=0;
            while RFS<=NRFS-1 do
            begin
                if TerminateNextTime then
                begin
                    Break;
                end;
                
                //
                // generate right part
                //
                SmallErr := True;
                APVMove(@XB[0], 0, N-1, @XC[0], 0, N-1);
                I:=0;
                while I<=N-1 do
                begin
                    APVMove(@XA[0], 0, N-1, @A[I][0], 0, N-1, ScaleA);
                    XA[N] := -1;
                    XB[N] := ScaleRight*BC[I];
                    XDot(XA, XB, N+1, TX, V, VErr);
                    Y[I] := -V;
                    SmallErr := SmallErr and AP_FP_Less(AbsReal(V),4*VErr);
                    Inc(I);
                end;
                if SmallErr then
                begin
                    TerminateNextTime := True;
                end;
                
                //
                // solve and update
                //
                RBasicLUSolve(LUA, P, ScaleA, N, Y, TX);
                APVAdd(@XC[0], 0, N-1, @Y[0], 0, N-1);
                Inc(RFS);
            end;
        end;
        
        //
        // Store xc.
        // Post-scale result.
        //
        V := ScaleA*MXB;
        for i_ := 0 to N-1 do
        begin
            X[i_,K] := V*XC[i_];
        end;
        Inc(K);
    end;
end;


(*************************************************************************
Internal Cholesky solver

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
procedure SPDMatrixCholeskySolveInternal(const CHA : TReal2DArray;
     const SqrtScaleA : Double;
     N : AlglibInteger;
     IsUpper : Boolean;
     const A : TReal2DArray;
     HaveA : Boolean;
     const B : TReal2DArray;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    RFS : AlglibInteger;
    NRFS : AlglibInteger;
    XC : TReal1DArray;
    Y : TReal1DArray;
    BC : TReal1DArray;
    XA : TReal1DArray;
    XB : TReal1DArray;
    TX : TReal1DArray;
    V : Double;
    VErr : Double;
    MXB : Double;
    ScaleRight : Double;
    SmallErr : Boolean;
    TerminateNextTime : Boolean;
    i_ : AlglibInteger;
begin
    Assert(AP_FP_Greater(SqrtScaleA,0));
    
    //
    // prepare: check inputs, allocate space...
    //
    if (N<=0) or (M<=0) then
    begin
        Info := -1;
        Exit;
    end;
    SetLength(X, N, M);
    SetLength(Y, N);
    SetLength(XC, N);
    SetLength(BC, N);
    SetLength(TX, N+1);
    SetLength(XA, N+1);
    SetLength(XB, N+1);
    
    //
    // estimate condition number, test for near singularity
    //
    Rep.R1 := SPDMatrixCholeskyRCond(CHA, N, IsUpper);
    Rep.RInf := Rep.R1;
    if AP_FP_Less(Rep.R1,RCondThreshold) then
    begin
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                X[I,J] := 0;
                Inc(J);
            end;
            Inc(I);
        end;
        Rep.R1 := 0;
        Rep.RInf := 0;
        Info := -3;
        Exit;
    end;
    Info := 1;
    
    //
    // solve
    //
    K:=0;
    while K<=M-1 do
    begin
        
        //
        // copy B to contiguous storage
        //
        for i_ := 0 to N-1 do
        begin
            BC[i_] := B[i_,K];
        end;
        
        //
        // Scale right part:
        // * MX stores max(|Bi|)
        // * ScaleRight stores actual scaling applied to B when solving systems
        //   it is chosen to make |scaleRight*b| close to 1.
        //
        MXB := 0;
        I:=0;
        while I<=N-1 do
        begin
            MXB := Max(MXB, AbsReal(BC[I]));
            Inc(I);
        end;
        if AP_FP_Eq(MXB,0) then
        begin
            MXB := 1;
        end;
        ScaleRight := 1/MXB;
        
        //
        // First, non-iterative part of solution process.
        // We use separate code for this task because
        // XDot is quite slow and we want to save time.
        //
        APVMove(@XC[0], 0, N-1, @BC[0], 0, N-1, ScaleRight);
        SPDBasicCholeskySolve(CHA, SqrtScaleA, N, IsUpper, XC, TX);
        
        //
        // Store xc.
        // Post-scale result.
        //
        V := AP_Sqr(SqrtScaleA)*MXB;
        for i_ := 0 to N-1 do
        begin
            X[i_,K] := V*XC[i_];
        end;
        Inc(K);
    end;
end;


(*************************************************************************
Internal LU solver

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
procedure CMatrixLUSolveInternal(const LUA : TComplex2DArray;
     const P : TInteger1DArray;
     const ScaleA : Double;
     N : AlglibInteger;
     const A : TComplex2DArray;
     HaveA : Boolean;
     const B : TComplex2DArray;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TComplex2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    RFS : AlglibInteger;
    NRFS : AlglibInteger;
    XC : TComplex1DArray;
    Y : TComplex1DArray;
    BC : TComplex1DArray;
    XA : TComplex1DArray;
    XB : TComplex1DArray;
    TX : TComplex1DArray;
    TmpBuf : TReal1DArray;
    V : Complex;
    VErr : Double;
    MXB : Double;
    ScaleRight : Double;
    SmallErr : Boolean;
    TerminateNextTime : Boolean;
    i_ : AlglibInteger;
begin
    Assert(AP_FP_Greater(ScaleA,0));
    
    //
    // prepare: check inputs, allocate space...
    //
    if (N<=0) or (M<=0) then
    begin
        Info := -1;
        Exit;
    end;
    I:=0;
    while I<=N-1 do
    begin
        if (P[I]>N-1) or (P[I]<I) then
        begin
            Info := -1;
            Exit;
        end;
        Inc(I);
    end;
    SetLength(X, N, M);
    SetLength(Y, N);
    SetLength(XC, N);
    SetLength(BC, N);
    SetLength(TX, N);
    SetLength(XA, N+1);
    SetLength(XB, N+1);
    SetLength(TmpBuf, 2*N+2);
    
    //
    // estimate condition number, test for near singularity
    //
    Rep.R1 := CMatrixLURCond1(LUA, N);
    Rep.RInf := CMatrixLURCondInf(LUA, N);
    if AP_FP_Less(Rep.R1,RCondThreshold) or AP_FP_Less(Rep.RInf,RCondThreshold) then
    begin
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                X[I,J] := C_Complex(0);
                Inc(J);
            end;
            Inc(I);
        end;
        Rep.R1 := 0;
        Rep.RInf := 0;
        Info := -3;
        Exit;
    end;
    Info := 1;
    
    //
    // solve
    //
    K:=0;
    while K<=M-1 do
    begin
        
        //
        // copy B to contiguous storage
        //
        for i_ := 0 to N-1 do
        begin
            BC[i_] := B[i_,K];
        end;
        
        //
        // Scale right part:
        // * MX stores max(|Bi|)
        // * ScaleRight stores actual scaling applied to B when solving systems
        //   it is chosen to make |scaleRight*b| close to 1.
        //
        MXB := 0;
        I:=0;
        while I<=N-1 do
        begin
            MXB := Max(MXB, AbsComplex(BC[I]));
            Inc(I);
        end;
        if AP_FP_Eq(MXB,0) then
        begin
            MXB := 1;
        end;
        ScaleRight := 1/MXB;
        
        //
        // First, non-iterative part of solution process.
        // We use separate code for this task because
        // XDot is quite slow and we want to save time.
        //
        for i_ := 0 to N-1 do
        begin
            XC[i_] := C_MulR(BC[i_],ScaleRight);
        end;
        CBasicLUSolve(LUA, P, ScaleA, N, XC, TX);
        
        //
        // Iterative refinement of xc:
        // * calculate r = bc-A*xc using extra-precise dot product
        // * solve A*y = r
        // * update x:=x+r
        //
        // This cycle is executed until one of two things happens:
        // 1. maximum number of iterations reached
        // 2. last iteration decreased error to the lower limit
        //
        if HaveA then
        begin
            NRFS := DenseSolverRFSMax(N, Rep.R1, Rep.RInf);
            TerminateNextTime := False;
            RFS:=0;
            while RFS<=NRFS-1 do
            begin
                if TerminateNextTime then
                begin
                    Break;
                end;
                
                //
                // generate right part
                //
                SmallErr := True;
                for i_ := 0 to N-1 do
                begin
                    XB[i_] := XC[i_];
                end;
                I:=0;
                while I<=N-1 do
                begin
                    for i_ := 0 to N-1 do
                    begin
                        XA[i_] := C_MulR(A[I,i_],ScaleA);
                    end;
                    XA[N] := C_Complex(-1);
                    XB[N] := C_MulR(BC[I],ScaleRight);
                    XCDot(XA, XB, N+1, TmpBuf, V, VErr);
                    Y[I] := C_Opposite(V);
                    SmallErr := SmallErr and AP_FP_Less(AbsComplex(V),4*VErr);
                    Inc(I);
                end;
                if SmallErr then
                begin
                    TerminateNextTime := True;
                end;
                
                //
                // solve and update
                //
                CBasicLUSolve(LUA, P, ScaleA, N, Y, TX);
                for i_ := 0 to N-1 do
                begin
                    XC[i_] := C_Add(XC[i_], Y[i_]);
                end;
                Inc(RFS);
            end;
        end;
        
        //
        // Store xc.
        // Post-scale result.
        //
        V := C_Complex(ScaleA*MXB);
        for i_ := 0 to N-1 do
        begin
            X[i_,K] := C_Mul(V, XC[i_]);
        end;
        Inc(K);
    end;
end;


(*************************************************************************
Internal Cholesky solver

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
procedure HPDMatrixCholeskySolveInternal(const CHA : TComplex2DArray;
     const SqrtScaleA : Double;
     N : AlglibInteger;
     IsUpper : Boolean;
     const A : TComplex2DArray;
     HaveA : Boolean;
     const B : TComplex2DArray;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TComplex2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    RFS : AlglibInteger;
    NRFS : AlglibInteger;
    XC : TComplex1DArray;
    Y : TComplex1DArray;
    BC : TComplex1DArray;
    XA : TComplex1DArray;
    XB : TComplex1DArray;
    TX : TComplex1DArray;
    V : Double;
    VErr : Double;
    MXB : Double;
    ScaleRight : Double;
    SmallErr : Boolean;
    TerminateNextTime : Boolean;
    i_ : AlglibInteger;
begin
    Assert(AP_FP_Greater(SqrtScaleA,0));
    
    //
    // prepare: check inputs, allocate space...
    //
    if (N<=0) or (M<=0) then
    begin
        Info := -1;
        Exit;
    end;
    SetLength(X, N, M);
    SetLength(Y, N);
    SetLength(XC, N);
    SetLength(BC, N);
    SetLength(TX, N+1);
    SetLength(XA, N+1);
    SetLength(XB, N+1);
    
    //
    // estimate condition number, test for near singularity
    //
    Rep.R1 := HPDMatrixCholeskyRCond(CHA, N, IsUpper);
    Rep.RInf := Rep.R1;
    if AP_FP_Less(Rep.R1,RCondThreshold) then
    begin
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                X[I,J] := C_Complex(0);
                Inc(J);
            end;
            Inc(I);
        end;
        Rep.R1 := 0;
        Rep.RInf := 0;
        Info := -3;
        Exit;
    end;
    Info := 1;
    
    //
    // solve
    //
    K:=0;
    while K<=M-1 do
    begin
        
        //
        // copy B to contiguous storage
        //
        for i_ := 0 to N-1 do
        begin
            BC[i_] := B[i_,K];
        end;
        
        //
        // Scale right part:
        // * MX stores max(|Bi|)
        // * ScaleRight stores actual scaling applied to B when solving systems
        //   it is chosen to make |scaleRight*b| close to 1.
        //
        MXB := 0;
        I:=0;
        while I<=N-1 do
        begin
            MXB := Max(MXB, AbsComplex(BC[I]));
            Inc(I);
        end;
        if AP_FP_Eq(MXB,0) then
        begin
            MXB := 1;
        end;
        ScaleRight := 1/MXB;
        
        //
        // First, non-iterative part of solution process.
        // We use separate code for this task because
        // XDot is quite slow and we want to save time.
        //
        for i_ := 0 to N-1 do
        begin
            XC[i_] := C_MulR(BC[i_],ScaleRight);
        end;
        HPDBasicCholeskySolve(CHA, SqrtScaleA, N, IsUpper, XC, TX);
        
        //
        // Store xc.
        // Post-scale result.
        //
        V := AP_Sqr(SqrtScaleA)*MXB;
        for i_ := 0 to N-1 do
        begin
            X[i_,K] := C_MulR(XC[i_],V);
        end;
        Inc(K);
    end;
end;


(*************************************************************************
Internal subroutine.
Returns maximum count of RFS iterations as function of:
1. machine epsilon
2. task size.
3. condition number

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
function DenseSolverRFSMax(N : AlglibInteger;
     R1 : Double;
     RInf : Double):AlglibInteger;
begin
    Result := 5;
end;


(*************************************************************************
Internal subroutine.
Returns maximum count of RFS iterations as function of:
1. machine epsilon
2. task size.
3. norm-2 condition number

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
function DenseSolverRFSMaxV2(N : AlglibInteger; R2 : Double):AlglibInteger;
begin
    Result := DenseSolverRFSMax(N, 0, 0);
end;


(*************************************************************************
Basic LU solver for ScaleA*PLU*x = y.

This subroutine assumes that:
* L is well-scaled, and it is U which needs scaling by ScaleA.
* A=PLU is well-conditioned, so no zero divisions or overflow may occur

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
procedure RBasicLUSolve(const LUA : TReal2DArray;
     const P : TInteger1DArray;
     ScaleA : Double;
     N : AlglibInteger;
     var XB : TReal1DArray;
     var Tmp : TReal1DArray);
var
    I : AlglibInteger;
    V : Double;
begin
    I:=0;
    while I<=N-1 do
    begin
        if P[I]<>I then
        begin
            V := XB[I];
            XB[I] := XB[P[I]];
            XB[P[I]] := V;
        end;
        Inc(I);
    end;
    I:=1;
    while I<=N-1 do
    begin
        V := APVDotProduct(@LUA[I][0], 0, I-1, @XB[0], 0, I-1);
        XB[I] := XB[I]-V;
        Inc(I);
    end;
    XB[N-1] := XB[N-1]/(ScaleA*LUA[N-1,N-1]);
    I:=N-2;
    while I>=0 do
    begin
        APVMove(@Tmp[0], I+1, N-1, @LUA[I][0], I+1, N-1, ScaleA);
        V := APVDotProduct(@Tmp[0], I+1, N-1, @XB[0], I+1, N-1);
        XB[I] := (XB[I]-V)/(ScaleA*LUA[I,I]);
        Dec(I);
    end;
end;


(*************************************************************************
Basic Cholesky solver for ScaleA*Cholesky(A)'*x = y.

This subroutine assumes that:
* A*ScaleA is well scaled
* A is well-conditioned, so no zero divisions or overflow may occur

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
procedure SPDBasicCholeskySolve(const CHA : TReal2DArray;
     SqrtScaleA : Double;
     N : AlglibInteger;
     IsUpper : Boolean;
     var XB : TReal1DArray;
     var Tmp : TReal1DArray);
var
    I : AlglibInteger;
    V : Double;
begin
    
    //
    // A = L*L' or A=U'*U
    //
    if IsUpper then
    begin
        
        //
        // Solve U'*y=b first.
        //
        I:=0;
        while I<=N-1 do
        begin
            XB[I] := XB[I]/(SqrtScaleA*CHA[I,I]);
            if I<N-1 then
            begin
                V := XB[I];
                APVMove(@Tmp[0], I+1, N-1, @CHA[I][0], I+1, N-1, SqrtScaleA);
                APVSub(@XB[0], I+1, N-1, @Tmp[0], I+1, N-1, V);
            end;
            Inc(I);
        end;
        
        //
        // Solve U*x=y then.
        //
        I:=N-1;
        while I>=0 do
        begin
            if I<N-1 then
            begin
                APVMove(@Tmp[0], I+1, N-1, @CHA[I][0], I+1, N-1, SqrtScaleA);
                V := APVDotProduct(@Tmp[0], I+1, N-1, @XB[0], I+1, N-1);
                XB[I] := XB[I]-V;
            end;
            XB[I] := XB[I]/(SqrtScaleA*CHA[I,I]);
            Dec(I);
        end;
    end
    else
    begin
        
        //
        // Solve L*y=b first
        //
        I:=0;
        while I<=N-1 do
        begin
            if I>0 then
            begin
                APVMove(@Tmp[0], 0, I-1, @CHA[I][0], 0, I-1, SqrtScaleA);
                V := APVDotProduct(@Tmp[0], 0, I-1, @XB[0], 0, I-1);
                XB[I] := XB[I]-V;
            end;
            XB[I] := XB[I]/(SqrtScaleA*CHA[I,I]);
            Inc(I);
        end;
        
        //
        // Solve L'*x=y then.
        //
        I:=N-1;
        while I>=0 do
        begin
            XB[I] := XB[I]/(SqrtScaleA*CHA[I,I]);
            if I>0 then
            begin
                V := XB[I];
                APVMove(@Tmp[0], 0, I-1, @CHA[I][0], 0, I-1, SqrtScaleA);
                APVSub(@XB[0], 0, I-1, @Tmp[0], 0, I-1, V);
            end;
            Dec(I);
        end;
    end;
end;


(*************************************************************************
Basic LU solver for ScaleA*PLU*x = y.

This subroutine assumes that:
* L is well-scaled, and it is U which needs scaling by ScaleA.
* A=PLU is well-conditioned, so no zero divisions or overflow may occur

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
procedure CBasicLUSolve(const LUA : TComplex2DArray;
     const P : TInteger1DArray;
     ScaleA : Double;
     N : AlglibInteger;
     var XB : TComplex1DArray;
     var Tmp : TComplex1DArray);
var
    I : AlglibInteger;
    V : Complex;
    i_ : AlglibInteger;
begin
    I:=0;
    while I<=N-1 do
    begin
        if P[I]<>I then
        begin
            V := XB[I];
            XB[I] := XB[P[I]];
            XB[P[I]] := V;
        end;
        Inc(I);
    end;
    I:=1;
    while I<=N-1 do
    begin
        V := C_Complex(0.0);
        for i_ := 0 to I-1 do
        begin
            V := C_Add(V,C_Mul(LUA[I,i_],XB[i_]));
        end;
        XB[I] := C_Sub(XB[I],V);
        Inc(I);
    end;
    XB[N-1] := C_Div(XB[N-1],C_MulR(LUA[N-1,N-1],ScaleA));
    I:=N-2;
    while I>=0 do
    begin
        for i_ := I+1 to N-1 do
        begin
            Tmp[i_] := C_MulR(LUA[I,i_],ScaleA);
        end;
        V := C_Complex(0.0);
        for i_ := I+1 to N-1 do
        begin
            V := C_Add(V,C_Mul(Tmp[i_],XB[i_]));
        end;
        XB[I] := C_Div(C_Sub(XB[I],V),C_MulR(LUA[I,I],ScaleA));
        Dec(I);
    end;
end;


(*************************************************************************
Basic Cholesky solver for ScaleA*Cholesky(A)'*x = y.

This subroutine assumes that:
* A*ScaleA is well scaled
* A is well-conditioned, so no zero divisions or overflow may occur

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************)
procedure HPDBasicCholeskySolve(const CHA : TComplex2DArray;
     SqrtScaleA : Double;
     N : AlglibInteger;
     IsUpper : Boolean;
     var XB : TComplex1DArray;
     var Tmp : TComplex1DArray);
var
    I : AlglibInteger;
    V : Complex;
    i_ : AlglibInteger;
begin
    
    //
    // A = L*L' or A=U'*U
    //
    if IsUpper then
    begin
        
        //
        // Solve U'*y=b first.
        //
        I:=0;
        while I<=N-1 do
        begin
            XB[I] := C_Div(XB[I],C_MulR(Conj(CHA[I,I]),SqrtScaleA));
            if I<N-1 then
            begin
                V := XB[I];
                for i_ := I+1 to N-1 do
                begin
                    Tmp[i_] := C_MulR(Conj(CHA[I,i_]),SqrtScaleA);
                end;
                for i_ := I+1 to N-1 do
                begin
                    XB[i_] := C_Sub(XB[i_], C_Mul(V, Tmp[i_]));
                end;
            end;
            Inc(I);
        end;
        
        //
        // Solve U*x=y then.
        //
        I:=N-1;
        while I>=0 do
        begin
            if I<N-1 then
            begin
                for i_ := I+1 to N-1 do
                begin
                    Tmp[i_] := C_MulR(CHA[I,i_],SqrtScaleA);
                end;
                V := C_Complex(0.0);
                for i_ := I+1 to N-1 do
                begin
                    V := C_Add(V,C_Mul(Tmp[i_],XB[i_]));
                end;
                XB[I] := C_Sub(XB[I],V);
            end;
            XB[I] := C_Div(XB[I],C_MulR(CHA[I,I],SqrtScaleA));
            Dec(I);
        end;
    end
    else
    begin
        
        //
        // Solve L*y=b first
        //
        I:=0;
        while I<=N-1 do
        begin
            if I>0 then
            begin
                for i_ := 0 to I-1 do
                begin
                    Tmp[i_] := C_MulR(CHA[I,i_],SqrtScaleA);
                end;
                V := C_Complex(0.0);
                for i_ := 0 to I-1 do
                begin
                    V := C_Add(V,C_Mul(Tmp[i_],XB[i_]));
                end;
                XB[I] := C_Sub(XB[I],V);
            end;
            XB[I] := C_Div(XB[I],C_MulR(CHA[I,I],SqrtScaleA));
            Inc(I);
        end;
        
        //
        // Solve L'*x=y then.
        //
        I:=N-1;
        while I>=0 do
        begin
            XB[I] := C_Div(XB[I],C_MulR(Conj(CHA[I,I]),SqrtScaleA));
            if I>0 then
            begin
                V := XB[I];
                for i_ := 0 to I-1 do
                begin
                    Tmp[i_] := C_MulR(Conj(CHA[I,i_]),SqrtScaleA);
                end;
                for i_ := 0 to I-1 do
                begin
                    XB[i_] := C_Sub(XB[i_], C_Mul(V, Tmp[i_]));
                end;
            end;
            Dec(I);
        end;
    end;
end;


end.