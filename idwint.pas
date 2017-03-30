{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Copyright (c) 2007-2010, Sergey Bochkanov (ALGLIB project).

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
unit idwint;
interface
uses Math, Sysutils, Ap, tsort, nearestneighbor, reflections, hblas, creflections, sblas, ablasf, ablas, ortfac, blas, rotations, bdsvd, svd, hqrnd, matgen, trfac, trlinsolve, safesolve, rcond, xblas, densesolver;

type
(*************************************************************************
IDW interpolant.
*************************************************************************)
IDWInterpolant = record
    N : AlglibInteger;
    NX : AlglibInteger;
    D : AlglibInteger;
    R : Double;
    NW : AlglibInteger;
    Tree : KDTree;
    ModelType : AlglibInteger;
    Q : TReal2DArray;
    XBuf : TReal1DArray;
    TBuf : TInteger1DArray;
    RBuf : TReal1DArray;
    XYBuf : TReal2DArray;
    DebugSolverFailures : AlglibInteger;
    DebugWorstRCond : Double;
    DebugBestRCond : Double;
end;



function IDWCalc(var Z : IDWInterpolant; const X : TReal1DArray):Double;
procedure IDWBuildModifiedShepard(const XY : TReal2DArray;
     N : AlglibInteger;
     NX : AlglibInteger;
     D : AlglibInteger;
     NQ : AlglibInteger;
     NW : AlglibInteger;
     var Z : IDWInterpolant);
procedure IDWBuildModifiedShepardR(const XY : TReal2DArray;
     N : AlglibInteger;
     NX : AlglibInteger;
     R : Double;
     var Z : IDWInterpolant);
procedure IDWBuildNoisy(const XY : TReal2DArray;
     N : AlglibInteger;
     NX : AlglibInteger;
     D : AlglibInteger;
     NQ : AlglibInteger;
     NW : AlglibInteger;
     var Z : IDWInterpolant);

implementation

const
    IDWQFactor = Double(1.5);
    IDWKMin = 5;

function IDWCalcQ(var Z : IDWInterpolant;
     const X : TReal1DArray;
     K : AlglibInteger):Double;forward;
procedure IDWInit1(N : AlglibInteger;
     NX : AlglibInteger;
     D : AlglibInteger;
     NQ : AlglibInteger;
     NW : AlglibInteger;
     var Z : IDWInterpolant);forward;
procedure IDWInternalSolver(var Y : TReal1DArray;
     var W : TReal1DArray;
     var FMatrix : TReal2DArray;
     var Temp : TReal1DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var X : TReal1DArray;
     var TaskRCond : Double);forward;


(*************************************************************************
IDW interpolation

INPUT PARAMETERS:
    Z   -   IDW interpolant built with one of model building
            subroutines.
    X   -   array[0..NX-1], interpolation point

Result:
    IDW interpolant Z(X)

  -- ALGLIB --
     Copyright 02.03.2010 by Bochkanov Sergey
*************************************************************************)
function IDWCalc(var Z : IDWInterpolant; const X : TReal1DArray):Double;
var
    NX : AlglibInteger;
    I : AlglibInteger;
    K : AlglibInteger;
    R : Double;
    S : Double;
    W : Double;
    V1 : Double;
    V2 : Double;
    D0 : Double;
    DI : Double;
    V : Double;
begin
    if Z.ModelType=0 then
    begin
        
        //
        // NQ/NW-based model
        //
        NX := Z.NX;
        KDTreeQueryKNN(Z.Tree, X, Z.NW, True);
        KDTreeQueryResultsDistances(Z.Tree, Z.RBuf, K);
        KDTreeQueryResultsTags(Z.Tree, Z.TBuf, K);
    end;
    if Z.ModelType=1 then
    begin
        
        //
        // R-based model
        //
        NX := Z.NX;
        KDTreeQueryRNN(Z.Tree, X, Z.R, True);
        KDTreeQueryResultsDistances(Z.Tree, Z.RBuf, K);
        KDTreeQueryResultsTags(Z.Tree, Z.TBuf, K);
        if K<IDWKMin then
        begin
            
            //
            // we need at least IDWKMin points
            //
            KDTreeQueryKNN(Z.Tree, X, IDWKMin, True);
            KDTreeQueryResultsDistances(Z.Tree, Z.RBuf, K);
            KDTreeQueryResultsTags(Z.Tree, Z.TBuf, K);
        end;
    end;
    
    //
    // initialize weights for linear/quadratic members calculation.
    //
    // NOTE 1: weights are calculated using NORMALIZED modified
    // Shepard's formula. Original formula gives w(i) = sqr((R-di)/(R*di)),
    // where di is i-th distance, R is max(di). Modified formula have
    // following form:
    //     w_mod(i) = 1, if di=d0
    //     w_mod(i) = w(i)/w(0), if di<>d0
    //
    // NOTE 2: self-match is USED for this query
    //
    // NOTE 3: last point almost always gain zero weight, but it MUST
    // be used for fitting because sometimes it will gain NON-ZERO
    // weight - for example, when all distances are equal.
    //
    R := Z.RBuf[K-1];
    D0 := Z.RBuf[0];
    Result := 0;
    S := 0;
    I:=0;
    while I<=K-1 do
    begin
        DI := Z.RBuf[I];
        if AP_FP_Eq(DI,D0) then
        begin
            
            //
            // distance is equal to shortest, set it 1.0
            // without explicitly calculating (which would give
            // us same result, but 'll expose us to the risk of
            // division by zero).
            //
            W := 1;
        end
        else
        begin
            
            //
            // use normalized formula
            //
            V1 := (R-DI)/(R-D0);
            V2 := D0/DI;
            W := AP_Sqr(V1*V2);
        end;
        Result := Result+W*IDWCalcQ(Z, X, Z.TBuf[I]);
        S := S+W;
        Inc(I);
    end;
    Result := Result/S;
end;


(*************************************************************************
IDW interpolant using modified Shepard method for uniform point
distributions.

INPUT PARAMETERS:
    XY  -   X and Y values, array[0..N-1,0..NX].
            First NX columns contain X-values, last column contain
            Y-values.
    N   -   number of nodes, N>0.
    NX  -   space dimension, NX>=1.
    D   -   nodal function type, either:
            * 0     constant  model.  Just  for  demonstration only, worst
                    model ever.
            * 1     linear model, least squares fitting. Simpe  model  for
                    datasets too small for quadratic models
            * 2     quadratic  model,  least  squares  fitting. Best model
                    available (if your dataset is large enough).
            * -1    "fast"  linear  model,  use  with  caution!!!   It  is
                    significantly  faster than linear/quadratic and better
                    than constant model. But it is less robust (especially
                    in the presence of noise).
    NQ  -   number of points used to calculate  nodal  functions  (ignored
            for constant models). NQ should be LARGER than:
            * max(1.5*(1+NX),2^NX+1) for linear model,
            * max(3/4*(NX+2)*(NX+1),2^NX+1) for quadratic model.
            Values less than this threshold will be silently increased.
    NW  -   number of points used to calculate weights and to interpolate.
            Required: >=2^NX+1, values less than this  threshold  will  be
            silently increased.
            Recommended value: about 2*NQ

OUTPUT PARAMETERS:
    Z   -   IDW interpolant.
    
NOTES:
  * best results are obtained with quadratic models, worst - with constant
    models
  * when N is large, NQ and NW must be significantly smaller than  N  both
    to obtain optimal performance and to obtain optimal accuracy. In 2  or
    3-dimensional tasks NQ=15 and NW=25 are good values to start with.
  * NQ  and  NW  may  be  greater  than  N.  In  such  cases  they will be
    automatically decreased.
  * this subroutine is always succeeds (as long as correct parameters  are
    passed).
  * see  'Multivariate  Interpolation  of Large Sets of Scattered Data' by
    Robert J. Renka for more information on this algorithm.
  * this subroutine assumes that point distribution is uniform at the small
    scales.  If  it  isn't  -  for  example,  points are concentrated along
    "lines", but "lines" distribution is uniform at the larger scale - then
    you should use IDWBuildModifiedShepardR()


  -- ALGLIB PROJECT --
     Copyright 02.03.2010 by Bochkanov Sergey
*************************************************************************)
procedure IDWBuildModifiedShepard(const XY : TReal2DArray;
     N : AlglibInteger;
     NX : AlglibInteger;
     D : AlglibInteger;
     NQ : AlglibInteger;
     NW : AlglibInteger;
     var Z : IDWInterpolant);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    J2 : AlglibInteger;
    J3 : AlglibInteger;
    V : Double;
    R : Double;
    S : Double;
    D0 : Double;
    DI : Double;
    V1 : Double;
    V2 : Double;
    NC : AlglibInteger;
    Offs : AlglibInteger;
    X : TReal1DArray;
    QRBuf : TReal1DArray;
    QXYBuf : TReal2DArray;
    Y : TReal1DArray;
    FMatrix : TReal2DArray;
    W : TReal1DArray;
    QSol : TReal1DArray;
    Temp : TReal1DArray;
    Tags : TInteger1DArray;
    Info : AlglibInteger;
    TaskRCond : Double;
begin
    
    //
    // assertions
    //
    Assert(N>0, 'IDWBuildModifiedShepard: N<=0!');
    Assert(NX>=1, 'IDWBuildModifiedShepard: NX<1!');
    Assert((D>=-1) and (D<=2), 'IDWBuildModifiedShepard: D<>-1 and D<>0 and D<>1 and D<>2!');
    
    //
    // Correct parameters if needed
    //
    if D=1 then
    begin
        NQ := Max(NQ, Ceil(IDWQFactor*(1+NX))+1);
        NQ := Max(NQ, Round(Power(2, NX))+1);
    end;
    if D=2 then
    begin
        NQ := Max(NQ, Ceil(IDWQFactor*(NX+2)*(NX+1)/2)+1);
        NQ := Max(NQ, Round(Power(2, NX))+1);
    end;
    NW := Max(NW, Round(Power(2, NX))+1);
    NQ := Min(NQ, N);
    NW := Min(NW, N);
    
    //
    // primary initialization of Z
    //
    IDWInit1(N, NX, D, NQ, NW, Z);
    Z.ModelType := 0;
    
    //
    // Create KD-tree
    //
    SetLength(Tags, N);
    I:=0;
    while I<=N-1 do
    begin
        Tags[I] := I;
        Inc(I);
    end;
    KDTreeBuildTagged(XY, Tags, N, NX, 1, 2, Z.Tree);
    
    //
    // build nodal functions
    //
    SetLength(Temp, NQ+1);
    SetLength(X, NX);
    SetLength(QRBuf, NQ);
    SetLength(QXYBuf, NQ, NX+1);
    if D=-1 then
    begin
        SetLength(W, NQ);
    end;
    if D=1 then
    begin
        SetLength(Y, NQ);
        SetLength(W, NQ);
        SetLength(QSol, NX);
        
        //
        // NX for linear members,
        // 1 for temporary storage
        //
        SetLength(FMatrix, NQ, NX+1);
    end;
    if D=2 then
    begin
        SetLength(Y, NQ);
        SetLength(W, NQ);
        SetLength(QSol, NX+Round(NX*(NX+1)*Double(0.5)));
        
        //
        // NX for linear members,
        // Round(NX*(NX+1)*0.5) for quadratic model,
        // 1 for temporary storage
        //
        SetLength(FMatrix, NQ, NX+Round(NX*(NX+1)*Double(0.5))+1);
    end;
    I:=0;
    while I<=N-1 do
    begin
        
        //
        // Initialize center and function value.
        // If D=0 it is all what we need
        //
        APVMove(@Z.Q[I][0], 0, NX, @XY[I][0], 0, NX);
        if D=0 then
        begin
            Inc(I);
            Continue;
        end;
        
        //
        // calculate weights for linear/quadratic members calculation.
        //
        // NOTE 1: weights are calculated using NORMALIZED modified
        // Shepard's formula. Original formula is w(i) = sqr((R-di)/(R*di)),
        // where di is i-th distance, R is max(di). Modified formula have
        // following form:
        //     w_mod(i) = 1, if di=d0
        //     w_mod(i) = w(i)/w(0), if di<>d0
        //
        // NOTE 2: self-match is NOT used for this query
        //
        // NOTE 3: last point almost always gain zero weight, but it MUST
        // be used for fitting because sometimes it will gain NON-ZERO
        // weight - for example, when all distances are equal.
        //
        APVMove(@X[0], 0, NX-1, @XY[I][0], 0, NX-1);
        KDTreeQueryKNN(Z.Tree, X, NQ, False);
        KDTreeQueryResultsXY(Z.Tree, QXYBuf, K);
        KDTreeQueryResultsDistances(Z.Tree, QRBuf, K);
        R := QRBuf[K-1];
        D0 := QRBuf[0];
        J:=0;
        while J<=K-1 do
        begin
            DI := QRBuf[J];
            if AP_FP_Eq(DI,D0) then
            begin
                
                //
                // distance is equal to shortest, set it 1.0
                // without explicitly calculating (which would give
                // us same result, but 'll expose us to the risk of
                // division by zero).
                //
                W[J] := 1;
            end
            else
            begin
                
                //
                // use normalized formula
                //
                V1 := (R-DI)/(R-D0);
                V2 := D0/DI;
                W[J] := AP_Sqr(V1*V2);
            end;
            Inc(J);
        end;
        
        //
        // calculate linear/quadratic members
        //
        if D=-1 then
        begin
            
            //
            // "Fast" linear nodal function calculated using
            // inverse distance weighting
            //
            J:=0;
            while J<=NX-1 do
            begin
                X[J] := 0;
                Inc(J);
            end;
            S := 0;
            J:=0;
            while J<=K-1 do
            begin
                
                //
                // calculate J-th inverse distance weighted gradient:
                //     grad_k = (y_j-y_k)*(x_j-x_k)/sqr(norm(x_j-x_k))
                //     grad   = sum(wk*grad_k)/sum(w_k)
                //
                V := 0;
                J2:=0;
                while J2<=NX-1 do
                begin
                    V := V+AP_Sqr(QXYBuf[J,J2]-XY[I,J2]);
                    Inc(J2);
                end;
                
                //
                // Although x_j<>x_k, sqr(norm(x_j-x_k)) may be zero due to
                // underflow. If it is, we assume than J-th gradient is zero
                // (i.e. don't add anything)
                //
                if AP_FP_Neq(V,0) then
                begin
                    J2:=0;
                    while J2<=NX-1 do
                    begin
                        X[J2] := X[J2]+W[J]*(QXYBuf[J,NX]-XY[I,NX])*(QXYBuf[J,J2]-XY[I,J2])/V;
                        Inc(J2);
                    end;
                end;
                S := S+W[J];
                Inc(J);
            end;
            J:=0;
            while J<=NX-1 do
            begin
                Z.Q[I,NX+1+J] := X[J]/S;
                Inc(J);
            end;
        end
        else
        begin
            
            //
            // Least squares models: build
            //
            if D=1 then
            begin
                
                //
                // Linear nodal function calculated using
                // least squares fitting to its neighbors
                //
                J:=0;
                while J<=K-1 do
                begin
                    J2:=0;
                    while J2<=NX-1 do
                    begin
                        FMatrix[J,J2] := QXYBuf[J,J2]-XY[I,J2];
                        Inc(J2);
                    end;
                    Y[J] := QXYBuf[J,NX]-XY[I,NX];
                    Inc(J);
                end;
                NC := NX;
            end;
            if D=2 then
            begin
                
                //
                // Quadratic nodal function calculated using
                // least squares fitting to its neighbors
                //
                J:=0;
                while J<=K-1 do
                begin
                    Offs := 0;
                    J2:=0;
                    while J2<=NX-1 do
                    begin
                        FMatrix[J,Offs] := QXYBuf[J,J2]-XY[I,J2];
                        Offs := Offs+1;
                        Inc(J2);
                    end;
                    J2:=0;
                    while J2<=NX-1 do
                    begin
                        J3:=J2;
                        while J3<=NX-1 do
                        begin
                            FMatrix[J,Offs] := (QXYBuf[J,J2]-XY[I,J2])*(QXYBuf[J,J3]-XY[I,J3]);
                            Offs := Offs+1;
                            Inc(J3);
                        end;
                        Inc(J2);
                    end;
                    Y[J] := QXYBuf[J,NX]-XY[I,NX];
                    Inc(J);
                end;
                NC := NX+Round(NX*(NX+1)*Double(0.5));
            end;
            IDWInternalSolver(Y, W, FMatrix, Temp, K, NC, Info, QSol, TaskRCond);
            
            //
            // Least squares models: copy results
            //
            if Info>0 then
            begin
                
                //
                // LLS task is solved, copy results
                //
                Z.DebugWorstRCond := Min(Z.DebugWorstRCond, TaskRCond);
                Z.DebugBestRCond := Max(Z.DebugBestRCond, TaskRCond);
                J:=0;
                while J<=NC-1 do
                begin
                    Z.Q[I,NX+1+J] := QSol[J];
                    Inc(J);
                end;
            end
            else
            begin
                
                //
                // Solver failure, very strange, but we will use
                // zero values to handle it.
                //
                Z.DebugSolverFailures := Z.DebugSolverFailures+1;
                J:=0;
                while J<=NC-1 do
                begin
                    Z.Q[I,NX+1+J] := 0;
                    Inc(J);
                end;
            end;
        end;
        Inc(I);
    end;
end;


(*************************************************************************
IDW interpolant using modified Shepard method for non-uniform datasets.

This type of model uses  constant  nodal  functions and interpolates using
all nodes which are closer than user-specified radius R. It  may  be  used
when points distribution is non-uniform at the small scale, but it  is  at
the distances as large as R.

INPUT PARAMETERS:
    XY  -   X and Y values, array[0..N-1,0..NX].
            First NX columns contain X-values, last column contain
            Y-values.
    N   -   number of nodes, N>0.
    NX  -   space dimension, NX>=1.
    R   -   radius, R>0

OUTPUT PARAMETERS:
    Z   -   IDW interpolant.

NOTES:
* if there is less than IDWKMin points within  R-ball,  algorithm  selects
  IDWKMin closest ones, so that continuity properties of  interpolant  are
  preserved even far from points.

  -- ALGLIB PROJECT --
     Copyright 11.04.2010 by Bochkanov Sergey
*************************************************************************)
procedure IDWBuildModifiedShepardR(const XY : TReal2DArray;
     N : AlglibInteger;
     NX : AlglibInteger;
     R : Double;
     var Z : IDWInterpolant);
var
    I : AlglibInteger;
    Tags : TInteger1DArray;
begin
    
    //
    // assertions
    //
    Assert(N>0, 'IDWBuildModifiedShepardR: N<=0!');
    Assert(NX>=1, 'IDWBuildModifiedShepardR: NX<1!');
    Assert(AP_FP_Greater(R,0), 'IDWBuildModifiedShepardR: R<=0!');
    
    //
    // primary initialization of Z
    //
    IDWInit1(N, NX, 0, 0, N, Z);
    Z.ModelType := 1;
    Z.R := R;
    
    //
    // Create KD-tree
    //
    SetLength(Tags, N);
    I:=0;
    while I<=N-1 do
    begin
        Tags[I] := I;
        Inc(I);
    end;
    KDTreeBuildTagged(XY, Tags, N, NX, 1, 2, Z.Tree);
    
    //
    // build nodal functions
    //
    I:=0;
    while I<=N-1 do
    begin
        APVMove(@Z.Q[I][0], 0, NX, @XY[I][0], 0, NX);
        Inc(I);
    end;
end;


(*************************************************************************
IDW model for noisy data.

This subroutine may be used to handle noisy data, i.e. data with noise  in
OUTPUT values.  It differs from IDWBuildModifiedShepard() in the following
aspects:
* nodal functions are not constrained to pass through  nodes:  Qi(xi)<>yi,
  i.e. we have fitting  instead  of  interpolation.
* weights which are used during least  squares fitting stage are all equal
  to 1.0 (independently of distance)
* "fast"-linear or constant nodal functions are not supported (either  not
  robust enough or too rigid)

This problem require far more complex tuning than interpolation  problems.
Below you can find some recommendations regarding this problem:
* focus on tuning NQ; it controls noise reduction. As for NW, you can just
  make it equal to 2*NQ.
* you can use cross-validation to determine optimal NQ.
* optimal NQ is a result of complex tradeoff  between  noise  level  (more
  noise = larger NQ required) and underlying  function  complexity  (given
  fixed N, larger NQ means smoothing of compex features in the data).  For
  example, NQ=N will reduce noise to the minimum level possible,  but  you
  will end up with just constant/linear/quadratic (depending on  D)  least
  squares model for the whole dataset.

INPUT PARAMETERS:
    XY  -   X and Y values, array[0..N-1,0..NX].
            First NX columns contain X-values, last column contain
            Y-values.
    N   -   number of nodes, N>0.
    NX  -   space dimension, NX>=1.
    D   -   nodal function degree, either:
            * 1     linear model, least squares fitting. Simpe  model  for
                    datasets too small for quadratic models (or  for  very
                    noisy problems).
            * 2     quadratic  model,  least  squares  fitting. Best model
                    available (if your dataset is large enough).
    NQ  -   number of points used to calculate nodal functions.  NQ should
            be  significantly   larger   than  1.5  times  the  number  of
            coefficients in a nodal function to overcome effects of noise:
            * larger than 1.5*(1+NX) for linear model,
            * larger than 3/4*(NX+2)*(NX+1) for quadratic model.
            Values less than this threshold will be silently increased.
    NW  -   number of points used to calculate weights and to interpolate.
            Required: >=2^NX+1, values less than this  threshold  will  be
            silently increased.
            Recommended value: about 2*NQ or larger

OUTPUT PARAMETERS:
    Z   -   IDW interpolant.

NOTES:
  * best results are obtained with quadratic models, linear models are not
    recommended to use unless you are pretty sure that it is what you want
  * this subroutine is always succeeds (as long as correct parameters  are
    passed).
  * see  'Multivariate  Interpolation  of Large Sets of Scattered Data' by
    Robert J. Renka for more information on this algorithm.


  -- ALGLIB PROJECT --
     Copyright 02.03.2010 by Bochkanov Sergey
*************************************************************************)
procedure IDWBuildNoisy(const XY : TReal2DArray;
     N : AlglibInteger;
     NX : AlglibInteger;
     D : AlglibInteger;
     NQ : AlglibInteger;
     NW : AlglibInteger;
     var Z : IDWInterpolant);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    J2 : AlglibInteger;
    J3 : AlglibInteger;
    V : Double;
    NC : AlglibInteger;
    Offs : AlglibInteger;
    TaskRCond : Double;
    X : TReal1DArray;
    QRBuf : TReal1DArray;
    QXYBuf : TReal2DArray;
    Y : TReal1DArray;
    W : TReal1DArray;
    FMatrix : TReal2DArray;
    QSol : TReal1DArray;
    Tags : TInteger1DArray;
    Temp : TReal1DArray;
    Info : AlglibInteger;
begin
    
    //
    // assertions
    //
    Assert(N>0, 'IDWBuildNoisy: N<=0!');
    Assert(NX>=1, 'IDWBuildNoisy: NX<1!');
    Assert((D>=1) and (D<=2), 'IDWBuildNoisy: D<>1 and D<>2!');
    
    //
    // Correct parameters if needed
    //
    if D=1 then
    begin
        NQ := Max(NQ, Ceil(IDWQFactor*(1+NX))+1);
    end;
    if D=2 then
    begin
        NQ := Max(NQ, Ceil(IDWQFactor*(NX+2)*(NX+1)/2)+1);
    end;
    NW := Max(NW, Round(Power(2, NX))+1);
    NQ := Min(NQ, N);
    NW := Min(NW, N);
    
    //
    // primary initialization of Z
    //
    IDWInit1(N, NX, D, NQ, NW, Z);
    Z.ModelType := 0;
    
    //
    // Create KD-tree
    //
    SetLength(Tags, N);
    I:=0;
    while I<=N-1 do
    begin
        Tags[I] := I;
        Inc(I);
    end;
    KDTreeBuildTagged(XY, Tags, N, NX, 1, 2, Z.Tree);
    
    //
    // build nodal functions
    // (special algorithm for noisy data is used)
    //
    SetLength(Temp, NQ+1);
    SetLength(X, NX);
    SetLength(QRBuf, NQ);
    SetLength(QXYBuf, NQ, NX+1);
    if D=1 then
    begin
        SetLength(Y, NQ);
        SetLength(W, NQ);
        SetLength(QSol, 1+NX);
        
        //
        // 1 for constant member,
        // NX for linear members,
        // 1 for temporary storage
        //
        SetLength(FMatrix, NQ, 1+NX+1);
    end;
    if D=2 then
    begin
        SetLength(Y, NQ);
        SetLength(W, NQ);
        SetLength(QSol, 1+NX+Round(NX*(NX+1)*Double(0.5)));
        
        //
        // 1 for constant member,
        // NX for linear members,
        // Round(NX*(NX+1)*0.5) for quadratic model,
        // 1 for temporary storage
        //
        SetLength(FMatrix, NQ, 1+NX+Round(NX*(NX+1)*Double(0.5))+1);
    end;
    I:=0;
    while I<=N-1 do
    begin
        
        //
        // Initialize center.
        //
        APVMove(@Z.Q[I][0], 0, NX-1, @XY[I][0], 0, NX-1);
        
        //
        // Calculate linear/quadratic members
        // using least squares fit
        // NOTE 1: all weight are equal to 1.0
        // NOTE 2: self-match is USED for this query
        //
        APVMove(@X[0], 0, NX-1, @XY[I][0], 0, NX-1);
        KDTreeQueryKNN(Z.Tree, X, NQ, True);
        KDTreeQueryResultsXY(Z.Tree, QXYBuf, K);
        KDTreeQueryResultsDistances(Z.Tree, QRBuf, K);
        if D=1 then
        begin
            
            //
            // Linear nodal function calculated using
            // least squares fitting to its neighbors
            //
            J:=0;
            while J<=K-1 do
            begin
                FMatrix[J,0] := Double(1.0);
                J2:=0;
                while J2<=NX-1 do
                begin
                    FMatrix[J,1+J2] := QXYBuf[J,J2]-XY[I,J2];
                    Inc(J2);
                end;
                Y[J] := QXYBuf[J,NX];
                W[J] := 1;
                Inc(J);
            end;
            NC := 1+NX;
        end;
        if D=2 then
        begin
            
            //
            // Quadratic nodal function calculated using
            // least squares fitting to its neighbors
            //
            J:=0;
            while J<=K-1 do
            begin
                FMatrix[J,0] := 1;
                Offs := 1;
                J2:=0;
                while J2<=NX-1 do
                begin
                    FMatrix[J,Offs] := QXYBuf[J,J2]-XY[I,J2];
                    Offs := Offs+1;
                    Inc(J2);
                end;
                J2:=0;
                while J2<=NX-1 do
                begin
                    J3:=J2;
                    while J3<=NX-1 do
                    begin
                        FMatrix[J,Offs] := (QXYBuf[J,J2]-XY[I,J2])*(QXYBuf[J,J3]-XY[I,J3]);
                        Offs := Offs+1;
                        Inc(J3);
                    end;
                    Inc(J2);
                end;
                Y[J] := QXYBuf[J,NX];
                W[J] := 1;
                Inc(J);
            end;
            NC := 1+NX+Round(NX*(NX+1)*Double(0.5));
        end;
        IDWInternalSolver(Y, W, FMatrix, Temp, K, NC, Info, QSol, TaskRCond);
        
        //
        // Least squares models: copy results
        //
        if Info>0 then
        begin
            
            //
            // LLS task is solved, copy results
            //
            Z.DebugWorstRCond := Min(Z.DebugWorstRCond, TaskRCond);
            Z.DebugBestRCond := Max(Z.DebugBestRCond, TaskRCond);
            J:=0;
            while J<=NC-1 do
            begin
                Z.Q[I,NX+J] := QSol[J];
                Inc(J);
            end;
        end
        else
        begin
            
            //
            // Solver failure, very strange, but we will use
            // zero values to handle it.
            //
            Z.DebugSolverFailures := Z.DebugSolverFailures+1;
            V := 0;
            J:=0;
            while J<=K-1 do
            begin
                V := V+QXYBuf[J,NX];
                Inc(J);
            end;
            Z.Q[I,NX] := V/K;
            J:=0;
            while J<=NC-2 do
            begin
                Z.Q[I,NX+1+J] := 0;
                Inc(J);
            end;
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Internal subroutine: K-th nodal function calculation

  -- ALGLIB --
     Copyright 02.03.2010 by Bochkanov Sergey
*************************************************************************)
function IDWCalcQ(var Z : IDWInterpolant;
     const X : TReal1DArray;
     K : AlglibInteger):Double;
var
    NX : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    Offs : AlglibInteger;
begin
    NX := Z.NX;
    
    //
    // constant member
    //
    Result := Z.Q[K,NX];
    
    //
    // linear members
    //
    if Z.D>=1 then
    begin
        I:=0;
        while I<=NX-1 do
        begin
            Result := Result+Z.Q[K,NX+1+I]*(X[I]-Z.Q[K,I]);
            Inc(I);
        end;
    end;
    
    //
    // quadratic members
    //
    if Z.D>=2 then
    begin
        Offs := NX+1+NX;
        I:=0;
        while I<=NX-1 do
        begin
            J:=I;
            while J<=NX-1 do
            begin
                Result := Result+Z.Q[K,Offs]*(X[I]-Z.Q[K,I])*(X[J]-Z.Q[K,J]);
                Offs := Offs+1;
                Inc(J);
            end;
            Inc(I);
        end;
    end;
end;


(*************************************************************************
Initialization of internal structures.

It assumes correctness of all parameters.

  -- ALGLIB --
     Copyright 02.03.2010 by Bochkanov Sergey
*************************************************************************)
procedure IDWInit1(N : AlglibInteger;
     NX : AlglibInteger;
     D : AlglibInteger;
     NQ : AlglibInteger;
     NW : AlglibInteger;
     var Z : IDWInterpolant);
begin
    Z.DebugSolverFailures := 0;
    Z.DebugWorstRCond := Double(1.0);
    Z.DebugBestRCond := 0;
    Z.N := N;
    Z.NX := NX;
    Z.D := 0;
    if D=1 then
    begin
        Z.D := 1;
    end;
    if D=2 then
    begin
        Z.D := 2;
    end;
    if D=-1 then
    begin
        Z.D := 1;
    end;
    Z.NW := NW;
    if D=-1 then
    begin
        SetLength(Z.Q, N, NX+1+NX);
    end;
    if D=0 then
    begin
        SetLength(Z.Q, N, NX+1);
    end;
    if D=1 then
    begin
        SetLength(Z.Q, N, NX+1+NX);
    end;
    if D=2 then
    begin
        SetLength(Z.Q, N, NX+1+NX+Round(NX*(NX+1)*Double(0.5)));
    end;
    SetLength(Z.TBuf, NW);
    SetLength(Z.RBuf, NW);
    SetLength(Z.XYBuf, NW, NX+1);
    SetLength(Z.XBuf, NX);
end;


(*************************************************************************
Linear least squares solver for small tasks.

Works faster than standard ALGLIB solver in non-degenerate cases  (due  to
absense of internal allocations and optimized row/colums).  In  degenerate
cases it calls standard solver, which results in small performance penalty
associated with preliminary steps.

INPUT PARAMETERS:
    Y           array[0..N-1]
    W           array[0..N-1]
    FMatrix     array[0..N-1,0..M], have additional column for temporary
                values
    Temp        array[0..N]
*************************************************************************)
procedure IDWInternalSolver(var Y : TReal1DArray;
     var W : TReal1DArray;
     var FMatrix : TReal2DArray;
     var Temp : TReal1DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var X : TReal1DArray;
     var TaskRCond : Double);
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    Tau : Double;
    B : TReal1DArray;
    SRep : DenseSolverLSReport;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    
    //
    // set up info
    //
    Info := 1;
    
    //
    // prepare matrix
    //
    I:=0;
    while I<=N-1 do
    begin
        FMatrix[I,M] := Y[I];
        V := W[I];
        APVMul(@FMatrix[I][0], 0, M, V);
        Inc(I);
    end;
    
    //
    // use either fast algorithm or general algorithm
    //
    if M<=N then
    begin
        
        //
        // QR decomposition
        // We assume that M<=N (we would have called LSFit() otherwise)
        //
        I:=0;
        while I<=M-1 do
        begin
            if I<N-1 then
            begin
                i1_ := (I) - (1);
                for i_ := 1 to N-I do
                begin
                    Temp[i_] := FMatrix[i_+i1_,I];
                end;
                GenerateReflection(Temp, N-I, Tau);
                FMatrix[I,I] := Temp[1];
                Temp[1] := 1;
                J:=I+1;
                while J<=M do
                begin
                    i1_ := (1)-(I);
                    V := 0.0;
                    for i_ := I to N-1 do
                    begin
                        V := V + FMatrix[i_,J]*Temp[i_+i1_];
                    end;
                    V := Tau*V;
                    i1_ := (1) - (I);
                    for i_ := I to N-1 do
                    begin
                        FMatrix[i_,J] := FMatrix[i_,J] - V*Temp[i_+i1_];
                    end;
                    Inc(J);
                end;
            end;
            Inc(I);
        end;
        
        //
        // Check condition number
        //
        TaskRCond := RMatrixTRRCondInf(FMatrix, M, True, False);
        
        //
        // use either fast algorithm for non-degenerate cases
        // or slow algorithm for degenerate cases
        //
        if AP_FP_Greater(TaskRCond,10000*N*MachineEpsilon) then
        begin
            
            //
            // solve triangular system R*x = FMatrix[0:M-1,M]
            // using fast algorithm, then exit
            //
            X[M-1] := FMatrix[M-1,M]/FMatrix[M-1,M-1];
            I:=M-2;
            while I>=0 do
            begin
                V := APVDotProduct(@FMatrix[I][0], I+1, M-1, @X[0], I+1, M-1);
                X[I] := (FMatrix[I,M]-V)/FMatrix[I,I];
                Dec(I);
            end;
        end
        else
        begin
            
            //
            // use more general algorithm
            //
            SetLength(B, M);
            I:=0;
            while I<=M-1 do
            begin
                J:=0;
                while J<=I-1 do
                begin
                    FMatrix[I,J] := Double(0.0);
                    Inc(J);
                end;
                B[I] := FMatrix[I,M];
                Inc(I);
            end;
            RMatrixSolveLS(FMatrix, M, M, B, 10000*MachineEpsilon, Info, SRep, X);
        end;
    end
    else
    begin
        
        //
        // use more general algorithm
        //
        SetLength(B, N);
        I:=0;
        while I<=N-1 do
        begin
            B[I] := FMatrix[I,M];
            Inc(I);
        end;
        RMatrixSolveLS(FMatrix, N, M, B, 10000*MachineEpsilon, Info, SRep, X);
        TaskRCond := SRep.R2;
    end;
end;


end.