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
unit nearestneighbor;
interface
uses Math, Sysutils, Ap, tsort;

type
KDTree = record
    N : AlglibInteger;
    NX : AlglibInteger;
    NY : AlglibInteger;
    NormType : AlglibInteger;
    DistMatrixType : AlglibInteger;
    XY : TReal2DArray;
    Tags : TInteger1DArray;
    BoxMin : TReal1DArray;
    BoxMax : TReal1DArray;
    CurBoxMin : TReal1DArray;
    CurBoxMax : TReal1DArray;
    CurDist : Double;
    Nodes : TInteger1DArray;
    Splits : TReal1DArray;
    X : TReal1DArray;
    KNeeded : AlglibInteger;
    RNeeded : Double;
    SelfMatch : Boolean;
    ApproxF : Double;
    KCur : AlglibInteger;
    Idx : TInteger1DArray;
    R : TReal1DArray;
    Buf : TReal1DArray;
    DebugCounter : AlglibInteger;
end;



procedure KDTreeBuild(const XY : TReal2DArray;
     N : AlglibInteger;
     NX : AlglibInteger;
     NY : AlglibInteger;
     NormType : AlglibInteger;
     var KDT : KDTree);
procedure KDTreeBuildTagged(const XY : TReal2DArray;
     const Tags : TInteger1DArray;
     N : AlglibInteger;
     NX : AlglibInteger;
     NY : AlglibInteger;
     NormType : AlglibInteger;
     var KDT : KDTree);
function KDTreeQueryKNN(var KDT : KDTree;
     const X : TReal1DArray;
     K : AlglibInteger;
     SelfMatch : Boolean):AlglibInteger;
function KDTreeQueryRNN(var KDT : KDTree;
     const X : TReal1DArray;
     R : Double;
     SelfMatch : Boolean):AlglibInteger;
function KDTreeQueryAKNN(var KDT : KDTree;
     const X : TReal1DArray;
     K : AlglibInteger;
     SelfMatch : Boolean;
     Eps : Double):AlglibInteger;
procedure KDTreeQueryResultsX(const KDT : KDTree;
     var X : TReal2DArray;
     var K : AlglibInteger);
procedure KDTreeQueryResultsXY(const KDT : KDTree;
     var XY : TReal2DArray;
     var K : AlglibInteger);
procedure KDTreeQueryResultsTags(const KDT : KDTree;
     var Tags : TInteger1DArray;
     var K : AlglibInteger);
procedure KDTreeQueryResultsDistances(const KDT : KDTree;
     var R : TReal1DArray;
     var K : AlglibInteger);

implementation

const
    SplitNodeSize = 6;

procedure KDTreeSplit(var KDT : KDTree;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     D : AlglibInteger;
     S : Double;
     var I3 : AlglibInteger);forward;
procedure KDTreeGenerateTreeRec(var KDT : KDTree;
     var NodesOffs : AlglibInteger;
     var SplitsOffs : AlglibInteger;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     MaxLeafSize : AlglibInteger);forward;
procedure KDTreeQueryNNRec(var KDT : KDTree; Offs : AlglibInteger);forward;
procedure KDTreeInitBox(var KDT : KDTree; const X : TReal1DArray);forward;
function VRootFreeNorm(const X : TReal1DArray;
     N : AlglibInteger;
     NormType : AlglibInteger):Double;forward;
function VRootFreeComponentNorm(X : Double;
     NormType : AlglibInteger):Double;forward;
function VRangeDist(X : Double; A : Double; B : Double):Double;forward;


(*************************************************************************
KD-tree creation

This subroutine creates KD-tree from set of X-values and optional Y-values

INPUT PARAMETERS
    XY      -   dataset, array[0..N-1,0..NX+NY-1].
                one row corresponds to one point.
                first NX columns contain X-values, next NY (NY may be zero)
                columns may contain associated Y-values
    N       -   number of points, N>=1
    NX      -   space dimension, NX>=1.
    NY      -   number of optional Y-values, NY>=0.
    NormType-   norm type:
                * 0 denotes infinity-norm
                * 1 denotes 1-norm
                * 2 denotes 2-norm (Euclidean norm)
                
OUTPUT PARAMETERS
    KDT     -   KD-tree
    
    
NOTES

1. KD-tree  creation  have O(N*logN) complexity and O(N*(2*NX+NY))  memory
   requirements.
2. Although KD-trees may be used with any combination of N  and  NX,  they
   are more efficient than brute-force search only when N >> 4^NX. So they
   are most useful in low-dimensional tasks (NX=2, NX=3). NX=1  is another
   inefficient case, because  simple  binary  search  (without  additional
   structures) is much more efficient in such tasks than KD-trees.

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************)
procedure KDTreeBuild(const XY : TReal2DArray;
     N : AlglibInteger;
     NX : AlglibInteger;
     NY : AlglibInteger;
     NormType : AlglibInteger;
     var KDT : KDTree);
var
    Tags : TInteger1DArray;
    I : AlglibInteger;
begin
    Assert(N>=1, 'KDTreeBuild: N<1!');
    Assert(NX>=1, 'KDTreeBuild: NX<1!');
    Assert(NY>=0, 'KDTreeBuild: NY<0!');
    Assert((NormType>=0) and (NormType<=2), 'KDTreeBuild: incorrect NormType!');
    SetLength(Tags, N);
    I:=0;
    while I<=N-1 do
    begin
        Tags[I] := 0;
        Inc(I);
    end;
    KDTreeBuildTagged(XY, Tags, N, NX, NY, NormType, KDT);
end;


(*************************************************************************
KD-tree creation

This  subroutine  creates  KD-tree  from set of X-values, integer tags and
optional Y-values

INPUT PARAMETERS
    XY      -   dataset, array[0..N-1,0..NX+NY-1].
                one row corresponds to one point.
                first NX columns contain X-values, next NY (NY may be zero)
                columns may contain associated Y-values
    Tags    -   tags, array[0..N-1], contains integer tags associated
                with points.
    N       -   number of points, N>=1
    NX      -   space dimension, NX>=1.
    NY      -   number of optional Y-values, NY>=0.
    NormType-   norm type:
                * 0 denotes infinity-norm
                * 1 denotes 1-norm
                * 2 denotes 2-norm (Euclidean norm)

OUTPUT PARAMETERS
    KDT     -   KD-tree

NOTES

1. KD-tree  creation  have O(N*logN) complexity and O(N*(2*NX+NY))  memory
   requirements.
2. Although KD-trees may be used with any combination of N  and  NX,  they
   are more efficient than brute-force search only when N >> 4^NX. So they
   are most useful in low-dimensional tasks (NX=2, NX=3). NX=1  is another
   inefficient case, because  simple  binary  search  (without  additional
   structures) is much more efficient in such tasks than KD-trees.

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************)
procedure KDTreeBuildTagged(const XY : TReal2DArray;
     const Tags : TInteger1DArray;
     N : AlglibInteger;
     NX : AlglibInteger;
     NY : AlglibInteger;
     NormType : AlglibInteger;
     var KDT : KDTree);
var
    I : AlglibInteger;
    J : AlglibInteger;
    MaxNodes : AlglibInteger;
    NodesOffs : AlglibInteger;
    SplitsOffs : AlglibInteger;
begin
    Assert(N>=1, 'KDTreeBuildTagged: N<1!');
    Assert(NX>=1, 'KDTreeBuildTagged: NX<1!');
    Assert(NY>=0, 'KDTreeBuildTagged: NY<0!');
    Assert((NormType>=0) and (NormType<=2), 'KDTreeBuildTagged: incorrect NormType!');
    
    //
    // initialize
    //
    KDT.N := N;
    KDT.NX := NX;
    KDT.NY := NY;
    KDT.NormType := NormType;
    KDT.DistMatrixType := 0;
    SetLength(KDT.XY, N, 2*NX+NY);
    SetLength(KDT.Tags, N);
    SetLength(KDT.Idx, N);
    SetLength(KDT.R, N);
    SetLength(KDT.X, NX);
    SetLength(KDT.Buf, Max(N, NX));
    
    //
    // Initial fill
    //
    I:=0;
    while I<=N-1 do
    begin
        APVMove(@KDT.XY[I][0], 0, NX-1, @XY[I][0], 0, NX-1);
        APVMove(@KDT.XY[I][0], NX, 2*NX+NY-1, @XY[I][0], 0, NX+NY-1);
        KDT.Tags[I] := Tags[I];
        Inc(I);
    end;
    
    //
    // Determine bounding box
    //
    SetLength(KDT.BoxMin, NX);
    SetLength(KDT.BoxMax, NX);
    SetLength(KDT.CurBoxMin, NX);
    SetLength(KDT.CurBoxMax, NX);
    APVMove(@KDT.BoxMin[0], 0, NX-1, @KDT.XY[0][0], 0, NX-1);
    APVMove(@KDT.BoxMax[0], 0, NX-1, @KDT.XY[0][0], 0, NX-1);
    I:=1;
    while I<=N-1 do
    begin
        J:=0;
        while J<=NX-1 do
        begin
            KDT.BoxMin[J] := Min(KDT.BoxMin[J], KDT.XY[I,J]);
            KDT.BoxMax[J] := Max(KDT.BoxMax[J], KDT.XY[I,J]);
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // prepare tree structure
    // * MaxNodes=N because we guarantee no trivial splits, i.e.
    //   every split will generate two non-empty boxes
    //
    MaxNodes := N;
    SetLength(KDT.Nodes, SplitNodeSize*2*MaxNodes);
    SetLength(KDT.Splits, 2*MaxNodes);
    NodesOffs := 0;
    SplitsOffs := 0;
    APVMove(@KDT.CurBoxMin[0], 0, NX-1, @KDT.BoxMin[0], 0, NX-1);
    APVMove(@KDT.CurBoxMax[0], 0, NX-1, @KDT.BoxMax[0], 0, NX-1);
    KDTreeGenerateTreeRec(KDT, NodesOffs, SplitsOffs, 0, N, 8);
    
    //
    // Set current query size to 0
    //
    KDT.KCur := 0;
end;


(*************************************************************************
K-NN query: K nearest neighbors

INPUT PARAMETERS
    KDT         -   KD-tree
    X           -   point, array[0..NX-1].
    K           -   number of neighbors to return, K>=1
    SelfMatch   -   whether self-matches are allowed:
                    * if True, nearest neighbor may be the point itself
                      (if it exists in original dataset)
                    * if False, then only points with non-zero distance
                      are returned

RESULT
    number of actual neighbors found (either K or N, if K>N).

This  subroutine  performs  query  and  stores  its result in the internal
structures of the KD-tree. You can use  following  subroutines  to  obtain
these results:
* KDTreeQueryResultsX() to get X-values
* KDTreeQueryResultsXY() to get X- and Y-values
* KDTreeQueryResultsTags() to get tag values
* KDTreeQueryResultsDistances() to get distances

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************)
function KDTreeQueryKNN(var KDT : KDTree;
     const X : TReal1DArray;
     K : AlglibInteger;
     SelfMatch : Boolean):AlglibInteger;
begin
    Result := KDTreeQueryAKNN(KDT, X, K, SelfMatch, Double(0.0));
end;


(*************************************************************************
R-NN query: all points within R-sphere centered at X

INPUT PARAMETERS
    KDT         -   KD-tree
    X           -   point, array[0..NX-1].
    R           -   radius of sphere (in corresponding norm), R>0
    SelfMatch   -   whether self-matches are allowed:
                    * if True, nearest neighbor may be the point itself
                      (if it exists in original dataset)
                    * if False, then only points with non-zero distance
                      are returned

RESULT
    number of neighbors found, >=0

This  subroutine  performs  query  and  stores  its result in the internal
structures of the KD-tree. You can use  following  subroutines  to  obtain
actual results:
* KDTreeQueryResultsX() to get X-values
* KDTreeQueryResultsXY() to get X- and Y-values
* KDTreeQueryResultsTags() to get tag values
* KDTreeQueryResultsDistances() to get distances

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************)
function KDTreeQueryRNN(var KDT : KDTree;
     const X : TReal1DArray;
     R : Double;
     SelfMatch : Boolean):AlglibInteger;
var
    I : AlglibInteger;
    J : AlglibInteger;
    VX : Double;
    VMin : Double;
    VMax : Double;
begin
    Assert(AP_FP_Greater(R,0), 'KDTreeQueryRNN: incorrect R!');
    
    //
    // Prepare parameters
    //
    KDT.KNeeded := 0;
    if KDT.NormType<>2 then
    begin
        KDT.RNeeded := R;
    end
    else
    begin
        KDT.RNeeded := AP_Sqr(R);
    end;
    KDT.SelfMatch := SelfMatch;
    KDT.ApproxF := 1;
    KDT.KCur := 0;
    
    //
    // calculate distance from point to current bounding box
    //
    KDTreeInitBox(KDT, X);
    
    //
    // call recursive search
    // results are returned as heap
    //
    KDTreeQueryNNRec(KDT, 0);
    
    //
    // pop from heap to generate ordered representation
    //
    // last element is non pop'ed because it is already in
    // its place
    //
    Result := KDT.KCur;
    J := KDT.KCur;
    I:=KDT.KCur;
    while I>=2 do
    begin
        TagHeapPopI(KDT.R, KDT.Idx, J);
        Dec(I);
    end;
end;


(*************************************************************************
K-NN query: approximate K nearest neighbors

INPUT PARAMETERS
    KDT         -   KD-tree
    X           -   point, array[0..NX-1].
    K           -   number of neighbors to return, K>=1
    SelfMatch   -   whether self-matches are allowed:
                    * if True, nearest neighbor may be the point itself
                      (if it exists in original dataset)
                    * if False, then only points with non-zero distance
                      are returned
    Eps         -   approximation factor, Eps>=0. eps-approximate  nearest
                    neighbor  is  a  neighbor  whose distance from X is at
                    most (1+eps) times distance of true nearest neighbor.

RESULT
    number of actual neighbors found (either K or N, if K>N).
    
NOTES
    significant performance gain may be achieved only when Eps  is  is  on
    the order of magnitude of 1 or larger.

This  subroutine  performs  query  and  stores  its result in the internal
structures of the KD-tree. You can use  following  subroutines  to  obtain
these results:
* KDTreeQueryResultsX() to get X-values
* KDTreeQueryResultsXY() to get X- and Y-values
* KDTreeQueryResultsTags() to get tag values
* KDTreeQueryResultsDistances() to get distances

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************)
function KDTreeQueryAKNN(var KDT : KDTree;
     const X : TReal1DArray;
     K : AlglibInteger;
     SelfMatch : Boolean;
     Eps : Double):AlglibInteger;
var
    I : AlglibInteger;
    J : AlglibInteger;
    VX : Double;
    VMin : Double;
    VMax : Double;
begin
    Assert(K>0, 'KDTreeQueryKNN: incorrect K!');
    Assert(AP_FP_Greater_Eq(Eps,0), 'KDTreeQueryKNN: incorrect Eps!');
    
    //
    // Prepare parameters
    //
    K := Min(K, KDT.N);
    KDT.KNeeded := K;
    KDT.RNeeded := 0;
    KDT.SelfMatch := SelfMatch;
    if KDT.NormType=2 then
    begin
        KDT.ApproxF := 1/AP_Sqr(1+Eps);
    end
    else
    begin
        KDT.ApproxF := 1/(1+Eps);
    end;
    KDT.KCur := 0;
    
    //
    // calculate distance from point to current bounding box
    //
    KDTreeInitBox(KDT, X);
    
    //
    // call recursive search
    // results are returned as heap
    //
    KDTreeQueryNNRec(KDT, 0);
    
    //
    // pop from heap to generate ordered representation
    //
    // last element is non pop'ed because it is already in
    // its place
    //
    Result := KDT.KCur;
    J := KDT.KCur;
    I:=KDT.KCur;
    while I>=2 do
    begin
        TagHeapPopI(KDT.R, KDT.Idx, J);
        Dec(I);
    end;
end;


(*************************************************************************
X-values from last query

INPUT PARAMETERS
    KDT     -   KD-tree
    X       -   pre-allocated array, at least K rows, at least NX columns
    
OUTPUT PARAMETERS
    X       -   K rows are filled with X-values
    K       -   number of points

NOTE
    points are ordered by distance from the query point (first = closest)

SEE ALSO
* KDTreeQueryResultsXY()            X- and Y-values
* KDTreeQueryResultsTags()          tag values
* KDTreeQueryResultsDistances()     distances

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************)
procedure KDTreeQueryResultsX(const KDT : KDTree;
     var X : TReal2DArray;
     var K : AlglibInteger);
var
    I : AlglibInteger;
begin
    K := KDT.KCur;
    I:=0;
    while I<=K-1 do
    begin
        APVMove(@X[I][0], 0, KDT.NX-1, @KDT.XY[KDT.Idx[I]][0], KDT.NX, 2*KDT.NX-1);
        Inc(I);
    end;
end;


(*************************************************************************
X- and Y-values from last query

INPUT PARAMETERS
    KDT     -   KD-tree
    XY      -   pre-allocated array, at least K rows, at least NX+NY columns

OUTPUT PARAMETERS
    X       -   K rows are filled with points: first NX columns with
                X-values, next NY columns - with Y-values.
    K       -   number of points

NOTE
    points are ordered by distance from the query point (first = closest)

SEE ALSO
* KDTreeQueryResultsX()             X-values
* KDTreeQueryResultsTags()          tag values
* KDTreeQueryResultsDistances()     distances

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************)
procedure KDTreeQueryResultsXY(const KDT : KDTree;
     var XY : TReal2DArray;
     var K : AlglibInteger);
var
    I : AlglibInteger;
begin
    K := KDT.KCur;
    I:=0;
    while I<=K-1 do
    begin
        APVMove(@XY[I][0], 0, KDT.NX+KDT.NY-1, @KDT.XY[KDT.Idx[I]][0], KDT.NX, 2*KDT.NX+KDT.NY-1);
        Inc(I);
    end;
end;


(*************************************************************************
point tags from last query

INPUT PARAMETERS
    KDT     -   KD-tree
    Tags    -   pre-allocated array, at least K elements

OUTPUT PARAMETERS
    Tags    -   first K elements are filled with tags associated with points,
                or, when no tags were supplied, with zeros
    K       -   number of points

NOTE
    points are ordered by distance from the query point (first = closest)

SEE ALSO
* KDTreeQueryResultsX()             X-values
* KDTreeQueryResultsXY()            X- and Y-values
* KDTreeQueryResultsDistances()     distances

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************)
procedure KDTreeQueryResultsTags(const KDT : KDTree;
     var Tags : TInteger1DArray;
     var K : AlglibInteger);
var
    I : AlglibInteger;
begin
    K := KDT.KCur;
    I:=0;
    while I<=K-1 do
    begin
        Tags[I] := KDT.Tags[KDT.Idx[I]];
        Inc(I);
    end;
end;


(*************************************************************************
Distances from last query

INPUT PARAMETERS
    KDT     -   KD-tree
    R       -   pre-allocated array, at least K elements

OUTPUT PARAMETERS
    R       -   first K elements are filled with distances
                (in corresponding norm)
    K       -   number of points

NOTE
    points are ordered by distance from the query point (first = closest)

SEE ALSO
* KDTreeQueryResultsX()             X-values
* KDTreeQueryResultsXY()            X- and Y-values
* KDTreeQueryResultsTags()          tag values

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************)
procedure KDTreeQueryResultsDistances(const KDT : KDTree;
     var R : TReal1DArray;
     var K : AlglibInteger);
var
    I : AlglibInteger;
begin
    K := KDT.KCur;
    
    //
    // unload norms
    //
    // Abs() call is used to handle cases with negative norms
    // (generated during KFN requests)
    //
    if KDT.NormType=0 then
    begin
        I:=0;
        while I<=K-1 do
        begin
            R[I] := AbsReal(KDT.R[I]);
            Inc(I);
        end;
    end;
    if KDT.NormType=1 then
    begin
        I:=0;
        while I<=K-1 do
        begin
            R[I] := AbsReal(KDT.R[I]);
            Inc(I);
        end;
    end;
    if KDT.NormType=2 then
    begin
        I:=0;
        while I<=K-1 do
        begin
            R[I] := Sqrt(AbsReal(KDT.R[I]));
            Inc(I);
        end;
    end;
end;


(*************************************************************************
Rearranges nodes [I1,I2) using partition in D-th dimension with S as threshold.
Returns split position I3: [I1,I3) and [I3,I2) are created as result.

This subroutine doesn't create tree structures, just rearranges nodes.
*************************************************************************)
procedure KDTreeSplit(var KDT : KDTree;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     D : AlglibInteger;
     S : Double;
     var I3 : AlglibInteger);
var
    I : AlglibInteger;
    J : AlglibInteger;
    ILeft : AlglibInteger;
    IRight : AlglibInteger;
    V : Double;
begin
    
    //
    // split XY/Tags in two parts:
    // * [ILeft,IRight] is non-processed part of XY/Tags
    //
    // After cycle is done, we have Ileft=IRight. We deal with
    // this element separately.
    //
    // After this, [I1,ILeft) contains left part, and [ILeft,I2)
    // contains right part.
    //
    ILeft := I1;
    IRight := I2-1;
    while ILeft<IRight do
    begin
        if AP_FP_Less_Eq(KDT.XY[ILeft,D],S) then
        begin
            
            //
            // XY[ILeft] is on its place.
            // Advance ILeft.
            //
            ILeft := ILeft+1;
        end
        else
        begin
            
            //
            // XY[ILeft,..] must be at IRight.
            // Swap and advance IRight.
            //
            I:=0;
            while I<=2*KDT.NX+KDT.NY-1 do
            begin
                V := KDT.XY[ILeft,I];
                KDT.XY[ILeft,I] := KDT.XY[IRight,I];
                KDT.XY[IRight,I] := V;
                Inc(I);
            end;
            J := KDT.Tags[ILeft];
            KDT.Tags[ILeft] := KDT.Tags[IRight];
            KDT.Tags[IRight] := J;
            IRight := IRight-1;
        end;
    end;
    if AP_FP_Less_Eq(KDT.XY[ILeft,D],S) then
    begin
        ILeft := ILeft+1;
    end
    else
    begin
        IRight := IRight-1;
    end;
    I3 := ILeft;
end;


(*************************************************************************
Recursive kd-tree generation subroutine.

PARAMETERS
    KDT         tree
    NodesOffs   unused part of Nodes[] which must be filled by tree
    SplitsOffs  unused part of Splits[]
    I1, I2      points from [I1,I2) are processed
    
NodesOffs[] and SplitsOffs[] must be large enough.

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************)
procedure KDTreeGenerateTreeRec(var KDT : KDTree;
     var NodesOffs : AlglibInteger;
     var SplitsOffs : AlglibInteger;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     MaxLeafSize : AlglibInteger);
var
    N : AlglibInteger;
    NX : AlglibInteger;
    NY : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    OldOffs : AlglibInteger;
    I3 : AlglibInteger;
    CntLess : AlglibInteger;
    CntGreater : AlglibInteger;
    MinV : Double;
    MaxV : Double;
    MinIdx : AlglibInteger;
    MaxIdx : AlglibInteger;
    D : AlglibInteger;
    DS : Double;
    S : Double;
    V : Double;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    Assert(I2>I1, 'KDTreeGenerateTreeRec: internal error');
    
    //
    // Generate leaf if needed
    //
    if I2-I1<=MaxLeafSize then
    begin
        KDT.Nodes[NodesOffs+0] := I2-I1;
        KDT.Nodes[NodesOffs+1] := I1;
        NodesOffs := NodesOffs+2;
        Exit;
    end;
    
    //
    // Load values for easier access
    //
    NX := KDT.NX;
    NY := KDT.NY;
    
    //
    // select dimension to split:
    // * D is a dimension number
    //
    D := 0;
    DS := KDT.CurBoxMax[0]-KDT.CurBoxMin[0];
    I:=1;
    while I<=NX-1 do
    begin
        V := KDT.CurBoxMax[I]-KDT.CurBoxMin[I];
        if AP_FP_Greater(V,DS) then
        begin
            DS := V;
            D := I;
        end;
        Inc(I);
    end;
    
    //
    // Select split position S using sliding midpoint rule,
    // rearrange points into [I1,I3) and [I3,I2)
    //
    S := KDT.CurBoxMin[D]+Double(0.5)*DS;
    i1_ := (I1) - (0);
    for i_ := 0 to I2-I1-1 do
    begin
        KDT.Buf[i_] := KDT.XY[i_+i1_,D];
    end;
    N := I2-I1;
    CntLess := 0;
    CntGreater := 0;
    MinV := KDT.Buf[0];
    MaxV := KDT.Buf[0];
    MinIdx := I1;
    MaxIdx := I1;
    I:=0;
    while I<=N-1 do
    begin
        V := KDT.Buf[I];
        if AP_FP_Less(V,MinV) then
        begin
            MinV := V;
            MinIdx := I1+I;
        end;
        if AP_FP_Greater(V,MaxV) then
        begin
            MaxV := V;
            MaxIdx := I1+I;
        end;
        if AP_FP_Less(V,S) then
        begin
            CntLess := CntLess+1;
        end;
        if AP_FP_Greater(V,S) then
        begin
            CntGreater := CntGreater+1;
        end;
        Inc(I);
    end;
    if (CntLess>0) and (CntGreater>0) then
    begin
        
        //
        // normal midpoint split
        //
        KDTreeSplit(KDT, I1, I2, D, S, I3);
    end
    else
    begin
        
        //
        // sliding midpoint
        //
        if CntLess=0 then
        begin
            
            //
            // 1. move split to MinV,
            // 2. place one point to the left bin (move to I1),
            //    others - to the right bin
            //
            S := MinV;
            if MinIdx<>I1 then
            begin
                I:=0;
                while I<=2*KDT.NX+KDT.NY-1 do
                begin
                    V := KDT.XY[MinIdx,I];
                    KDT.XY[MinIdx,I] := KDT.XY[I1,I];
                    KDT.XY[I1,I] := V;
                    Inc(I);
                end;
                J := KDT.Tags[MinIdx];
                KDT.Tags[MinIdx] := KDT.Tags[I1];
                KDT.Tags[I1] := J;
            end;
            I3 := I1+1;
        end
        else
        begin
            
            //
            // 1. move split to MaxV,
            // 2. place one point to the right bin (move to I2-1),
            //    others - to the left bin
            //
            S := MaxV;
            if MaxIdx<>I2-1 then
            begin
                I:=0;
                while I<=2*KDT.NX+KDT.NY-1 do
                begin
                    V := KDT.XY[MaxIdx,I];
                    KDT.XY[MaxIdx,I] := KDT.XY[I2-1,I];
                    KDT.XY[I2-1,I] := V;
                    Inc(I);
                end;
                J := KDT.Tags[MaxIdx];
                KDT.Tags[MaxIdx] := KDT.Tags[I2-1];
                KDT.Tags[I2-1] := J;
            end;
            I3 := I2-1;
        end;
    end;
    
    //
    // Generate 'split' node
    //
    KDT.Nodes[NodesOffs+0] := 0;
    KDT.Nodes[NodesOffs+1] := D;
    KDT.Nodes[NodesOffs+2] := SplitsOffs;
    KDT.Splits[SplitsOffs+0] := S;
    OldOffs := NodesOffs;
    NodesOffs := NodesOffs+SplitNodeSize;
    SplitsOffs := SplitsOffs+1;
    
    //
    // Recirsive generation:
    // * update CurBox
    // * call subroutine
    // * restore CurBox
    //
    KDT.Nodes[OldOffs+3] := NodesOffs;
    V := KDT.CurBoxMax[D];
    KDT.CurBoxMax[D] := S;
    KDTreeGenerateTreeRec(KDT, NodesOffs, SplitsOffs, I1, I3, MaxLeafSize);
    KDT.CurBoxMax[D] := V;
    KDT.Nodes[OldOffs+4] := NodesOffs;
    V := KDT.CurBoxMin[D];
    KDT.CurBoxMin[D] := S;
    KDTreeGenerateTreeRec(KDT, NodesOffs, SplitsOffs, I3, I2, MaxLeafSize);
    KDT.CurBoxMin[D] := V;
end;


(*************************************************************************
Recursive subroutine for NN queries.

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************)
procedure KDTreeQueryNNRec(var KDT : KDTree; Offs : AlglibInteger);
var
    PtDist : Double;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    TI : AlglibInteger;
    NX : AlglibInteger;
    I1 : AlglibInteger;
    I2 : AlglibInteger;
    K1 : AlglibInteger;
    K2 : AlglibInteger;
    R1 : Double;
    R2 : Double;
    D : AlglibInteger;
    S : Double;
    V : Double;
    T1 : Double;
    ChildBestOffs : AlglibInteger;
    ChildWorstOffs : AlglibInteger;
    ChildOffs : AlglibInteger;
    PrevDist : Double;
    ToDive : Boolean;
    BestIsLeft : Boolean;
    UpdateMin : Boolean;
begin
    
    //
    // Leaf node.
    // Process points.
    //
    if KDT.Nodes[Offs]>0 then
    begin
        I1 := KDT.Nodes[Offs+1];
        I2 := I1+KDT.Nodes[Offs];
        I:=I1;
        while I<=I2-1 do
        begin
            
            //
            // Calculate distance
            //
            PtDist := 0;
            NX := KDT.NX;
            if KDT.NormType=0 then
            begin
                J:=0;
                while J<=NX-1 do
                begin
                    PtDist := Max(PtDist, AbsReal(KDT.XY[I,J]-KDT.X[J]));
                    Inc(J);
                end;
            end;
            if KDT.NormType=1 then
            begin
                J:=0;
                while J<=NX-1 do
                begin
                    PtDist := PtDist+AbsReal(KDT.XY[I,J]-KDT.X[J]);
                    Inc(J);
                end;
            end;
            if KDT.NormType=2 then
            begin
                J:=0;
                while J<=NX-1 do
                begin
                    PtDist := PtDist+AP_Sqr(KDT.XY[I,J]-KDT.X[J]);
                    Inc(J);
                end;
            end;
            
            //
            // Skip points with zero distance if self-matches are turned off
            //
            if AP_FP_Eq(PtDist,0) and  not KDT.SelfMatch then
            begin
                Inc(I);
                Continue;
            end;
            
            //
            // We CAN'T process point if R-criterion isn't satisfied,
            // i.e. (RNeeded<>0) AND (PtDist>R).
            //
            if AP_FP_Eq(KDT.RNeeded,0) or AP_FP_Less_Eq(PtDist,KDT.RNeeded) then
            begin
                
                //
                // R-criterion is satisfied, we must either:
                // * replace worst point, if (KNeeded<>0) AND (KCur=KNeeded)
                //   (or skip, if worst point is better)
                // * add point without replacement otherwise
                //
                if (KDT.KCur<KDT.KNeeded) or (KDT.KNeeded=0) then
                begin
                    
                    //
                    // add current point to heap without replacement
                    //
                    TagHeapPushI(KDT.R, KDT.Idx, KDT.KCur, PtDist, I);
                end
                else
                begin
                    
                    //
                    // New points are added or not, depending on their distance.
                    // If added, they replace element at the top of the heap
                    //
                    if AP_FP_Less(PtDist,KDT.R[0]) then
                    begin
                        if KDT.KNeeded=1 then
                        begin
                            KDT.Idx[0] := I;
                            KDT.R[0] := PtDist;
                        end
                        else
                        begin
                            TagHeapReplaceTopI(KDT.R, KDT.Idx, KDT.KNeeded, PtDist, I);
                        end;
                    end;
                end;
            end;
            Inc(I);
        end;
        Exit;
    end;
    
    //
    // Simple split
    //
    if KDT.Nodes[Offs]=0 then
    begin
        
        //
        // Load:
        // * D  dimension to split
        // * S  split position
        //
        D := KDT.Nodes[Offs+1];
        S := KDT.Splits[KDT.Nodes[Offs+2]];
        
        //
        // Calculate:
        // * ChildBestOffs      child box with best chances
        // * ChildWorstOffs     child box with worst chances
        //
        if AP_FP_Less_Eq(KDT.X[D],S) then
        begin
            ChildBestOffs := KDT.Nodes[Offs+3];
            ChildWorstOffs := KDT.Nodes[Offs+4];
            BestIsLeft := True;
        end
        else
        begin
            ChildBestOffs := KDT.Nodes[Offs+4];
            ChildWorstOffs := KDT.Nodes[Offs+3];
            BestIsLeft := False;
        end;
        
        //
        // Navigate through childs
        //
        I:=0;
        while I<=1 do
        begin
            
            //
            // Select child to process:
            // * ChildOffs      current child offset in Nodes[]
            // * UpdateMin      whether minimum or maximum value
            //                  of bounding box is changed on update
            //
            if I=0 then
            begin
                ChildOffs := ChildBestOffs;
                UpdateMin :=  not BestIsLeft;
            end
            else
            begin
                UpdateMin := BestIsLeft;
                ChildOffs := ChildWorstOffs;
            end;
            
            //
            // Update bounding box and current distance
            //
            if UpdateMin then
            begin
                PrevDist := KDT.CurDist;
                T1 := KDT.X[D];
                V := KDT.CurBoxMin[D];
                if AP_FP_Less_Eq(T1,S) then
                begin
                    if KDT.NormType=0 then
                    begin
                        KDT.CurDist := Max(KDT.CurDist, S-T1);
                    end;
                    if KDT.NormType=1 then
                    begin
                        KDT.CurDist := KDT.CurDist-Max(V-T1, 0)+S-T1;
                    end;
                    if KDT.NormType=2 then
                    begin
                        KDT.CurDist := KDT.CurDist-AP_Sqr(Max(V-T1, 0))+AP_Sqr(S-T1);
                    end;
                end;
                KDT.CurBoxMin[D] := S;
            end
            else
            begin
                PrevDist := KDT.CurDist;
                T1 := KDT.X[D];
                V := KDT.CurBoxMax[D];
                if AP_FP_Greater_Eq(T1,S) then
                begin
                    if KDT.NormType=0 then
                    begin
                        KDT.CurDist := Max(KDT.CurDist, T1-S);
                    end;
                    if KDT.NormType=1 then
                    begin
                        KDT.CurDist := KDT.CurDist-Max(T1-V, 0)+T1-S;
                    end;
                    if KDT.NormType=2 then
                    begin
                        KDT.CurDist := KDT.CurDist-AP_Sqr(Max(T1-V, 0))+AP_Sqr(T1-S);
                    end;
                end;
                KDT.CurBoxMax[D] := S;
            end;
            
            //
            // Decide: to dive into cell or not to dive
            //
            if AP_FP_Neq(KDT.RNeeded,0) and AP_FP_Greater(KDT.CurDist,KDT.RNeeded) then
            begin
                ToDive := False;
            end
            else
            begin
                if (KDT.KCur<KDT.KNeeded) or (KDT.KNeeded=0) then
                begin
                    
                    //
                    // KCur<KNeeded (i.e. not all points are found)
                    //
                    ToDive := True;
                end
                else
                begin
                    
                    //
                    // KCur=KNeeded, decide to dive or not to dive
                    // using point position relative to bounding box.
                    //
                    ToDive := AP_FP_Less_Eq(KDT.CurDist,KDT.R[0]*KDT.ApproxF);
                end;
            end;
            if ToDive then
            begin
                KDTreeQueryNNRec(KDT, ChildOffs);
            end;
            
            //
            // Restore bounding box and distance
            //
            if UpdateMin then
            begin
                KDT.CurBoxMin[D] := V;
            end
            else
            begin
                KDT.CurBoxMax[D] := V;
            end;
            KDT.CurDist := PrevDist;
            Inc(I);
        end;
        Exit;
    end;
end;


(*************************************************************************
Copies X[] to KDT.X[]
Loads distance from X[] to bounding box.
Initializes CurBox[].

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************)
procedure KDTreeInitBox(var KDT : KDTree; const X : TReal1DArray);
var
    I : AlglibInteger;
    VX : Double;
    VMin : Double;
    VMax : Double;
begin
    
    //
    // calculate distance from point to current bounding box
    //
    KDT.CurDist := 0;
    if KDT.NormType=0 then
    begin
        I:=0;
        while I<=KDT.NX-1 do
        begin
            VX := X[I];
            VMin := KDT.BoxMin[I];
            VMax := KDT.BoxMax[I];
            KDT.X[I] := VX;
            KDT.CurBoxMin[I] := VMin;
            KDT.CurBoxMax[I] := VMax;
            if AP_FP_Less(VX,VMin) then
            begin
                KDT.CurDist := Max(KDT.CurDist, VMin-VX);
            end
            else
            begin
                if AP_FP_Greater(VX,VMax) then
                begin
                    KDT.CurDist := Max(KDT.CurDist, VX-VMax);
                end;
            end;
            Inc(I);
        end;
    end;
    if KDT.NormType=1 then
    begin
        I:=0;
        while I<=KDT.NX-1 do
        begin
            VX := X[I];
            VMin := KDT.BoxMin[I];
            VMax := KDT.BoxMax[I];
            KDT.X[I] := VX;
            KDT.CurBoxMin[I] := VMin;
            KDT.CurBoxMax[I] := VMax;
            if AP_FP_Less(VX,VMin) then
            begin
                KDT.CurDist := KDT.CurDist+VMin-VX;
            end
            else
            begin
                if AP_FP_Greater(VX,VMax) then
                begin
                    KDT.CurDist := KDT.CurDist+VX-VMax;
                end;
            end;
            Inc(I);
        end;
    end;
    if KDT.NormType=2 then
    begin
        I:=0;
        while I<=KDT.NX-1 do
        begin
            VX := X[I];
            VMin := KDT.BoxMin[I];
            VMax := KDT.BoxMax[I];
            KDT.X[I] := VX;
            KDT.CurBoxMin[I] := VMin;
            KDT.CurBoxMax[I] := VMax;
            if AP_FP_Less(VX,VMin) then
            begin
                KDT.CurDist := KDT.CurDist+AP_Sqr(VMin-VX);
            end
            else
            begin
                if AP_FP_Greater(VX,VMax) then
                begin
                    KDT.CurDist := KDT.CurDist+AP_Sqr(VX-VMax);
                end;
            end;
            Inc(I);
        end;
    end;
end;


(*************************************************************************
Returns norm_k(x)^k (root-free = faster, but preserves ordering)

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************)
function VRootFreeNorm(const X : TReal1DArray;
     N : AlglibInteger;
     NormType : AlglibInteger):Double;
var
    I : AlglibInteger;
begin
    Result := 0;
    if NormType=0 then
    begin
        Result := 0;
        I:=0;
        while I<=N-1 do
        begin
            Result := Max(Result, AbsReal(X[I]));
            Inc(I);
        end;
        Exit;
    end;
    if NormType=1 then
    begin
        Result := 0;
        I:=0;
        while I<=N-1 do
        begin
            Result := Result+AbsReal(X[I]);
            Inc(I);
        end;
        Exit;
    end;
    if NormType=2 then
    begin
        Result := 0;
        I:=0;
        while I<=N-1 do
        begin
            Result := Result+AP_Sqr(X[I]);
            Inc(I);
        end;
        Exit;
    end;
end;


(*************************************************************************
Returns norm_k(x)^k (root-free = faster, but preserves ordering)

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************)
function VRootFreeComponentNorm(X : Double; NormType : AlglibInteger):Double;
begin
    Result := 0;
    if NormType=0 then
    begin
        Result := AbsReal(X);
    end;
    if NormType=1 then
    begin
        Result := AbsReal(X);
    end;
    if NormType=2 then
    begin
        Result := AP_Sqr(X);
    end;
end;


(*************************************************************************
Returns range distance: distance from X to [A,B]

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************)
function VRangeDist(X : Double; A : Double; B : Double):Double;
begin
    if AP_FP_Less(X,A) then
    begin
        Result := A-X;
    end
    else
    begin
        if AP_FP_Greater(X,B) then
        begin
            Result := X-B;
        end
        else
        begin
            Result := 0;
        end;
    end;
end;


end.