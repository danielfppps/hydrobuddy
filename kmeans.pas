{$MODESWITCH RESULT+}
{$GOTO ON}
(*************************************************************************
Copyright (c) 2008, Sergey Bochkanov (ALGLIB project).

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
unit kmeans;
interface
uses Math, Sysutils, Ap, blas;

procedure KMeansGenerate(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     K : AlglibInteger;
     Restarts : AlglibInteger;
     var Info : AlglibInteger;
     var C : TReal2DArray;
     var XYC : TInteger1DArray);

implementation

function SelectCenterPP(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     var Centers : TReal2DArray;
     BusyCenters : TBoolean1DArray;
     CCnt : AlglibInteger;
     var D2 : TReal1DArray;
     var P : TReal1DArray;
     var Tmp : TReal1DArray):Boolean;forward;


(*************************************************************************
k-means++ clusterization

INPUT PARAMETERS:
    XY          -   dataset, array [0..NPoints-1,0..NVars-1].
    NPoints     -   dataset size, NPoints>=K
    NVars       -   number of variables, NVars>=1
    K           -   desired number of clusters, K>=1
    Restarts    -   number of restarts, Restarts>=1

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -3, if task is degenerate (number of distinct points is
                          less than K)
                    * -1, if incorrect NPoints/NFeatures/K/Restarts was passed
                    *  1, if subroutine finished successfully
    C           -   array[0..NVars-1,0..K-1].matrix whose columns store
                    cluster's centers
    XYC         -   array which contains number of clusters dataset points
                    belong to.

  -- ALGLIB --
     Copyright 21.03.2009 by Bochkanov Sergey
*************************************************************************)
procedure KMeansGenerate(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     K : AlglibInteger;
     Restarts : AlglibInteger;
     var Info : AlglibInteger;
     var C : TReal2DArray;
     var XYC : TInteger1DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    CT : TReal2DArray;
    CTBest : TReal2DArray;
    XYCBest : TInteger1DArray;
    E : Double;
    EBest : Double;
    X : TReal1DArray;
    Tmp : TReal1DArray;
    D2 : TReal1DArray;
    P : TReal1DArray;
    CSizes : TInteger1DArray;
    CBusy : TBoolean1DArray;
    V : Double;
    CClosest : AlglibInteger;
    DClosest : Double;
    WORK : TReal1DArray;
    WasChanges : Boolean;
    ZeroSizeClusters : Boolean;
    Pass : AlglibInteger;
begin
    
    //
    // Test parameters
    //
    if (NPoints<K) or (NVars<1) or (K<1) or (Restarts<1) then
    begin
        Info := -1;
        Exit;
    end;
    
    //
    // TODO: special case K=1
    // TODO: special case K=NPoints
    //
    Info := 1;
    
    //
    // Multiple passes of k-means++ algorithm
    //
    SetLength(CT, K, NVars);
    SetLength(CTBest, K, NVars);
    SetLength(XYC, NPoints);
    SetLength(XYCBest, NPoints);
    SetLength(D2, NPoints);
    SetLength(P, NPoints);
    SetLength(Tmp, NVars);
    SetLength(CSizes, K);
    SetLength(CBusy, K);
    EBest := MaxRealNumber;
    Pass:=1;
    while Pass<=Restarts do
    begin
        
        //
        // Select initial centers  using k-means++ algorithm
        // 1. Choose first center at random
        // 2. Choose next centers using their distance from centers already chosen
        //
        // Note that for performance reasons centers are stored in ROWS of CT, not
        // in columns. We'll transpose CT in the end and store it in the C.
        //
        I := RandomInteger(NPoints);
        APVMove(@CT[0][0], 0, NVars-1, @XY[I][0], 0, NVars-1);
        CBusy[0] := True;
        I:=1;
        while I<=K-1 do
        begin
            CBusy[I] := False;
            Inc(I);
        end;
        if  not SelectCenterPP(XY, NPoints, NVars, CT, CBusy, K, D2, P, Tmp) then
        begin
            Info := -3;
            Exit;
        end;
        
        //
        // Update centers:
        // 2. update center positions
        //
        while True do
        begin
            
            //
            // fill XYC with center numbers
            //
            WasChanges := False;
            I:=0;
            while I<=NPoints-1 do
            begin
                CClosest := -1;
                DClosest := MaxRealNumber;
                J:=0;
                while J<=K-1 do
                begin
                    APVMove(@Tmp[0], 0, NVars-1, @XY[I][0], 0, NVars-1);
                    APVSub(@Tmp[0], 0, NVars-1, @CT[J][0], 0, NVars-1);
                    V := APVDotProduct(@Tmp[0], 0, NVars-1, @Tmp[0], 0, NVars-1);
                    if AP_FP_Less(V,DClosest) then
                    begin
                        CClosest := J;
                        DClosest := V;
                    end;
                    Inc(J);
                end;
                if XYC[I]<>CClosest then
                begin
                    WasChanges := True;
                end;
                XYC[I] := CClosest;
                Inc(I);
            end;
            
            //
            // Update centers
            //
            J:=0;
            while J<=K-1 do
            begin
                CSizes[J] := 0;
                Inc(J);
            end;
            I:=0;
            while I<=K-1 do
            begin
                J:=0;
                while J<=NVars-1 do
                begin
                    CT[I,J] := 0;
                    Inc(J);
                end;
                Inc(I);
            end;
            I:=0;
            while I<=NPoints-1 do
            begin
                CSizes[XYC[I]] := CSizes[XYC[I]]+1;
                APVAdd(@CT[XYC[I]][0], 0, NVars-1, @XY[I][0], 0, NVars-1);
                Inc(I);
            end;
            ZeroSizeClusters := False;
            I:=0;
            while I<=K-1 do
            begin
                CBusy[I] := CSizes[I]<>0;
                ZeroSizeClusters := ZeroSizeClusters or (CSizes[I]=0);
                Inc(I);
            end;
            if ZeroSizeClusters then
            begin
                
                //
                // Some clusters have zero size - rare, but possible.
                // We'll choose new centers for such clusters using k-means++ rule
                // and restart algorithm
                //
                if  not SelectCenterPP(XY, NPoints, NVars, CT, CBusy, K, D2, P, Tmp) then
                begin
                    Info := -3;
                    Exit;
                end;
                Continue;
            end;
            J:=0;
            while J<=K-1 do
            begin
                V := AP_Double(1)/CSizes[J];
                APVMul(@CT[J][0], 0, NVars-1, V);
                Inc(J);
            end;
            
            //
            // if nothing has changed during iteration
            //
            if  not WasChanges then
            begin
                Break;
            end;
        end;
        
        //
        // 3. Calculate E, compare with best centers found so far
        //
        E := 0;
        I:=0;
        while I<=NPoints-1 do
        begin
            APVMove(@Tmp[0], 0, NVars-1, @XY[I][0], 0, NVars-1);
            APVSub(@Tmp[0], 0, NVars-1, @CT[XYC[I]][0], 0, NVars-1);
            V := APVDotProduct(@Tmp[0], 0, NVars-1, @Tmp[0], 0, NVars-1);
            E := E+V;
            Inc(I);
        end;
        if AP_FP_Less(E,EBest) then
        begin
            
            //
            // store partition.
            //
            EBest := E;
            CopyMatrix(CT, 0, K-1, 0, NVars-1, CTBest, 0, K-1, 0, NVars-1);
            I:=0;
            while I<=NPoints-1 do
            begin
                XYCBest[I] := XYC[I];
                Inc(I);
            end;
        end;
        Inc(Pass);
    end;
    
    //
    // Copy and transpose
    //
    SetLength(C, NVars-1+1, K-1+1);
    CopyAndTranspose(CTBest, 0, K-1, 0, NVars-1, C, 0, NVars-1, 0, K-1);
    I:=0;
    while I<=NPoints-1 do
    begin
        XYC[I] := XYCBest[I];
        Inc(I);
    end;
end;


(*************************************************************************
Select center for a new cluster using k-means++ rule
*************************************************************************)
function SelectCenterPP(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NVars : AlglibInteger;
     var Centers : TReal2DArray;
     BusyCenters : TBoolean1DArray;
     CCnt : AlglibInteger;
     var D2 : TReal1DArray;
     var P : TReal1DArray;
     var Tmp : TReal1DArray):Boolean;
var
    I : AlglibInteger;
    J : AlglibInteger;
    CC : AlglibInteger;
    V : Double;
    S : Double;
begin
    BusyCenters := DynamicArrayCopy(BusyCenters);
    Result := True;
    CC:=0;
    while CC<=CCnt-1 do
    begin
        if  not BusyCenters[CC] then
        begin
            
            //
            // fill D2
            //
            I:=0;
            while I<=NPoints-1 do
            begin
                D2[I] := MaxRealNumber;
                J:=0;
                while J<=CCnt-1 do
                begin
                    if BusyCenters[J] then
                    begin
                        APVMove(@Tmp[0], 0, NVars-1, @XY[I][0], 0, NVars-1);
                        APVSub(@Tmp[0], 0, NVars-1, @Centers[J][0], 0, NVars-1);
                        V := APVDotProduct(@Tmp[0], 0, NVars-1, @Tmp[0], 0, NVars-1);
                        if AP_FP_Less(V,D2[I]) then
                        begin
                            D2[I] := V;
                        end;
                    end;
                    Inc(J);
                end;
                Inc(I);
            end;
            
            //
            // calculate P (non-cumulative)
            //
            S := 0;
            I:=0;
            while I<=NPoints-1 do
            begin
                S := S+D2[I];
                Inc(I);
            end;
            if AP_FP_Eq(S,0) then
            begin
                Result := False;
                Exit;
            end;
            S := 1/S;
            APVMove(@P[0], 0, NPoints-1, @D2[0], 0, NPoints-1, S);
            
            //
            // choose one of points with probability P
            // random number within (0,1) is generated and
            // inverse empirical CDF is used to randomly choose a point.
            //
            S := 0;
            V := RandomReal;
            I:=0;
            while I<=NPoints-1 do
            begin
                S := S+P[I];
                if AP_FP_Less_Eq(V,S) or (I=NPoints-1) then
                begin
                    APVMove(@Centers[CC][0], 0, NVars-1, @XY[I][0], 0, NVars-1);
                    BusyCenters[CC] := True;
                    Break;
                end;
                Inc(I);
            end;
        end;
        Inc(CC);
    end;
end;


end.