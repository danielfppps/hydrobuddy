{ Unit to solve sparse linear systems, by Mark Horridge, Monash University. 
  You are free to use or adapt this code for any purpose. }

{ DEFINE DBUG}
{$IFDEF DBUG}
{$R+,S+,Q+,X+}
{$ELSE}
{$R-,S+,Q-,X+}
{$ENDIF}
{$X+}

{you may have to adjust these defines .........}
{$DEFINE DELPHI}   {if using Delphi}
{$DEFINE DELPHI32} {if using Delphi 2 or later or any 32-bit pascal}
{$DEFINE BIT32}    {if using any 32-bit pascal}
{x$DEFINE BIT16}   {if using 16-bit pascal; just ONE of BIT32 or BIT16 must be defined}


{ $DEFINE ALIGNEDBLOCKS} {probably only useful for 16bit pascal}


unit SparSolv;
interface

 {Unit consists of the following 5 Boolean Functions and 3 procedures.
  Each function returns True if operation was successfully completed -
  otherwise False.  If False is returned, call the procedure
  GetErrorMessage to find the reason.  Finally, call the procedure
  ReleaseStruc to return memory to the heap.}

function InitStruc(NumEq: Integer): Boolean;
 {creates and initializes sparse matrix structure - call this first.
  NumEq is number of equations/variables.}

function AddLHS(ThisEqu, ThisVar: Integer; ThisVal: Double): Boolean;
 {add an element to sparse matrix for equation ThisEqu and variable ThisVar;
  if such an entry already exists, ThisVal will be added to existing value}

function AddRHS(ThisEqu: Integer; ThisVal: Double): Boolean;
 {Set RHS for equation ThisEqu; if RHS has already been set, ThisVal will be
  added to existing value}

function Solve1: Boolean;
 {calculate solutions; sparse matrix is destroyed}

function GetAnswer(ThisVar: Integer; var ThisVal: Double): Boolean;
 {read solution for variable ThisVar - probably called for each variable in
  turn}

procedure ReleaseStruc;
 {releases remaining memory used by sparse matrix structure - call this
  last}

procedure GetErrorMsg(var S: string; var N1, N2, N3: Integer);
 {N1: error number; S: Error Description; N2, N3 : other possibly useful
  numbers}

procedure ShowMat; { displays small sparse matrix }

procedure WriteMatrixAsMTX(const filename: string);
 {writes matrix in MTX format}

var
  SparMemUsed: LongInt; {no of bytes of heap currently used by routines}
  MaxMemUsed: LongInt; {highest value reached by SparMemUsed}
  FillIns: LongInt; {no of elements added during solve}
const
{$IFDEF BIT16}
  MaxMemToUse: LongInt = 4 * 16 * High(Word); {upper limit to heap use: default 4MB}
{$ELSE}
  MaxMemToUse: LongInt = High(LongInt); {upper limit to heap use: default no limit}
{$ENDIF}

implementation // see also line 166

const
  Scaling = True; {if true, matrix is preconditioned: see below}
  Msg = False;
{$IFDEF ALIGNEDBLOCKS}
  ParaAlign = 16; {suggest for paragraph alignment - try 1 to see slowdown}
{$ENDIF}
{$IFDEF BIT32}
  MaxEq = 45000; {too big doesnt matter}
{$ELSE} {16-bit, so 64k structure limit applies}
  MaxSize = 65520; {largest variable size allowed by Turbo Pascal}
  MaxUsable = MaxSize - ParaAlign; {allows for paragraph alignment of data}
  MaxEq = (MaxUsable div SizeOf(Double)) - 1; {allows 1 extra for 0-based arrays}
 {based on above, MaxEq = 8187}
{$ENDIF}
  UValue = 0.1; {number 0<U<1; if larger, time and memory usage are
                          increased at the gain, maybe, of accuracy}

type
  RealArr = array[0..MaxEq] of Double;
  IntArr = array[0..MaxEq] of Integer;
  PRealArr = ^RealArr;
  PIntArr = ^IntArr;

  { The ElmRec record type stores a single entry in the sparse matrix.
  For efficiency reasons, it is 16 bytes long, and should be paragraph-
  aligned.  A Padding field makes up the 16 bytes.  The 'preload
  PrevPtr' trick, see below, requires that Next field be at the start
  of the record.  Hence, order of fields is important.  [Note:
  Integers are 4 bytes in Delphi 2 and above, otherwise 2 bytes.]
    By converting the Value field to a Single, the Column Field to a
  Word, and eliminating the Padding field, the record size could be
  reduced to 10 bytes.  However, accuracy would be greatly reduced.
    The nodes in use are arranged in a series of linked lists, one for each
  equation (or row).  Within each list variables (columns) always appear in
  ascending order.  Each list is terminated by a RHS or constant entry; this
  is treated as though it corresponded to an extra variable, numbered (Neq+1).}

  PElmRec = ^ElmRec;

  ElmRec = record
    Next: PElmRec; {pointer to next node in row}
    Column: Integer; {variable no. of this node}
{$IFDEF BIT16}
    Padding: Integer; {only needed if Integer is 2 bytes long}
{$ENDIF}
    Value: Double; {coefficient value}
  end;

  PtrArr = array[0..MaxEq] of PElmRec;
  PPtrArr = ^PtrArr;

  { Spare nodes, not currently in use, are linked together in one list,
  pointed to by FreePtr.  FreeCount is the number of spare nodes.  Because
  there may be very many nodes, all of the same size, they are not allocated
  directly on the heap by the System Heap Manager.  Instead a more efficient
  scheme is used.  When the list of free nodes needs to be expanded, a large
  block of memory is requested, sufficient for many nodes. Another linked list
  is used to store the addresses of these blocks, for later disposal to the
  system heap}

const
  ElmSize = SizeOf(ElmRec);
{$IFDEF ALIGNEDBLOCKS}    {see GetMemParaAlign below}
  ElmsPerBlock = (65520 - ParaAlign) div ElmSize;
{$ELSE}
  ElmsPerBlock = 65520 div ElmSize;
{$ENDIF}
  {This sets block size at just under 64k - the largest size allowed
   by 16-bit Pascal. Hence each block holds around 16000 nodes. This
   should be fine for 32-bit Pascal too}

type
  TNodeBlock = array[1..ElmsPerBlock] of ElmRec;
  PNodeBlock = ^TNodeBlock;
  PHeapRec = ^HeapRec;
  HeapRec = record
    BlockPtr: PNodeBlock; {original address of block}
    NewBlockPtr: PNodeBlock; {adjusted (aligned) address of block}
    BlockSiz: LongInt; {size of block}
    NextRec: PHeapRec; {address of next list item}
  end;

var {this list contains variables which must have unit-wide scope}
  Reason: string; {for error messages}
  ErrNo1, ErrNo2, ErrNo3: Integer; {for error messages}
  Answer: PRealArr; {holds solution}
  FreePtr: PElmRec; {points to list of free nodes}
  BlockList: PHeapRec; {points to list of allocated node blocks}
  FreeCount: Integer; {number of free nodes}
  Neq: Integer; {number of equations}
  FirstElm: PPtrArr; {array of pointers to first node in each equation}
  LastElm: PPtrArr; {array of pointers to last node in each equation}

{$IFDEF DELPHI32}{$J+}{$ENDIF} {allow writeable constants}
const
  Initialized: Boolean = False; {marks whether unit is in use}
  Solved: Boolean = False; {marks whether solution was successful}

//implementation  // with "implementation" here, all is exposed !

procedure Assert(P: Boolean; S: string);
  {used for debugging}
begin
  if P then Exit;
  WriteLn('assertion failed: ', S);
  WriteLn('hit enter to finish');
  ReadLn;
  Halt(1);
end; {Assert}

procedure RecordAllocation(const OrigP, NewP: Pointer; const Siz: LongInt);
  {prepend a new heap record to the list}
var NewRec: PHeapRec;
begin
  New(NewRec);
  NewRec^.NextRec := BlockList;
  NewRec^.BlockSiz := Siz;
  NewRec^.BlockPtr := OrigP;
  NewRec^.NewBlockPtr := NewP;
  BlockList := NewRec;
  Inc(SparMemUsed, Siz + SizeOf(HeapRec)); {tally of memory used}
  if (SparMemUsed > MaxMemUsed) then MaxMemUsed := SparMemUsed;
end;

{$IFDEF ALIGNEDBLOCKS}

{routine will be faster if all node records start on a paragraph boundary.
 Unfortunately the heap manager in Borland 16-bit Pascals did not
 always issue paragraph-aligned addresses. The GetMemParaAlign function
 corrects this problem.

 It appears that the heap manager in 32-bit Delphi and Free Pascal does
 issue paragraph-aligned addresses. So GetMemParaAlign is replaced
 by a simpler version...see below. }

function GetMemParaAlign(var P: Pointer; Size: LongInt): Boolean;
  { Returns a pointer P which is paragraph aligned (a multiple of ParaAlign)
  and which points to a new block of memory of which at least Size bytes can
  be used.  The 486 cache is well suited to paragraph aligned data.
    In more detail, GetMemParaAlign obtains a block of memory from the system
  heap which is ParaAlign bytes larger than Size.  OrigP points to this
  original block.  P is the first address after OrigP which is a multiple of
  ParaAlign.  OrigP and the allocated size are saved in a linked list for
  later release.}
var
  OrigP: Pointer;
  AllocatedSize: LongInt;
{$IFDEF BIT16}
type {used only for typecasting}
  PtrRec = record OfsWord, SegWord: Word; end;
var
  Ofset: Word;
{$ENDIF}
begin
{$IFDEF TRACE}WriteLn('Entering GetMemParaAlign'); {$ENDIF}
  P := nil;
  GetMemParaAlign := False;
  AllocatedSize := Size + ParaAlign;
{$IFNDEF DELPHI32}
  if (AllocatedSize > MaxAvail) then Exit;
{$ENDIF}
  if ((MaxMemToUse <> -1) and
    ((SparMemUsed + AllocatedSize) > MaxMemToUse)) then Exit;
  GetMem(OrigP, AllocatedSize);
  if (OrigP = nil) then Exit;
{$IFDEF BIT32}
  P := Pointer(ParaAlign + ParaAlign * (LongInt(OrigP) div ParaAlign));
{$ELSE}
  P := OrigP; {to load segment}
  Ofset := PtrRec(P).OfsWord;
  {adjust offset to paragraph boundary}
  PtrRec(P).OfsWord := ParaAlign + ParaAlign * (Ofset div ParaAlign);
{$ENDIF}
  RecordAllocation(OrigP, P, AllocatedSize);
  GetMemParaAlign := True;
end; {GetMemParaAlign}

{$ELSE}

function GetMemParaAlign(var P: Pointer; Size: LongInt): Boolean;
  { replacement GetMemParaAlign function if original gives problems.
    Returns a pointer P which points to a new block of memory.
  be used.  P and the allocated size are saved in a linked list for
  later release.}
var
  OrigP: Pointer;
begin
  P := nil;
  GetMemParaAlign := False;
  if ((MaxMemToUse <> -1) and
    ((SparMemUsed + Size) > MaxMemToUse)) then Exit;
  GetMem(OrigP, Size);
  if (OrigP = nil) then Exit;
  P := OrigP;
  RecordAllocation(OrigP, P, Size);
  GetMemParaAlign := True;
end; {GetMemParaAlign}

{$ENDIF}

function TopUpFreeList: Boolean;
type {used only for typecasting}
  PtrRec = record OfsWord, SegWord: Word; end;
var
  NewBlock: PNodeBlock;
{$IFDEF BIT16}
  Ofset: Word;
{$ENDIF}
  Count: Integer;
  P: PElmRec;
begin
{$IFDEF TRACE}WriteLn('Entering TopUpFreeList'); {$ENDIF}
  TopUpFreeList := False;
  {get new block}
  if not GetMemParaAlign(Pointer(NewBlock), SizeOf(TNodeBlock)) then Exit;
  {fill new block with linked nodes}
  P := PElmRec(NewBlock);
{$IFDEF BIT32}
  for Count := 1 to ElmsPerBlock do begin
    Inc(P); {increments P by size of ElmRec}
    NewBlock^[Count].Next := P;
  end;
{$ELSE}
  Ofset := PtrRec(P).OfsWord;
  for Count := 1 to ElmsPerBlock do begin
    Inc(Ofset, ElmSize);
    PtrRec(P).OfsWord := Ofset;
    NewBlock^[Count].Next := P;
  end;
{$ENDIF}

  {attach new block to free list}
  NewBlock^[ElmsPerBlock].Next := FreePtr; {point end of new block to FreePtr}
  FreePtr := PElmRec(NewBlock); {point FreePtr at start of new block}
  Inc(FreeCount, ElmsPerBlock);

  TopUpFreeList := True;
end; {TopUpFreeList}

procedure SetErrorMsg(S: string; N1, N2, N3: Integer);
begin
  Reason := S; ErrNo1 := N1; ErrNo2 := N2; ErrNo3 := N3;
end; {SetErrorMsg}

procedure GetErrorMsg(var S: string; var N1, N2, N3: Integer);
begin
  S := Reason; N1 := ErrNo1; N2 := ErrNo2; N3 := ErrNo3;
end; {GetErrorMsg}

procedure NotInit;
begin
  SetErrorMsg('InitStruc must be called', 1, 0, 0);
end;

function InitStruc(NumEq: Integer): Boolean;
begin
{$IFDEF TRACE}WriteLn('Entering InitStruc'); {$ENDIF}
  InitStruc := False;
  Neq := NumEq;
  Solved := False;
  SetErrorMsg('', 0, 0, 0);
  BlockList := nil;
  SparMemUsed := 0;
  FillIns := 0;
  MaxMemUsed := 0;
  if Initialized then begin
    SetErrorMsg('Initialize without releasing ', 2, 0, 0);
    Exit;
  end else Initialized := True;
  if (Neq > MaxEq) then begin
    SetErrorMsg('Too many equations ', 3, Neq, 0);
    Exit;
  end;
  if (Neq < 1) then begin
    SetErrorMsg('Too few equations ', 4, Neq, 0);
    Exit;
  end;

  if not GetMemParaAlign(Pointer(FirstElm), (1 + Neq) * SizeOf(PElmRec)) then begin
    SetErrorMsg('Out of Space', 5, 0, 0);
    Exit;
  end else FillChar(FirstElm^, (1 + Neq) * SizeOf(FreePtr), 0);

  if not GetMemParaAlign(Pointer(LastElm), (1 + Neq) * SizeOf(PElmRec)) then begin
    SetErrorMsg('Out of Space', 6, 0, 0);
    Exit;
  end else FillChar(LastElm^, (1 + Neq) * SizeOf(FreePtr), 0);

  FreePtr := nil;
  FreeCount := 0;
  InitStruc := True;
end; {InitStruc}

procedure ReleaseItem(var P: Pointer);
  {release one item from user heap}
  {note no error if P is nil or is not on user heap}
var NextPtr: PHeapRec;
begin
{$IFDEF TRACE}WriteLn('Entering ReleaseItem'); {$ENDIF}
  if (P = nil) then Exit;
  NextPtr := BlockList;
  while (NextPtr <> nil) do with NextPtr^ do begin
      if (NewBlockPtr = P) then begin
        FreeMem(BlockPtr, BlockSiz);
        Dec(SparMemUsed, BlockSiz);
        BlockPtr := nil;
        NewBlockPtr := nil;
        P := nil;
{$IFDEF TRACE}WriteLn('released an item'); {$ENDIF}
        Exit;
      end;
      NextPtr := NextPtr^.NextRec;
    end; {while}
{$IFDEF DBUG}Assert(False, '#3941'); {$ENDIF}
end; {ReleaseItem}

procedure ReleaseStruc;
var NextPtr: PHeapRec;
begin
{$IFDEF TRACE}WriteLn('Entering ReleaseStruc'); {$ENDIF}
  {get rid of user heap}
  while (BlockList <> nil) do with BlockList^ do begin
      if (BlockPtr <> nil) then begin
        FreeMem(BlockPtr, BlockSiz);
        Dec(SparMemUsed, BlockSiz);
{$IFDEF TRACE}WriteLn('releasing a block'); {$ENDIF}
      end;
      NextPtr := BlockList^.NextRec;
      Dispose(BlockList); Dec(SparMemUsed, SizeOf(HeapRec));
      BlockList := NextPtr;
    end; {while}
  Initialized := False;
end; {ReleaseStruc}

function AddElement(ThisEqu, ThisVar: Integer; ThisVal: Double): Boolean;
  {AddElement will run quicker if: for each row, you first set LHS elements in ascending
  order, THEN set the RHS}
var
  PrevPtr, ElmPtr, NewPtr: PElmRec;
begin {PivRow[Row] points to last element }
  AddElement := False;
  if (ThisEqu < 1) then begin SetErrorMsg('Row < 1', 7, ThisEqu, ThisVar); Exit; end;
  if (ThisEqu > Neq) then begin SetErrorMsg('Row > Neq', 8, ThisEqu, Neq); Exit; end;
  if (FreeCount < Neq) then begin
    if not TopUpFreeList then begin
      SetErrorMsg('Out of Space', 9, 0, 0); Exit;
    end;
  end;

  {get node from free list and set its values}
  NewPtr := FreePtr; FreePtr := FreePtr^.Next; Dec(FreeCount);
  NewPtr^.Value := ThisVal;
  NewPtr^.Column := ThisVar;
  NewPtr^.Next := nil;

  {insert node in proper place in linked list}
  if FirstElm^[ThisEqu] = nil then begin {this is first entry}
    FirstElm^[ThisEqu] := NewPtr;
    LastElm^[ThisEqu] := NewPtr;
  end
  else if (LastElm^[ThisEqu]^.Column < ThisVar) then begin
   {attach at end of row}
    LastElm^[ThisEqu]^.Next := NewPtr;
    LastElm^[ThisEqu] := NewPtr;
  end
  else begin {traverse the row till we find the right place}
    ElmPtr := FirstElm^[ThisEqu];
    PrevPtr := nil;
    while (ElmPtr^.Column < ThisVar) do begin
      PrevPtr := ElmPtr;
      ElmPtr := ElmPtr^.Next;
    end;
    if (ElmPtr^.Column = ThisVar) then begin
      ElmPtr^.Value := ElmPtr^.Value + ThisVal;
    {return node to free list}
      NewPtr^.Next := FreePtr; FreePtr := NewPtr; Inc(FreeCount);
    end
    else begin {ElmPtr^.Column > ThisVar}
    {insert new node before elmptr}
      if PrevPtr = nil then FirstElm^[ThisEqu] := NewPtr
      else PrevPtr^.Next := NewPtr;
      NewPtr^.Next := ElmPtr;
    end;
  end;
  AddElement := True;
end; {AddElement}

function AddLHS(ThisEqu, ThisVar: Integer; ThisVal: Double): Boolean;
begin
  if (ThisVal = 0.0) then begin AddLHS := True; Exit; end;
  AddLHS := False;
  if not Initialized then begin NotInit; Exit; end;
  if (ThisVar < 1) then begin SetErrorMsg('Col < 1', 10, ThisEqu, ThisVar); Exit; end;
  if (ThisVar > Neq) then begin SetErrorMsg('Col > Neq', 11, ThisEqu, ThisVar); Exit; end;
  AddLHS := AddElement(ThisEqu, ThisVar, ThisVal);
end; {AddLHS}

function AddRHS(ThisEqu: Integer; ThisVal: Double): Boolean;
begin
  AddRHS := False;
  if not Initialized then begin NotInit; Exit; end;
  AddRHS := AddElement(ThisEqu, 1 + Neq, ThisVal);
end; {AddRHS}

function GetAnswer(ThisVar: Integer; var ThisVal: Double): Boolean;
  {should fail if solve not called}
begin
  GetAnswer := False;
  if not Initialized then begin NotInit; Exit; end;
  if not Solved then begin SetErrorMsg('System not solved', 12, 0, 0); Exit; end;
  if (ThisVar < 1) then begin SetErrorMsg('VarNo < 1', 13, ThisVar, 0); Exit; end;
  if (ThisVar > Neq) then begin SetErrorMsg('VarNo >Neq', 14, ThisVar, Neq); Exit; end;
  ThisVal := Answer^[ThisVar];
  GetAnswer := True;
end; {GetAnswer}

procedure ShowMat;
var
  Row, Col, c, LastCol: Integer;
  ElmPtr: PElmRec;
begin
  for Row := 1 to Neq do begin
    ElmPtr := FirstElm^[Row];
    LastCol := 0;
    while ElmPtr <> nil do begin
      Col := ElmPtr^.Column;
      for c := (LastCol + 1) to (Col - 1) do Write('nil': 6);
      Write(ElmPtr^.Value: 6: 2);
      LastCol := Col;
      ElmPtr := ElmPtr^.Next;
    end;
    for c := (LastCol + 1) to (Neq + 1) do Write('nil': 6);
    WriteLn;
  end; {for row}
end; {showmat}

procedure WriteMatrixAsMTX(const filename: string);
var
  Row, Col, Count: Integer;
  ElmPtr: PElmRec;
  Outfile: text;
begin
  Count := 0;
  for Row := 1 to Neq do begin {count elements}
    ElmPtr := FirstElm^[Row];
    while ElmPtr <> nil do begin
      Col := ElmPtr^.Column;
      if (Col <= Neq) then Inc(Count);
      ElmPtr := ElmPtr^.Next;
    end;
  end; {for row}

  Assign(outfile, filename); Rewrite(Outfile);
  Writeln(Outfile, '%%MatrixMarket matrix coordinate real general');
  Writeln(Outfile, Neq, ' ', Neq, ' ', Count);
  for Row := 1 to Neq do begin
    ElmPtr := FirstElm^[Row];
    while ElmPtr <> nil do begin
      Col := ElmPtr^.Column;
      if (Col <= Neq) then Writeln(Outfile, Row, ' ', Col, ' ', ElmPtr^.Value: 0);
      ElmPtr := ElmPtr^.Next;
    end;
  end; {for row}
  Close(OutFile);
end; {showmat}

{All these variables are used only by Solve1, but have been declared
here so that their address is known at compile time.  This can add to
speed}
var
  PrevPtr: PElmRec;
  ElmPtr: PElmRec;
  SumTerm: Extended;
  Factor: Extended;
  RHS: Extended;
  Biggest: Extended;
  Coeff: Extended;
  PivotValue: Extended;
  BestPtr: PElmRec;
  Next_Pivot: PElmRec;
  Next_Tar: PElmRec;
  NewPtr: PElmRec;
  NextActiveRow: PIntArr;
  ColCount: PIntArr;
  RowCount: PIntArr; {for active rows=no of LHS elements; for solved rows=the variable solved for}
  ColScale: PRealArr;
  PivRow: PIntArr; {used to remember in which order pivot rows were chosen}
  MinRowCount: Integer;
  Best_Addelm: Integer;
  NextTarCol: Integer;
  NumToFind: Integer;
  AddElm: Integer;
  LastCol: Integer;
  SingleCount: Integer;
  PivotStep: Integer;
  PivotCol: Integer;
  NextPivotCol: Integer;
  LastRow: Integer;
  PrevRow, NextRow: Integer;
  PivotRow: Integer;

function Solve1: Boolean;
var Row, Col: Integer; {to avoid delphi 2 warnings}
label Fail, OutOfSpace;
begin {solve1}
{$IFDEF TRACE}WriteLn('Entering Solve1'); {$ENDIF}
  Solve1 := False;
  PivotStep := 0;
  if not Initialized then begin NotInit; goto Fail; end;
  if Solved then begin SetErrorMsg('System already solved', 15, 0, 0); goto Fail; end;

  ReleaseItem(Pointer(LastElm));

  if not GetMemParaAlign(Pointer(NextActiveRow), (1 + Neq) * SizeOf(Integer)) then goto OutOfSpace;
  if not GetMemParaAlign(Pointer(ColCount), (1 + Neq) * SizeOf(Integer)) then goto OutOfSpace;
  if not GetMemParaAlign(Pointer(RowCount), (1 + Neq) * SizeOf(Integer)) then goto OutOfSpace;
  if not GetMemParaAlign(Pointer(PivRow), (1 + Neq) * SizeOf(Integer)) then goto OutOfSpace;

  {set vectors to zero}
  FillChar(RowCount^, (1 + Neq) * SizeOf(Integer), 0);
  FillChar(PivRow^, (1 + Neq) * SizeOf(Integer), 0);
  FillChar(ColCount^, (1 + Neq) * SizeOf(Integer), 0);

{$IFDEF TRACE}WriteLn('about to set up row and column counts'); {$ENDIF}
  for Row := 1 to Neq do begin
    ElmPtr := FirstElm^[Row];
    if (ElmPtr = nil) then begin
      SetErrorMsg('Empty Row', 16, Row, 0); goto Fail;
    end;
    LastCol := 0;
    while ElmPtr <> nil do begin
      Col := ElmPtr^.Column;
      if (Col <= LastCol) then begin
        SetErrorMsg('Cols out of Order', 17, Row, Col); goto Fail;
      end
      else LastCol := Col;
      Inc(RowCount^[Row]);
{$IFDEF DBUG}Assert(Col > 0, '#2133'); {$ENDIF}
{$IFDEF DBUG}Assert(Col <= (1 + Neq), '#2134'); {$ENDIF}
      if (Col <= Neq) then Inc(ColCount^[Col]);
      ElmPtr := ElmPtr^.Next;
    end;
    if (LastCol <> (1 + Neq)) then begin
      SetErrorMsg('No RHS', 18, LastCol, 0); goto Fail;
    end;
  end; {for row}

  for Col := 1 to Neq do if (ColCount^[Col] = 0) then begin
      SetErrorMsg('Empty Col', 19, Col, 0); goto Fail;
    end;

  for Row := 1 to Neq do begin
    NextActiveRow^[Row] := Row + 1;
    if RowCount^[Row] <= 1 then begin
    {get this if you have RHS but no LHS}
      SetErrorMsg('Row without Variables', 20, Row, 0); goto Fail;
    end;
  end; {for Row:=1 to neq}
  {end setup}
{$IFDEF TRACE}WriteLn('completed setup'); {$ENDIF}

  if Scaling then begin
{$IFDEF TRACE}WriteLn('about to scale'); {$ENDIF}
   {for each row, scale all elements, including RHS, so that the largest LHS coefficent is 1.00}
    for Row := 1 to Neq do begin
    {find biggest element}
      ElmPtr := FirstElm^[Row];
      Biggest := 0.0;
      while ElmPtr <> nil do begin
        if (ElmPtr^.Column > Neq) then Break;
        if (Abs(ElmPtr^.Value) > Biggest) then Biggest := Abs(ElmPtr^.Value);
        ElmPtr := ElmPtr^.Next;
      end;
      if (Biggest = 0.0) then begin
        SetErrorMsg('All-Zero Row', 21, Row, 0); goto Fail;
      end;
    {divide each element by biggest}
      ElmPtr := FirstElm^[Row];
      while ElmPtr <> nil do begin
        ElmPtr^.Value := ElmPtr^.Value / Biggest;
        ElmPtr := ElmPtr^.Next;
      end;
    end; {for row}

    if not GetMemParaAlign(Pointer(ColScale), (1 + Neq) * SizeOf(Double)) then goto OutOfSpace;
    FillChar(ColScale^, (1 + Neq) * SizeOf(Double), 0);
   {for each LHS Col, set ColScale^[Col] to largest element in col}
    for Row := 1 to Neq do begin
      ElmPtr := FirstElm^[Row];
      while ElmPtr <> nil do begin
        Col := ElmPtr^.Column; if (Col > Neq) then Break;
        if (Abs(ElmPtr^.Value) > ColScale^[Col]) then
          ColScale^[Col] := Abs(ElmPtr^.Value);
        ElmPtr := ElmPtr^.Next;
      end;
    end; {for row}
    for Col := 1 to Neq do if (ColScale^[Col] = 0.0) then begin
        SetErrorMsg('All-Zero Column', 22, Col, 0); goto Fail;
      end;

   {for each LHS Col, divide all elements by largest element in col}
    for Row := 1 to Neq do begin
      ElmPtr := FirstElm^[Row];
      while ElmPtr <> nil do begin
        Col := ElmPtr^.Column; if (Col > Neq) then Break;
        ElmPtr^.Value := ElmPtr^.Value / ColScale^[Col];
        ElmPtr := ElmPtr^.Next;
      end;
    end; {for row}

  end; {if scaling}

  {begin pivoting}
  NextActiveRow^[0] := 1;
  NextActiveRow^[Neq] := 0;

  repeat {pivot on variables which are mentioned only once}
    PivotCol := 0;
    PrevRow := 0;
    Row := NextActiveRow^[0];
    while Row <> 0 do begin
      NextRow := NextActiveRow^[Row];
{$IFDEF DBUG}Assert(Row > 0, '#8033'); {$ENDIF}
{$IFDEF DBUG}Assert(Row <= Neq, '#9033'); {$ENDIF}
      SingleCount := 0;
      ElmPtr := FirstElm^[Row];
{$IFDEF DBUG}Assert(ElmPtr <> nil, '#34'); {$ENDIF}
{$IFDEF DBUG}Assert(ElmPtr^.Column <= Neq, '#77'); {$ENDIF}
      Col := ElmPtr^.Column;
      while (Col <= Neq) do begin
        if (ColCount^[Col] = 1) then begin
          PivotCol := Col;
          Inc(SingleCount);
        end;
        ElmPtr := ElmPtr^.Next;
{$IFDEF DBUG}Assert(ElmPtr <> nil, '#35'); {$ENDIF}
        Col := ElmPtr^.Column;
      end;
      if (SingleCount > 1) then begin
        SetErrorMsg('Two Singles', 23, Row, PivotCol); goto Fail;
      end
      else if (SingleCount = 1) then begin
        Inc(PivotStep);
        PivotRow := Row;
        ElmPtr := FirstElm^[PivotRow];
{$IFDEF DBUG}Assert(ElmPtr <> nil, '#34'); {$ENDIF}
{$IFDEF DBUG}Assert(ElmPtr^.Column <= Neq, '#77'); {$ENDIF}
        Col := ElmPtr^.Column;
        while (Col <= Neq) do begin
{$IFDEF DBUG}Assert(ColCount^[Col] > 0, '#4177'); {$ENDIF}
          Dec(ColCount^[Col]);
          ElmPtr := ElmPtr^.Next;
{$IFDEF DBUG}Assert(ElmPtr <> nil, '#35'); {$ENDIF}
          Col := ElmPtr^.Column;
        end;
        PivRow^[PivotStep] := PivotRow;
        NextActiveRow^[PrevRow] := NextActiveRow^[PivotRow];
        NextActiveRow^[PivotRow] := -1; {useful ?}
        RowCount^[PivotRow] := PivotCol; {change of meaning}
        ColCount^[PivotCol] := -1; {mark as done}
      end {if (SingleCount=1) }
      else {no Singles} PrevRow := Row;
      Row := NextRow;
    end;

  until (PivotCol = 0);

  {*************main loop}
  while PivotStep < Neq do begin
    Inc(PivotStep);
{$IFDEF TRACE}WriteLn('starting step ', PivotStep); {$ENDIF}

   { Find shortest row (PivotRow) and the preceding active row (LastRow)  }
    MinRowCount := High(Integer);
    PrevRow := 0; LastRow := 0;
    Row := NextActiveRow^[0];
{$IFDEF DBUG}Assert(Row <> 0, '#33'); {$ENDIF}
    while Row <> 0 do begin
      if RowCount^[Row] < MinRowCount then begin
        MinRowCount := RowCount^[Row];
        PivotRow := Row;
        LastRow := PrevRow;
      end;
      PrevRow := Row;
      Row := NextActiveRow^[Row];
    end;

{$IFDEF TRACE}WriteLn('Pivotrow: ', PivotRow, ' Rowcount ', MinRowCount); {$ENDIF}
   {find Biggest Element in Pivot Row}
    Biggest := -1;
    ElmPtr := FirstElm^[PivotRow];
    while (ElmPtr^.Column <= Neq) do begin
      if (Abs(ElmPtr^.Value) > Biggest) then Biggest := Abs(ElmPtr^.Value);
      ElmPtr := ElmPtr^.Next;
    end;
    if (Biggest < 0.0) then begin {row had no elements}
      SetErrorMsg('Structurally Singular', 26, PivotRow, 0); goto Fail;
    end;
    if (Biggest = 0.0) then begin
      SetErrorMsg('Numerically Singular', 24, PivotRow, 0); goto Fail;
    end;

{$IFDEF TRACE}WriteLn('Biggest was :', Biggest); {$ENDIF}

      {find element in pivotrow with sparsest column,
      as long as it has absolute value at least UValue*biggest element}
    Biggest := Biggest * UValue;
    BestPtr := nil;
    Best_Addelm := High(Integer);
    ElmPtr := FirstElm^[PivotRow];
{$IFDEF DBUG}Assert(ElmPtr <> nil, '#36'); {$ENDIF}
    while (ElmPtr^.Column <= Neq) do begin
      Col := ElmPtr^.Column;
      Dec(ColCount^[Col]);
      if (Abs(ElmPtr^.Value) >= Biggest) then begin
        AddElm := ColCount^[Col];
     {addelm*rowcount is the number of additional nonzeros which would be added}
     {if this Pivot were chosen}
        if AddElm < Best_Addelm then begin
          BestPtr := ElmPtr;
          Best_Addelm := AddElm;
        end;
      end; {if (....>=UValue)}
      ElmPtr := ElmPtr^.Next;
    end;
{$IFDEF DBUG}Assert(BestPtr <> nil, '#38'); {$ENDIF}

    PivotCol := BestPtr^.Column;
    PivotValue := BestPtr^.Value;
   {Mark Pivot Row as inactive}
    NextActiveRow^[LastRow] := NextActiveRow^[PivotRow];

   {Answer Values for use by backsub}
    PivRow^[PivotStep] := PivotRow;
   {note change of meaning}
    NextActiveRow^[PivotRow] := -1; {useful ?}
    RowCount^[PivotRow] := PivotCol; {change of meaning}
    NumToFind := ColCount^[PivotCol];
    ColCount^[PivotCol] := -1; {mark as done}

{$IFDEF TRACE}
    Writeln('Start Pivot, Pivotrow: ', PivotRow, ' Rowcount ', MinRowCount);
    WriteLn(' PivotCol: ', PivotCol, ' ColCount ', NumToFind);
{$ENDIF}

    Row := NextActiveRow^[0];
    while ((Row <> 0) and (NumToFind > 0)) do begin
    {check that FreeList has enough items for the maximum possible number of insertions}
      if (FreeCount < (Neq - PivotStep)) then if not TopUpFreeList then goto OutOfSpace;

        {preload PrevPtr so that: PrevPtr^.Next :=  FirstElm^[Row]
        this works because the Next field of an ElmRec is at the start of the record}
      PrevPtr := Addr(FirstElm^[Row]);
      ElmPtr := FirstElm^[Row];

        {search along row looking for pivot column:
         if we find a bigger column we have gone far enough}
      while (ElmPtr^.Column < PivotCol) do begin
        PrevPtr := ElmPtr;
        ElmPtr := ElmPtr^.Next;
      end;

      if (ElmPtr^.Column = PivotCol) then begin

          {current row  (called "Tar", after target)contains pivot col,
           so we must add to it a multiple of pivot row}
{$IFDEF TRACE}WriteLn('Altering Row ', Row); {$ENDIF}

        Factor := ElmPtr^.Value / PivotValue;
        Dec(NumToFind);

     {unlink pivot col from current row}
        PrevPtr^.Next := ElmPtr^.Next; Dec(RowCount^[Row]);
     {prepend discarded node to Free List}
        ElmPtr^.Next := FreePtr; FreePtr := ElmPtr; Inc(FreeCount);

        Next_Pivot := FirstElm^[PivotRow];
{$IFDEF DBUG}Assert(Next_Pivot <> nil, '#333'); {$ENDIF}
        PrevPtr := Addr(FirstElm^[Row]); {so that PrevPtr^.Next :=  FirstElm^[Row]}
        Next_Tar := FirstElm^[Row];
{$IFDEF DBUG}Assert(Next_Tar <> nil, '#334'); {$ENDIF}
        NextTarCol := Next_Tar^.Column;
        while Next_Pivot <> nil do begin
          NextPivotCol := Next_Pivot^.Column;
          if (NextPivotCol <> PivotCol) then begin

       {Skip along Tar Row until we find a column as big as NextPivotCol}
            while NextTarCol < NextPivotCol do begin
              PrevPtr := Next_Tar;
              Next_Tar := Next_Tar^.Next;
{$IFDEF DBUG}Assert(Next_Tar <> nil, '#99'); {$ENDIF}
              NextTarCol := Next_Tar^.Column
            end;

            if (NextTarCol > NextPivotCol) then begin
        {element in pivot row but not in current row: add in new element}
{$IFDEF DBUG}Assert(NextTarCol > NextPivotCol, '#69'); {$ENDIF}
        {get element from free list}
              NewPtr := FreePtr; FreePtr := FreePtr^.Next; Dec(FreeCount); Inc(Fillins);
        {set values for new item}
              NewPtr^.Value := -Factor * Next_Pivot^.Value;
              NewPtr^.Column := NextPivotCol;
              NewPtr^.Next := Next_Tar;
        {connect previous item in list to new item}
              PrevPtr^.Next := NewPtr;
        {update PrevPtr}
              PrevPtr := NewPtr;
        {update column and row counts}
              Inc(ColCount^[NextPivotCol]);
              Inc(RowCount^[Row]);
            end {if (NextTarCol > NextPivotCol)}
            else begin
        {element in pivot row and also in current row: adjust value}
{$IFDEF DBUG}Assert(NextTarCol = NextPivotCol, '#67'); {$ENDIF}
              Next_Tar^.Value := Next_Tar^.Value - Factor * Next_Pivot^.Value;
            end; {else begin}
          end; {if (NextPivotCol <> PivotCol)}
          Next_Pivot := Next_Pivot^.Next; {move along pivot row}
        end; {while Next_Pivot <>  Nil}
      end; {if (ElmPtr^.Column = PivotCol)}
      Row := NextActiveRow^[Row];
    end; {while row}
{$IFDEF DBUG}Assert(NumToFind = 0, '#66'); {$ENDIF}
  end; {main loop}

  {release un-needed vectors}
  ReleaseItem(Pointer(ColCount));
  ReleaseItem(Pointer(NextActiveRow));

  {create Answer vector}
  if not GetMemParaAlign(Pointer(Answer), (1 + Neq) * SizeOf(Double))
    then goto OutOfSpace;
{$IFDEF DBUG}for Row := 1 to Neq do Answer^[Row] := -99; {$ENDIF}

  PivotStep := 1 + Neq;
  while (PivotStep > 1) do begin
    Dec(PivotStep);
    Row := PivRow^[PivotStep];
    PivotCol := RowCount^[Row]; {note change of meaning}
    SumTerm := 0.0;
    Coeff := 0.0;
    ElmPtr := FirstElm^[Row];
{$IFDEF DBUG}Assert(ElmPtr <> nil, '#188'); {$ENDIF}
    while ElmPtr <> nil do begin
      Col := ElmPtr^.Column;
      if (Col = PivotCol) then
        Coeff := ElmPtr^.Value
      else if (Col <= Neq) then begin
{$IFDEF DBUG}Assert(Answer^[Col] <> -99, '#177'); {$ENDIF}
        SumTerm := SumTerm + Answer^[Col] * ElmPtr^.Value;
      end
      else RHS := ElmPtr^.Value;
      ElmPtr := ElmPtr^.Next;
    end; {until (elmptr=Nil)}
{$IFDEF DBUG}Assert(Answer^[PivotCol] = -99, '#77'); {$ENDIF}
    Answer^[PivotCol] := (RHS - SumTerm) / Coeff;
  end; {for PivotRow:=neq downto 1}

  if Scaling then for Col := 1 to Neq do Answer^[Col] := Answer^[Col] / ColScale^[Col];

  Solved := True;
  Solve1 := True;
  Exit; {normal exit}

  {exit for error conditions}
  OutOfSpace: SetErrorMsg('Out of Space', 25, PivotStep, 0); Exit;

  Fail:

end; {Solve1}

end.

