unit testmtx1;
{$IFNDEF WIN32}
You need 32 bit delphi to run this example !
{$ENDIF}

interface

uses
  Windows, Messages, SysUtils, Classes, Graphics, Controls, Forms, Dialogs,
  StdCtrls, Buttons;

type
  TForm1 = class(TForm)
    FileLabel: TLabel;
    Label3: TLabel;
    Label4: TLabel;
    Label2: TLabel;
    lbSize: TLabel;
    lbElements: TLabel;
    lbFillins: TLabel;
    lbError: TLabel;
    Label6: TLabel;
    Resolver: TBitBtn;
    OpenDialog1: TOpenDialog;
    ChooseBtn: TBitBtn;
    TimeLabel: TLabel;
    procedure ChooseBtnClick(Sender: TObject);
    procedure ResolverClick(Sender: TObject);
  private
    MTXFileName: string;
    { Private declarations }
  public
    { Public declarations }
  end;

var
  Form1: TForm1;

implementation
uses SparSolv;
{$R *.DFM}

function Seconds: Double;
var Hour, Min, Sec, MSec: Word;
begin
  DecodeTime(Now, Hour, Min, Sec, MSec);
  Result := 3600 * Hour + 60 * Min + sec + MSec / 1000;
end;

procedure TForm1.ChooseBtnClick(Sender: TObject);
begin
  FileLabel.Caption := '';
  TimeLabel.Caption := '';
  lbSize.Caption := '';
  lbElements.Caption := '';
  lbFillins.Caption := '';
  with OpenDialog1 do begin
    Filename := '*.mtx';
    if Execute then MTXFileName := FileName;
  end;
  Resolver.Enabled := (MTXFileName <> '') and FileExists(MTXFileName);
  FileLabel.Caption := MTXFileName;
end;

procedure TForm1.ResolverClick(Sender: TObject);
var
  Reason: string;
  ErrNo1, ErrNo2, ErrNo3: Integer;
  Col, Row, N,M: Integer;
  Count: LongInt;
  Value: Double;
  Start, Stop: Double;
  ErrSum: Double;
  S: string[200];
  Infile: TextFile;
  TextBuf: array[0..High(Word)] of Char;
  Ans, RHS: array[1..20000] of Double;
  elm: integer;
  symmetric: Boolean;
label Fail, EndProg;
begin

  System.Assign(infile, MTXfilename); Reset(Infile); System.SetTextBuf(Infile, TextBuf);
  Readln(Infile, s);
  if (Pos('real symmetric', s) > 0) then
    symmetric := true
  else if (Pos('real general', s) > 0) then
    symmetric := false
  else begin
    MessageDlg('Unrecognized matrix type', mtError, [mbOK], 0);
    System.Close(Infile);
    exit;
  end;
  Readln(Infile, N, M, Count);
  if (M<>N) then begin
    MessageDlg('Expected square matrix', mtError, [mbOK], 0);
    System.Close(Infile);
    exit;
  end;
  if not InitStruc(N) then goto Fail;
  for Col := 1 to N do Ans[Col] := Col;
  for Row := 1 to N do RHS[Row] := 0.0;

 {RHS is set so that variable V will have value V}
  for Elm := 1 to Count do begin
    Readln(Infile, row, col, value);
    if not AddLHS(Row, Col, Value) then goto Fail;
    if (symmetric and (Row<>COL)) then
     if not AddLHS(Row, Col, Value) then goto Fail;
    RHS[Row] := RHS[row] + Ans[col] * value;
  end;
  System.Close(Infile);
  for Row := 1 to N do
    if not AddRHS(Row, RHS[Row]) then goto Fail;

  TimeLabel.Caption := 'Solving...Please Wait';
  lbSize.Caption := IntToStr(N);
  lbElements.Caption := IntToStr(Count);
  Application.ProcessMessages;
  Start := Seconds;
  if not Solve1 then goto Fail;
  Stop := Seconds;
  ErrSum := 0.0;
  for Row := 1 to N do begin
    if not GetAnswer(Row, Value) then goto Fail;
    ErrSum := ErrSum + Abs(Value - Row);
  end;
  Str((Stop - Start): 0: 1, S);
  TimeLabel.Caption := 'Time: ' + S + ' seconds';
  Str(ErrSum: 0, S);
  lbError.Caption := S;
  lbFillins.Caption := IntToStr(Fillins);
  goto EndProg;

  Fail:
    TimeLabel.Caption := '';
    GetErrorMsg(Reason, ErrNo1, ErrNo2, ErrNo3);
   MessageDlg('Failed:  Error '+IntToStr(ErrNo1)+'   '+ Reason
   + '   '+IntToStr(ErrNo2)+ '   '+IntToStr(ErrNo3), mtError, [mbOK], 0);
  EndProg:
    ReleaseStruc;
end;

end.

