unit hb_ratios;

{$mode objfpc}

interface

uses
  Classes, SysUtils, FileUtil, LResources, Forms, Controls, Graphics, Dialogs,
  StdCtrls, Grids;

type

  { TForm14 }

  TForm14 = class(TForm)
    Button1: TButton;
    Button2: TButton;
    ComboBox1: TComboBox;
    ComboBox2: TComboBox;
    ComboBox3: TComboBox;
    Label3: TLabel;
    Label4: TLabel;
    StringGrid1: TStringGrid;
    procedure Button1Click(Sender: TObject);
    procedure Button2Click(Sender: TObject);
  private
    { private declarations }
  public
    { public declarations }
  end; 

var
  Form14: TForm14; 

implementation

{ TForm14 }

uses HB_Main ;

procedure TForm14.Button2Click(Sender: TObject);
begin
  Form14.Visible := false ;
end;


procedure TForm14.Button1Click(Sender: TObject);
var
i : integer ;
j : integer ;
names : array[0..2] of string ;
begin

  names[0] := ComboBox1.Items[ComboBox1.ItemIndex] ;
  names[1] := ComboBox2.Items[ComboBox2.ItemIndex] ;
  names[2] := ComboBox3.Items[ComboBox3.ItemIndex] ;

  StringGrid1.Cells[0, StringGrid1.RowCount-1] := (names[0] + ': ' + names[1] + ': ' + names[2]) ;
  StringGrid1.Cells[1, StringGrid1.RowCount-2] := (HB_Main.Form1.getratio(names[0], names[1], names[2], 3)) ;

end;

initialization
  {$I hb_ratios.lrs}

end.

