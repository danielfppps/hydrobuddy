unit hb_comparison;

{$mode objfpc}

interface

uses
  Classes, SysUtils, FileUtil, LResources, Forms, Controls, Graphics, Dialogs,
  Grids, StdCtrls;

type

  { TForm15 }

  TForm15 = class(TForm)
    Button1: TButton;
    Button2: TButton;
    StringGrid1: TStringGrid;
    procedure Button1Click(Sender: TObject);
    procedure Button2Click(Sender: TObject);
  private
    { private declarations }
  public
    { public declarations }
  end; 

var
  Form15: TForm15; 

implementation

{ TForm15 }

procedure TForm15.Button2Click(Sender: TObject);
begin
  StringGrid1.ColCount := 1 ;
end;

procedure TForm15.Button1Click(Sender: TObject);
begin
  Form15.Visible := false ;
end;

initialization
  {$I hb_comparison.lrs}

end.

