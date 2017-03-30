unit hb_persubstance;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, FileUtil, LResources, Forms, Controls, Graphics, Dialogs,
  StdCtrls, Grids;

type

  { TForm9 }

  TForm9 = class(TForm)
    Button1: TButton;
    StringGrid1: TStringGrid;
    procedure Button1Click(Sender: TObject);
  private
    { private declarations }
  public
    { public declarations }
  end; 

var
  Form9: TForm9; 

implementation

{ TForm9 }

procedure TForm9.Button1Click(Sender: TObject);
begin

Form9.Visible := false ;

end;


initialization
  {$I hb_persubstance.lrs}

end.

