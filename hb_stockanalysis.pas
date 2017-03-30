unit hb_stockanalysis;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, FileUtil, LResources, Forms, Controls, Graphics, Dialogs,
  StdCtrls, Grids;

type

  { TForm8 }

  TForm8 = class(TForm)
    Button1: TButton;
    Label3: TLabel;
    StringGrid1: TStringGrid;
    procedure Button1Click(Sender: TObject);
  private
    { private declarations }
  public
    { public declarations }
  end; 

var
  Form8: TForm8; 

implementation

{ TForm8 }

procedure TForm8.Button1Click(Sender: TObject);
begin

Form8.Visible := false ;

end;



initialization
  {$I hb_stockanalysis.lrs}

end.

