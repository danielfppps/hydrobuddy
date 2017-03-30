unit hb_analysis;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, FileUtil, LResources, Forms, Controls, Graphics, Dialogs,
  StdCtrls, Grids;

type

  { TForm11 }

  TForm11 = class(TForm)
    Button1: TButton;
    StringGrid1: TStringGrid;
    procedure Button1Click(Sender: TObject);
  private
    { private declarations }
  public
    { public declarations }
  end; 

var
  Form11: TForm11; 

implementation

{ TForm11 }

procedure TForm11.Button1Click(Sender: TObject);
begin

Form11.Visible := false ;

end;

initialization
  {$I hb_analysis.lrs}

end.

