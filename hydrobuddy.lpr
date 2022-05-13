program hydrobuddy;

{$mode objfpc}{$H+}

uses
  {$IFDEF UNIX}{$IFDEF UseCThreads}
  cthreads,
  {$ENDIF}{$ENDIF}
  Interfaces, // this includes the LCL widgetset
  Classes, SysUtils, FileUtil, Forms, Controls, Graphics, Dialogs, ComCtrls,
  StdCtrls, Menus, ExtCtrls, tachartlazaruspkg, HB_Main, hb_load_salts,
  hb_newcustomsalt, densesolver, hb_addweight, hb_commercialnutrient,
  hb_waterquality, hb_insprecision, hb_stockanalysis, Dbf, db, hb_persubstance,
  hb_datasetname, hb_analysis, hb_freedom, dbf_fields, hb_ph, hb_ratios,
  hb_comparison, hb_tissue_analysis;

procedure AssignValues ;

var
MyDbf: TDbf;
Units: String;
begin

DefaultFormatSettings.DecimalSeparator := '.'    ;

   MyDbf := TDbf.Create(nil) ;
   MyDbf.FilePathFull := '';
   MyDbf.TableName := Form1.water_quality_db;
   MyDbf.Open             ;
   MyDbf.Active := true ;

          while not MyDbf.EOF do
    begin

        if MyDbf.FieldByName('Default').AsInteger = 1 then

        begin

    Form6.Edit25.text := MyDbf.FieldByName('Name').AsString;
    Form6.Edit1.text := MyDbf.FieldByName('N (NO3-)').AsString ;
    Form6.Edit3.text := MyDbf.FieldByName('P').AsString ;
    Form6.Edit2.text := MyDbf.FieldByName('K').AsString ;
    Form6.Edit4.text := MyDbf.FieldByName('Mg').AsString ;
    Form6.Edit5.text := MyDbf.FieldByName('Ca').AsString ;
    Form6.Edit6.text := MyDbf.FieldByName('S').AsString ;
    Form6.Edit7.text := MyDbf.FieldByName('Fe').AsString ;
    Form6.Edit9.text := MyDbf.FieldByName('B').AsString ;
    Form6.Edit8.text := MyDbf.FieldByName('Zn').AsString ;
    Form6.Edit10.text := MyDbf.FieldByName('Cu').AsString ;
    Form6.Edit11.text := MyDbf.FieldByName('Mo').AsString ;
    Form6.Edit12.text := MyDbf.FieldByName('Na').AsString ;
    Form6.Edit15.text := MyDbf.FieldByName('Mn').AsString ;
    Form6.Edit13.text := MyDbf.FieldByName('Si').AsString ;
    Form6.Edit14.text := MyDbf.FieldByName('Cl').AsString ;
    Form6.Edit16.text := MyDbf.FieldByName('N (NH4+)').AsString ;

        end;

        MyDbf.next;
                                           // use .next here NOT .findnext!
    end;

    MyDbf.Close ;

    MyDbf.Free ;
end ;

procedure UpdateComboBoxes ;

begin

Form1.UpdateComboBox ;
Form6.UpdateComboBox ;

end ;

procedure SetActiveTab;
begin

    // set active tab and ensure all forms are visible
    Form1.PageControl1.ActivePage := Form1.TabSheet4 ;
    Form2.Position := poMainFormCenter;
    Form3.Position := poMainFormCenter;
    Form4.Position := poMainFormCenter;
    Form5.Position := poMainFormCenter;
    Form6.Position := poMainFormCenter;
    Form7.Position := poMainFormCenter;
    Form8.Position := poMainFormCenter;
    Form9.Position := poMainFormCenter;
    Form10.Position := poMainFormCenter;
    Form11.Position := poMainFormCenter;
    Form12.Position := poMainFormCenter;

end ;

procedure CheckDatabaseFiles;
begin

     if (FileExists(Form1.formulations_db) = false) or (FileExists(Form1.water_quality_db) = false) or (FileExists(Form1.substances_db) = false) or (FileExists(Form1.substances_used_db) = false) or (FileExists(Form1.tissue_analysis_db) = false) then
     begin
          ShowMessage('Database files not found, please select HydroBuddy''s installation folder.' );
          if Form1.SelectDirectoryDialog1.Execute then
          begin
               Form1.water_quality_db := Form1.SelectDirectoryDialog1.FileName + '/' + Form1.water_quality_db ;
               Form1.formulations_db := Form1.SelectDirectoryDialog1.FileName + '/' + Form1.formulations_db;
               Form1.substances_db := Form1.SelectDirectoryDialog1.FileName + '/' +  Form1.substances_db ;
               Form1.substances_used_db := Form1.SelectDirectoryDialog1.FileName + '/' + Form1.substances_used_db;
          end;
     end;

     if (FileExists(Form1.formulations_db) = false) or (FileExists(Form1.water_quality_db) = false) or (FileExists(Form1.substances_db) = false) or (FileExists(Form1.substances_used_db) = false) or (FileExists(Form1.tissue_analysis_db) = false) then
     begin
     ShowMessage('Selected folder does not contain database files, please reinstall HydroBuddy');
     Application.Terminate;
     end;


end;

{$R *.res}

begin
  Application.Title:='HydroBuddy - an Open source nutrient calculator';
  Application.Initialize;
  Application.CreateForm(TForm1, Form1);

  {$IFDEF UNIX}
  Form1.water_quality_db := 'waterquality_unix.dbf';
  Form1.formulations_db := 'formulations_unix.dbf' ;
  Form1.substances_db := 'substances_unix.dbf';
  Form1.substances_used_db := 'substances_used_unix.dbf';
  Form1.tissue_analysis_db := 'tissue_analysis_unix.dbf' ;
  {$ENDIF}

  {$IFDEF WINDOWS}
  Form1.water_quality_db := 'waterquality_win.dbf' ;
  Form1.formulations_db := 'formulations_win.dbf' ;
  Form1.substances_db := 'substances_win.dbf'    ;
  Form1.substances_used_db := 'substances_used_win.dbf' ;
  Form1.tissue_analysis_db := 'tissue_analysis_win.dbf' ;
  {$ENDIF}

  Application.CreateForm(TForm2, Form2);
  Application.CreateForm(TForm3, Form3);
  Application.CreateForm(TForm4, Form4);
  Application.CreateForm(TForm5, Form5);
  Application.CreateForm(TForm6, Form6);
  Application.CreateForm(TForm7, Form7);
  Application.CreateForm(TForm8, Form8);
  Application.CreateForm(TForm9, Form9);
  Application.CreateForm(TForm10, Form10);
  Application.CreateForm(TForm11, Form11);
  Application.CreateForm(TForm12, Form12);
  Application.CreateForm(TForm13, Form13);
  Application.CreateForm(TForm14, Form14);
  Application.CreateForm(TForm15, Form15);

  SetActiveTab  ;
  CheckDatabaseFiles;
  AssignValues ;
  UpdateComboBoxes ;
  Form1.LoadValues;
  Application.CreateForm(TForm16, Form16);
  Application.Run;
end.

