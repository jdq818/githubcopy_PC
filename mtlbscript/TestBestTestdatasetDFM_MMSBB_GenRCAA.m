testFold='F:\Coronary_0\Coronary_Niessen\ProcessByJDQ\testing\';
trainFold='F:\Coronary_0\Coronary_Niessen\ProcessByJDQ\training\';
HighResoLinefolder='F:\Coronary_0\ZXHJDQCAEDMP\DFM_BestResults\';
Resultfolder='F:\Coronary_0\ZXHJDQCAEDMP\DFM_BestResults\';
RotSubmfolder=Resultfolder;
for unseenNUM=8:1:31
    strunseenNUM=['0' num2str(unseenNUM)];
    if unseenNUM>9, strunseenNUM=num2str(unseenNUM); end;
    for vesselNUM=0:1:2
        strvesselNUM=num2str(vesselNUM);
        
        Radio=['F:\Coronary_0\Coronary_Niessen\ProcessByJDQ\RadIo\vessel' strvesselNUM '.txt'];
         disp(['multi-model begins for dataset' strunseenNUM ' -vessel' strvesselNUM]);
        if(unseenNUM<8),%training dataset
            FourREFPoinfolder=trainFold;
            
            strmodelNUM=['0' num2str(modelNUM)];
            if modelNUM>9, strmodelNUM=num2str(modelNUM); end;
            
            FourREFPoints=[FourREFPoinfolder 'dataset' strunseenNUM '\vessel' strvesselNUM '\point'];%dataset07/vessel1/point
            FourREFLinePoints=[FourREFPoinfolder 'dataset' strunseenNUM '\vessel' strvesselNUM '\BCFromRefPoint3.vtk'];%dataset07/vessel1/point
            Ref=[FourREFPoinfolder 'dataset' strunseenNUM '\vessel' strvesselNUM '\reference.txt' ];
            
            
            %% 1,SMD select
            HighResoLine=[HighResoLinefolder 'unseen' strunseenNUM '_resultsSMD\vessel' strvesselNUM];
            Rfold=[HighResoLine];
            if ~exist(Rfold,'dir')
                a=['mkdir ' ' ' Rfold];%创建命令
                system(a) %；创建文件夹
            end
            SMDtxt=[Resultfolder 'unseen' strunseenNUM '_resultsSMD\vessel' strvesselNUM '\SMD.txt' ];
            if exist(SMDtxt,'file')
                delete(SMDtxt);
            end
            ResultSMDLine=[Rfold '\SMDLine'];
            command1=['zxhcaeSoanchFL'];
            command1=[command1 ' ' Rfold ' ' ResultSMDLine ' ' FourREFPoints ' ' Radio ' ' '-SMD'];
            cd('F:\Coronary_0\code\AboutMinimalPathForExtractArtery\coronary_extraction\package_arteryextraction\package_arteryextraction\F58taosvncopy_buildx6410\bin\Release')
            system(command1);
            %Score
            SMDLine=[ResultSMDLine '.txt'];
            
            ScoreResult=[Rfold '\ScoreSMD.txt'];
            command6=['cat2008' ' ' num2str(unseenNUM) ' ' strvesselNUM ' ' Ref ' '  SMDLine ' ' ScoreResult];
            cd('F:\Coronary_0\Coronary_Niessen\software\Release')
            system(command6);
            %% 2,SMD5R
            HighResoLine=[HighResoLinefolder 'unseen' strunseenNUM '_resultsSMD5R\vessel' strvesselNUM];
            Rfold=[HighResoLine];
            if ~exist(Rfold,'dir')
                a=['mkdir ' ' ' Rfold];%创建命令
                system(a) %；创建文件夹
            end
            ResultSMD5RLine=[HighResoLine '\SMD5RLine'];
            command1=['zxhcaeSoanchFL'];
            command1=[command1 ' ' HighResoLine ' ' ResultSMD5RLine ' ' FourREFPoints ' ' Radio ' ' '-SMD5R'];
            cd('F:\Coronary_0\code\AboutMinimalPathForExtractArtery\coronary_extraction\package_arteryextraction\package_arteryextraction\F58taosvncopy_buildx6410\bin\Release')
            system(command1);
            %SMD5R Score
            SMD5RLine=[ResultSMD5RLine '.txt'];
            ScoreResult=[Rfold '\ScoreSMD5R.txt'];
            command6=['cat2008' ' ' num2str(unseenNUM) ' ' strvesselNUM ' ' Ref ' '  SMD5RLine ' ' ScoreResult];
            cd('F:\Coronary_0\Coronary_Niessen\software\Release')
            system(command6);
            
            %% 3,SMDIR5
            HighResoLine=[HighResoLinefolder 'unseen' strunseenNUM '_resultsSMDIR5\vessel' strvesselNUM];
            Rfold=[HighResoLine];
            if ~exist(Rfold,'dir')
                a=['mkdir ' ' ' Rfold];%创建命令
                system(a) %；创建文件夹
            end
            
            
            ResultSMDIR5Line=[HighResoLine '\SMDIR5Line'];
            command1=['zxhcaeSoanchFL'];
            command1=[command1 ' ' HighResoLine ' ' ResultSMDIR5Line ' ' FourREFPoints ' ' FourREFLinePoints ' ' Radio ' ' '-SMDIR5'];
            cd('F:\Coronary_0\code\AboutMinimalPathForExtractArtery\coronary_extraction\package_arteryextraction\package_arteryextraction\F58taosvncopy_buildx6410\bin\Release')
            system(command1);
            %SMDIR5 Score
            SMDIR5Line=[ResultSMDIR5Line '.txt'];
            ScoreResult=[Rfold '\ScoreSMDIR5.txt'];
            command6=['cat2008' ' ' num2str(unseenNUM) ' ' strvesselNUM ' ' Ref ' '  SMDIR5Line ' ' ScoreResult];
            cd('F:\Coronary_0\Coronary_Niessen\software\Release')
            system(command6);
            
            %% 4, SD5RSMD
            HighResoLine=[HighResoLinefolder 'unseen' strunseenNUM '_resultsSD5RSMD\vessel' strvesselNUM];
            Rfold=[HighResoLine];
            if ~exist(Rfold,'dir')
                a=['mkdir ' ' ' Rfold];%创建命令
                system(a) %；创建文件夹
            end
            ResultSD5RSMDLine=[HighResoLine '\SD5RSMDLine'];
            command1=['zxhcaeSoanchFL'];
            command1=[command1 ' ' HighResoLine ' ' ResultSD5RSMDLine ' ' FourREFPoints ' ' Radio ' ' '-SD5RSMD'];
            cd('F:\Coronary_0\code\AboutMinimalPathForExtractArtery\coronary_extraction\package_arteryextraction\package_arteryextraction\F58taosvncopy_buildx6410\bin\Release')
            system(command1);
            %SD5RSMD Score
            SD5RSMDLine=[ResultSD5RSMDLine '.txt'];
            ScoreResult=[Rfold '\ScoreSD5RSMD.txt'];
            command6=['cat2008' ' ' num2str(unseenNUM) ' ' strvesselNUM ' ' Ref ' '  SD5RSMDLine ' ' ScoreResult];
            cd('F:\Coronary_0\Coronary_Niessen\software\Release')
            system(command6);
            
            %% 5, SMDIL
            HighResoLine=[HighResoLinefolder 'unseen' strunseenNUM '_resultsSMDIL\vessel' strvesselNUM];
            Rfold=[HighResoLine];
            if ~exist(Rfold,'dir')
                a=['mkdir ' ' ' Rfold];%创建命令
                system(a) %；创建文件夹
            end
            ResultSMDILLine=[HighResoLine '\SMDILLine'];
            command1=['zxhcaeSoanchFL'];
            command1=[command1 ' ' HighResoLine ' ' ResultSMDILLine ' ' FourREFPoints ' ' FourREFLinePoints ' ' Radio ' ' '-SMDIL'];
            cd('F:\Coronary_0\code\AboutMinimalPathForExtractArtery\coronary_extraction\package_arteryextraction\package_arteryextraction\F58taosvncopy_buildx6410\bin\Release')
            system(command1);
            %SMDIL Score
            SMDILLine=[ResultSMDILLine '.txt'];
            ScoreResult=[Rfold '\ScoreSMDIL.txt'];
            command6=['cat2008' ' ' num2str(unseenNUM) ' ' strvesselNUM ' ' Ref ' '  SMDILLine ' ' ScoreResult];
            cd('F:\Coronary_0\Coronary_Niessen\software\Release')
            system(command6);
            
        else%testing dataset
            FourREFPoinfolder=testFold;
            FourREFPoints=[FourREFPoinfolder 'dataset' strunseenNUM '\vessel' strvesselNUM '\point'];%dataset07/vessel1/point
            FourREFLinePoints=[FourREFPoinfolder 'dataset' strunseenNUM '\vessel' strvesselNUM '\BCFromRefPoint3.vtk'];%dataset07/vessel1/point
            Ref=[FourREFPoinfolder 'dataset' strunseenNUM '\vessel' strvesselNUM '\reference.txt' ];
            
            
            %% 1,SMD select
            HighResoLine=[HighResoLinefolder 'unseen' strunseenNUM '_resultsSMD\vessel' strvesselNUM];
            Rfold=[HighResoLine];
            if ~exist(Rfold,'dir')
                a=['mkdir ' ' ' Rfold];%创建命令
                system(a) %；创建文件夹
            end
            SMDtxt=[Resultfolder 'unseen' strunseenNUM '_resultsSMD\vessel' strvesselNUM '\SMD.txt' ];
            if exist(SMDtxt,'file')
                delete(SMDtxt);
            end
            ResultSMDLine=[Rfold '\SMDLine'];
            command1=['zxhcaeSoanchFL'];
            command1=[command1 ' ' Rfold ' ' ResultSMDLine ' ' FourREFPoints ' ' Radio ' ' '-SMD'];
            cd('F:\Coronary_0\code\AboutMinimalPathForExtractArtery\coronary_extraction\package_arteryextraction\package_arteryextraction\F58taosvncopy_buildx6410\bin\Release')
            system(command1);
            %genterate the result line  for RCAA
                RotSubm_MM_SMD=[RotSubmfolder 'RCAA_mm_DFM\dataset' strunseenNUM '\vessel' strvesselNUM ];
            if ~exist(RotSubm_MM_SMD,'dir')
                a=['mkdir ' ' ' RotSubm_MM_SMD];%创建命令
                system(a) %；创建文件夹
            end
                        RotSubm_MM_SMD_Line=[RotSubm_MM_SMD '\result'];
            command1=['zxhcaeSoanchFL'];
            command1=[command1 ' ' Rfold ' ' RotSubm_MM_SMD_Line ' ' FourREFPoints ' ' Radio ' ' '-SMD'];
            system(command1);
            %delete .vtk
             resultvtk=[RotSubm_MM_SMD '\result.vtk'];
            if exist(resultvtk,'file')
                delete(resultvtk);
            end
            disp(['multi-model finished dataset' strunseenNUM ' -vessel' strvesselNUM]);
        end
    end
    
end


