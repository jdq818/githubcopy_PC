testFold='F:\Coronary_0\Coronary_Niessen\ProcessByJDQ\testing\';
trainFold='F:\Coronary_0\Coronary_Niessen\ProcessByJDQ\training\';
HighResoLinefolder='F:\Coronary_0\ZXHJDQCAEDMP\DFM_BestResults\';
Resultfolder='F:\Coronary_0\ZXHJDQCAEDMP\DFM_BestResults\';

for unseenNUM=16:1:16
    strunseenNUM=['0' num2str(unseenNUM)];
    if unseenNUM>9, strunseenNUM=num2str(unseenNUM); end;
    for vesselNUM=1:1:1
        strvesselNUM=num2str(vesselNUM);
        mNUM=0;
        Radio=['F:\Coronary_0\Coronary_Niessen\ProcessByJDQ\RadIo\vessel' strvesselNUM '.txt'];
        for modelNUM=0:1:7
            if modelNUM==unseenNUM, continue;end;
            if(unseenNUM<8),%training dataset
                FourREFPoinfolder=trainFold;
                mNUM=mNUM+1;
                strmodelNUM=['0' num2str(modelNUM)];
                if modelNUM>9, strmodelNUM=num2str(modelNUM); end;
                HighResoLine=[HighResoLinefolder 'mod' strmodelNUM '_to_unseen' strunseenNUM '_results\meanimg' strmodelNUM 'v' strvesselNUM 'model'];%unseen00_resultsDUF_new/vessel0/MCLine
                FourREFPoints=[FourREFPoinfolder 'dataset' strunseenNUM '\vessel' strvesselNUM '\point'];%dataset07/vessel1/point
                FourREFLinePoints=[FourREFPoinfolder 'dataset' strunseenNUM '\vessel' strvesselNUM '\BCFromRefPoint3.vtk'];%dataset07/vessel1/point
               Ref=[FourREFPoinfolder 'dataset' strunseenNUM '\vessel' strvesselNUM '\reference.txt' ];
                Rfold=[HighResoLine];
                if ~exist(Rfold,'dir')
                    a=['mkdir ' ' ' Rfold];%创建命令
                    system(a) %；创建文件夹
                end
                 
                %% 1,SMD select
                SMDtxt=[Resultfolder 'mod' strmodelNUM '_to_unseen' strunseenNUM '_results\meanimg' strmodelNUM 'v' strvesselNUM 'model\SMD.txt' ];
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
                %copy
                RfoldMM=[HighResoLinefolder 'unseen' strunseenNUM '_resultsSMD\vessel' strvesselNUM];
                if ~exist(RfoldMM,'dir')
                    a=['mkdir ' ' ' RfoldMM];%创建命令
                    system(a) %；创建文件夹
                end
                ResultSMDLineMM=[RfoldMM '\MCLine' num2str(mNUM-1)];
                command2=['zxhcaeSoanchFL'];
                command2=[command2 ' ' Rfold ' ' ResultSMDLineMM ' ' FourREFPoints ' ' Radio ' ' '-SMD'];
                cd('F:\Coronary_0\code\AboutMinimalPathForExtractArtery\coronary_extraction\package_arteryextraction\package_arteryextraction\F58taosvncopy_buildx6410\bin\Release')
                system(command2);
                %% 2,SMD5R
             
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
                 %copy
                RfoldMM=[HighResoLinefolder 'unseen' strunseenNUM '_resultsSMD5R\vessel' strvesselNUM];
                if ~exist(RfoldMM,'dir')
                    a=['mkdir ' ' ' RfoldMM];%创建命令
                    system(a) %；创建文件夹
                end
                ResultSMD5RLineMM=[RfoldMM '\MCLine' num2str(mNUM-1)];
                command2=['zxhcaeSoanchFL'];
                command2=[command2 ' ' HighResoLine ' ' ResultSMD5RLineMM ' ' FourREFPoints ' ' Radio ' ' '-SMD5R'];
                   cd('F:\Coronary_0\code\AboutMinimalPathForExtractArtery\coronary_extraction\package_arteryextraction\package_arteryextraction\F58taosvncopy_buildx6410\bin\Release')
                system(command2);
                %% 3,SMDIR5
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
                 %copy
                RfoldMM=[HighResoLinefolder 'unseen' strunseenNUM '_resultsSMDIR5\vessel' strvesselNUM];
                if ~exist(RfoldMM,'dir')
                    a=['mkdir ' ' ' RfoldMM];%创建命令
                    system(a) %；创建文件夹
                end
                ResultSMDIR5LineMM=[RfoldMM '\MCLine' num2str(mNUM-1)];
                command2=['zxhcaeSoanchFL'];
                command2=[command2 ' ' HighResoLine ' ' ResultSMDIR5LineMM ' ' FourREFPoints ' ' FourREFLinePoints ' ' Radio ' ' '-SMDIR5'];
                cd('F:\Coronary_0\code\AboutMinimalPathForExtractArtery\coronary_extraction\package_arteryextraction\package_arteryextraction\F58taosvncopy_buildx6410\bin\Release')
                system(command2);
                %% 4, SD5RSMD
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
                 %copy
                RfoldMM=[HighResoLinefolder 'unseen' strunseenNUM '_resultsSD5RSMD\vessel' strvesselNUM];
                if ~exist(RfoldMM,'dir')
                    a=['mkdir ' ' ' RfoldMM];%创建命令
                    system(a) %；创建文件夹
                end
                ResultSD5RSMDLineMM=[RfoldMM '\MCLine' num2str(mNUM-1)];
                command2=['zxhcaeSoanchFL'];
                command2=[command2 ' ' HighResoLine ' ' ResultSD5RSMDLineMM ' ' FourREFPoints ' ' Radio ' ' '-SD5RSMD'];
                cd('F:\Coronary_0\code\AboutMinimalPathForExtractArtery\coronary_extraction\package_arteryextraction\package_arteryextraction\F58taosvncopy_buildx6410\bin\Release')
                system(command2);
                 %% 5, SMDIL
                cd('F:\Coronary_0\code\AboutMinimalPathForExtractArtery\coronary_extraction\package_arteryextraction\package_arteryextraction\F58taosvncopy_buildx6410\bin\Release')
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
                 %copy
                RfoldMM=[HighResoLinefolder 'unseen' strunseenNUM '_resultsSMDIL\vessel' strvesselNUM];
                if ~exist(RfoldMM,'dir')
                    a=['mkdir ' ' ' RfoldMM];%创建命令
                    system(a) %；创建文件夹
                end
                ResultSMDILLineMM=[RfoldMM '\MCLine' num2str(mNUM-1)];
                command2=['zxhcaeSoanchFL'];
                command2=[command2 ' ' HighResoLine ' ' ResultSMDILLineMM ' ' FourREFPoints ' ' FourREFLinePoints ' ' Radio ' ' '-SMDIL'];
                cd('F:\Coronary_0\code\AboutMinimalPathForExtractArtery\coronary_extraction\package_arteryextraction\package_arteryextraction\F58taosvncopy_buildx6410\bin\Release')
                system(command2);
            else%testing dataset
                FourREFPoinfolder=testFold;
                mNUM=mNUM+1;
                strmodelNUM=['0' num2str(modelNUM)];
                if modelNUM>9, strmodelNUM=num2str(modelNUM); end;
                HighResoLine=[HighResoLinefolder 'mod' strmodelNUM '_to_unseen' strunseenNUM '_results\meanimg' strmodelNUM 'v' strvesselNUM 'model'];%unseen00_resultsDUF_new/vessel0/MCLine
                FourREFPoints=[FourREFPoinfolder 'dataset' strunseenNUM '\vessel' strvesselNUM '\point'];%dataset07/vessel1/point
                FourREFLinePoints=[FourREFPoinfolder 'dataset' strunseenNUM '\vessel' strvesselNUM '\BCFromRefPoint3.vtk'];%dataset07/vessel1/point
               Ref=[FourREFPoinfolder 'dataset' strunseenNUM '\vessel' strvesselNUM '\reference.txt' ];
                Rfold=[HighResoLine];
                if ~exist(Rfold,'dir')
                    a=['mkdir ' ' ' Rfold];%创建命令
                    system(a) %；创建文件夹
                end
                 
                %% 1,SMD select
                SMDtxt=[Resultfolder 'mod' strmodelNUM '_to_unseen' strunseenNUM '_results\meanimg' strmodelNUM 'v' strvesselNUM 'model\SMD.txt' ];
                 if exist(SMDtxt,'file')
                    delete(SMDtxt);
                end
                 ResultSMDLine=[Rfold '\SMDLine'];
                command1=['zxhcaeSoanchFL'];
                command1=[command1 ' ' Rfold ' ' ResultSMDLine ' ' FourREFPoints ' ' Radio ' ' '-SMD'];
                cd('F:\Coronary_0\code\AboutMinimalPathForExtractArtery\coronary_extraction\package_arteryextraction\package_arteryextraction\F58taosvncopy_buildx6410\bin\Release')
                system(command1);
               
                %copy
                RfoldMM=[HighResoLinefolder 'unseen' strunseenNUM '_resultsSMD\vessel' strvesselNUM];
                if ~exist(RfoldMM,'dir')
                    a=['mkdir ' ' ' RfoldMM];%创建命令
                    system(a) %；创建文件夹
                end
                ResultSMDLineMM=[RfoldMM '\MCLine' num2str(mNUM-1)];
                command2=['zxhcaeSoanchFL'];
                command2=[command2 ' ' Rfold ' ' ResultSMDLineMM ' ' FourREFPoints ' ' Radio ' ' '-SMD'];
                cd('F:\Coronary_0\code\AboutMinimalPathForExtractArtery\coronary_extraction\package_arteryextraction\package_arteryextraction\F58taosvncopy_buildx6410\bin\Release')
                system(command2);
%                 %% 2,SMD5R
%              
%                 ResultSMD5RLine=[HighResoLine '\SMD5RLine'];
%                 command1=['zxhcaeSoanchFL'];
%                 command1=[command1 ' ' HighResoLine ' ' ResultSMD5RLine ' ' FourREFPoints ' ' Radio ' ' '-SMD5R'];
%                    cd('F:\Coronary_0\code\AboutMinimalPathForExtractArtery\coronary_extraction\package_arteryextraction\package_arteryextraction\F58taosvncopy_buildx6410\bin\Release')
%                 system(command1);
%               
%                  %copy
%                 RfoldMM=[HighResoLinefolder 'unseen' strunseenNUM '_resultsSMD5R\vessel' strvesselNUM];
%                 if ~exist(RfoldMM,'dir')
%                     a=['mkdir ' ' ' RfoldMM];%创建命令
%                     system(a) %；创建文件夹
%                 end
%                 ResultSMD5RLineMM=[RfoldMM '\MCLine' num2str(mNUM-1)];
%                 command2=['zxhcaeSoanchFL'];
%                 command2=[command2 ' ' HighResoLine ' ' ResultSMD5RLineMM ' ' FourREFPoints ' ' Radio ' ' '-SMD5R'];
%                    cd('F:\Coronary_0\code\AboutMinimalPathForExtractArtery\coronary_extraction\package_arteryextraction\package_arteryextraction\F58taosvncopy_buildx6410\bin\Release')
%                 system(command2);
                %% 3,SMDIR5
                ResultSMDIR5Line=[HighResoLine '\SMDIR5Line'];
                command1=['zxhcaeSoanchFL'];
                command1=[command1 ' ' HighResoLine ' ' ResultSMDIR5Line ' ' FourREFPoints ' ' FourREFLinePoints ' ' Radio ' ' '-SMDIR5'];
                 cd('F:\Coronary_0\code\AboutMinimalPathForExtractArtery\coronary_extraction\package_arteryextraction\package_arteryextraction\F58taosvncopy_buildx6410\bin\Release')
                system(command1);
             
                 %copy
                RfoldMM=[HighResoLinefolder 'unseen' strunseenNUM '_resultsSMDIR5\vessel' strvesselNUM];
                if ~exist(RfoldMM,'dir')
                    a=['mkdir ' ' ' RfoldMM];%创建命令
                    system(a) %；创建文件夹
                end
                ResultSMDIR5LineMM=[RfoldMM '\MCLine' num2str(mNUM-1)];
                command2=['zxhcaeSoanchFL'];
                command2=[command2 ' ' HighResoLine ' ' ResultSMDIR5LineMM ' ' FourREFPoints ' ' FourREFLinePoints ' ' Radio ' ' '-SMDIR5'];
                cd('F:\Coronary_0\code\AboutMinimalPathForExtractArtery\coronary_extraction\package_arteryextraction\package_arteryextraction\F58taosvncopy_buildx6410\bin\Release')
                system(command2);
                %% 4, SD5RSMD
                ResultSD5RSMDLine=[HighResoLine '\SD5RSMDLine'];
                command1=['zxhcaeSoanchFL'];
                command1=[command1 ' ' HighResoLine ' ' ResultSD5RSMDLine ' ' FourREFPoints ' ' Radio ' ' '-SD5RSMD'];
                   cd('F:\Coronary_0\code\AboutMinimalPathForExtractArtery\coronary_extraction\package_arteryextraction\package_arteryextraction\F58taosvncopy_buildx6410\bin\Release')
                system(command1);
               
                 %copy
                RfoldMM=[HighResoLinefolder 'unseen' strunseenNUM '_resultsSD5RSMD\vessel' strvesselNUM];
                if ~exist(RfoldMM,'dir')
                    a=['mkdir ' ' ' RfoldMM];%创建命令
                    system(a) %；创建文件夹
                end
                ResultSD5RSMDLineMM=[RfoldMM '\MCLine' num2str(mNUM-1)];
                command2=['zxhcaeSoanchFL'];
                command2=[command2 ' ' HighResoLine ' ' ResultSD5RSMDLineMM ' ' FourREFPoints ' ' Radio ' ' '-SD5RSMD'];
                cd('F:\Coronary_0\code\AboutMinimalPathForExtractArtery\coronary_extraction\package_arteryextraction\package_arteryextraction\F58taosvncopy_buildx6410\bin\Release')
                system(command2);
%                  %% 5, SMDIL
%                 cd('F:\Coronary_0\code\AboutMinimalPathForExtractArtery\coronary_extraction\package_arteryextraction\package_arteryextraction\F58taosvncopy_buildx6410\bin\Release')
%                 ResultSMDILLine=[HighResoLine '\SMDILLine'];
%                 command1=['zxhcaeSoanchFL'];
%                 command1=[command1 ' ' HighResoLine ' ' ResultSMDILLine ' ' FourREFPoints ' ' FourREFLinePoints ' ' Radio ' ' '-SMDIL'];
%                     cd('F:\Coronary_0\code\AboutMinimalPathForExtractArtery\coronary_extraction\package_arteryextraction\package_arteryextraction\F58taosvncopy_buildx6410\bin\Release')
%                 system(command1);
%                
%                  %copy
%                 RfoldMM=[HighResoLinefolder 'unseen' strunseenNUM '_resultsSMDIL\vessel' strvesselNUM];
%                 if ~exist(RfoldMM,'dir')
%                     a=['mkdir ' ' ' RfoldMM];%创建命令
%                     system(a) %；创建文件夹
%                 end
%                 ResultSMDILLineMM=[RfoldMM '\MCLine' num2str(mNUM-1)];
%                 command2=['zxhcaeSoanchFL'];
%                 command2=[command2 ' ' HighResoLine ' ' ResultSMDILLineMM ' ' FourREFPoints ' ' FourREFLinePoints ' ' Radio ' ' '-SMDIL'];
%                 cd('F:\Coronary_0\code\AboutMinimalPathForExtractArtery\coronary_extraction\package_arteryextraction\package_arteryextraction\F58taosvncopy_buildx6410\bin\Release')
%                 system(command2);
                
            end
            disp(['single model finished dataset' strunseenNUM ' -vessel' strvesselNUM ' by model ' strmodelNUM]);
        end
    end
    
end


