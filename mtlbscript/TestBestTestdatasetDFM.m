% here only for traindataset and testing dataset totally 24 cases
% generate results centerlines
%zxhcaeDMP_1.1_MLLMN(input: imageRaw (.nii),vesselness(.nii),atlas_vesselx,results,atlas_vessel,detected-ostial,atlas_mean,unseen_mean,ModelPointsIntensity,segment_num; mode;output:u1 image file)
%zxhjdqCAE(testFold,trainFold,HighImgRawfold,HighImgVslsfold,HighImgVslsu1fold,AtlasVXfold,Resultfold,trainDetOstialfold,testingDetOstialfold,Meanfold,MolPoiIntfold,DetOstialfold,SegNUM,unseenNUM,vessleNUM,modelNUM);
Bestresource='F:\Coronary_0\Coronary_Niessen\ProcessByJDQ\testcases.txt';
bestcase=load(Bestresource);


testFold='F:\Coronary_0\Coronary_Niessen\ProcessByJDQ\testing\';
trainFold='F:\Coronary_0\Coronary_Niessen\PriorModel\';
HighImgRawfold='F:\Coronary_0\Coronary_Niessen\ProcessByJDQ\mol_image\' ;

HighImgRawnewfold='F:\Coronary_0\Coronary_Niessen\ProcessByJDQ\mol_image_new\' ;
AtlasVXfold='F:\Coronary_0\Coronary_Niessen\mean_centerline_2014_01_17\';
Resultfold='F:\Coronary_0\ZXHJDQCAEDMP\DFM_BestResults\';
trainDetOstialfold='F:\Coronary_0\Coronary_Niessen\ProcessByJDQ\training\';
testingDetOstialfold='F:\Coronary_0\Exp6_detect_ostia\Ostia0007To0831\Ostia0007To0831\';
Meanfold='F:\Coronary_0\meanstd\';
MolPoiIntfold='F:\Coronary_0\Coronary_Niessen\mean_centerline_intensity\';
DetOstialfold=trainDetOstialfold;
VpOut='-N';
SegNUM=30;
% badvtxtfolder='F:\Coronary_0\ZXHJDQCAEDMP\';
% badvtxt=[badvtxtfolder '\new1results.txt'];
% badv=load(badvtxt);
if ~exist(Resultfold,'dir')
    a=['mkdir ' ' ' Resultfold];%创建命令
    system(a) %；创建文件夹
end
dataNumber=size(bestcase,1);
SMNUM=zeros(dataNumber,6);
%for i=1:1:dataNumber/2
%for i=dataNumber/2+1:1:dataNumber
for i=25:1:25
    
    unseenNUM=bestcase(i,1);
    vesselNUM=bestcase(i,2);
    if(unseenNUM>7),DetOstialfold=testingDetOstialfold;end;
    
    strunseenNUM=['0' num2str(unseenNUM)];
    if unseenNUM>9, strunseenNUM=num2str(unseenNUM); end;
    strvesselNUM=num2str(vesselNUM);
    mNUM=0;
    for modelNUM=0:1:7
        if bestcase(i,3)==1%OP
            HighImgVslsfold='F:\Coronary_0\Coronary_Niessen\ProcessByJDQ\moeb_vesselness_TestBR\OP\';
            HighImgVslsu1fold='F:\Coronary_0\Coronary_Niessen\ProcessByJDQ\moeb_u1_vesselness_TestBR\OP\';
            cd('F:\Coronary_0\ZXHJDQCAEDMP');
            zxhcaeDFM_lv(testFold,trainFold,HighImgRawfold,HighImgVslsfold,HighImgVslsu1fold,AtlasVXfold,Resultfold,trainDetOstialfold,testingDetOstialfold,Meanfold,MolPoiIntfold,DetOstialfold,SegNUM,VpOut,unseenNUM,vesselNUM,modelNUM);
        end
        if bestcase(i,3)==2%tj4
            HighImgVslsfold='F:\Coronary_0\Coronary_Niessen\ProcessByJDQ\moeb_vesselness_TestBR\tj4\';
            HighImgVslsu1fold='F:\Coronary_0\Coronary_Niessen\ProcessByJDQ\moeb_u1_vesselness_TestBR\tj4\';
            cd('F:\Coronary_0\ZXHJDQCAEDMP');
            zxhcaeDFM(testFold,trainFold,HighImgRawnewfold,HighImgVslsfold,HighImgVslsu1fold,AtlasVXfold,Resultfold,trainDetOstialfold,testingDetOstialfold,Meanfold,MolPoiIntfold,DetOstialfold,SegNUM,VpOut,unseenNUM,vesselNUM,modelNUM);
        end
        if bestcase(i,3)==3%tm6
            HighImgVslsfold='F:\Coronary_0\Coronary_Niessen\ProcessByJDQ\moeb_vesselness_TestBR\tm6\';
            HighImgVslsu1fold='F:\Coronary_0\Coronary_Niessen\ProcessByJDQ\moeb_u1_vesselness_TestBR\tm6\';
            cd('F:\Coronary_0\ZXHJDQCAEDMP');
            zxhcaeDFM(testFold,trainFold,HighImgRawnewfold,HighImgVslsfold,HighImgVslsu1fold,AtlasVXfold,Resultfold,trainDetOstialfold,testingDetOstialfold,Meanfold,MolPoiIntfold,DetOstialfold,SegNUM,VpOut,unseenNUM,vesselNUM,modelNUM);
        end
        if bestcase(i,3)==4%new1
            HighImgVslsfold='F:\Coronary_0\Coronary_Niessen\ProcessByJDQ\moeb_vesselness_TestBR\new1\';
            HighImgVslsu1fold='F:\Coronary_0\Coronary_Niessen\ProcessByJDQ\moeb_u1_vesselness_TestBR\new1\';
            %             HighImgVslsfold='F:\Coronary_0\Coronary_Niessen\ProcessByJDQ\moeb_vesselness_new\';
            % HighImgVslsu1fold='F:\Coronary_0\Coronary_Niessen\ProcessByJDQ\moeb_u1_vesselness_new\';
            
            cd('F:\Coronary_0\ZXHJDQCAEDMP');
            zxhcaeDFM(testFold,trainFold,HighImgRawnewfold,HighImgVslsfold,HighImgVslsu1fold,AtlasVXfold,Resultfold,trainDetOstialfold,testingDetOstialfold,Meanfold,MolPoiIntfold,DetOstialfold,SegNUM,VpOut,unseenNUM,vesselNUM,modelNUM);
        end
        
        if bestcase(i,3)==5%OP
            HighImgVslsfold='F:\Coronary_0\Coronary_Niessen\ProcessByJDQ\moeb_vesselness_TestBR\OP\';
            HighImgVslsu1fold='F:\Coronary_0\Coronary_Niessen\ProcessByJDQ\moeb_u1_vesselness_TestBR\OP\';
            cd('F:\Coronary_0\ZXHJDQCAEDMP');
            zxhcaeDFM_Test_FBV(testFold,trainFold,HighImgRawfold,HighImgVslsfold,HighImgVslsu1fold,AtlasVXfold,Resultfold,trainDetOstialfold,testingDetOstialfold,Meanfold,MolPoiIntfold,DetOstialfold,SegNUM,unseenNUM,vesselNUM,modelNUM);
        end
        %single model
        % mNUM=mNUM+1;
        % cd('F:\Coronary_0\ZXHJDQCAEDMP');
        %cd('F:\Coronary_0\ZXHJDQCAEDMP');%zxhjdqSMSBB(testFold,trainFold,Resultfold,Resultfold,Resultfold,unseenNUM,vesselNUM,modelNUM,mNUM);
    end
    
end
%     %multi model
%     cd('F:\Coronary_0\ZXHJDQCAEDMP');
%     zxhjdqMMSBB(testFold,trainFold,Resultfold,Resultfold,Resultfold,unseenNUM,vesselNUM);
%
%     SMNUM(i,1)=unseenNUM;
%     SMNUM(i,2)=vesselNUM;
%     %read SMD in MM-DFM
%     MMSMDtxt=[Resultfold 'unseen' strunseenNUM '_resultsSMD\vessel' strvesselNUM '\SMD.txt'];
%     MMSMD=load(MMSMDtxt);
%     smmodelNUM=MMSMD(1,1);
%     SMNUM(i,3)=MMSMD(1,1);
%     %read SMD in SM-DFM
%     strsmmodelNUM=['0' num2str(smmodelNUM)];
%     if modelNUM>9, strsmmodelNUM=num2str(smmodelNUM); end;
%     BBSMDtxt=[Resultfold 'mod' strsmmodelNUM '_to_unseen' strunseenNUM '_results\meanimg' strsmmodelNUM 'v' strvesselNUM 'model\SMD.txt' ];
%     BBSMD=load(BBSMDtxt);
%     SMNUM(i,4)=BBSMD(1,1);
%     SMNUM(i,5)=BBSMD(1,5);
%     SMNUM(i,6)=MMSMD(1,5);
