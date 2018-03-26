





<2015-9-22> 
first version that is compilable
svn co  http://code.taobao.org/svn/zxhprojcae/trunk/zxhprojCAE  F:\Coronary_0\code\AboutMinimalPathForExtractArtery\coronary_extraction\package_arteryextraction\package_arteryextraction\taosvncopy ... --username jdq818  (pwd: jdq818123 )

<2015-9-30> 
second version that is compilable  
1,Use the World coordinate to correct model,update the model point which is the endpoint of the direction vector.
2,Anywhere to calculate distance all use world coordinate.
svn co  http://code.taobao.org/svn/zxhprojcae/trunk/zxhprojCAE  C:\Users\zxh\work\0_projects\myfold ... --username jdq818  (pwd: jdq818123 )

<2015-10-16> 
third version that is compilable
1,Add the FindModelPointCbasedOnGeodesicDistance function in miiMinPathModel;
2,Cout the NMPosi and VPosi to the Posi.txt which means the model vector indexs(start and end of the vector).

<2015-10-16> 
fourth version that is compilable
1,Create an 3D array to export the Alive Points(extract from Narraw Band),Active(Narraw Band)Points, and BackTrack Points of minimal path into a "nii.gz" file.Intensity of the three different kinds of points are set to  50, 20 and 100 respectively.
2,Creat an app called zxhcaeDMPSimulation.cpp to simulate the reference as different models.
<2015-11-3>
1,modify the FastMarchingEvolution UpWind FastMarchinginit into new ones.
2,resave intensity and vesselness into m_sUnseenImgIntensity and m_sNormVes respectively and calculate the potential value in the Upwind function instead of multiplying them before FastMarchingRevolution.
<2015-11-10>
1,Add the zxhcaeDMPMultiResolution filter into the code which is complied 
2,Update the three methods called fastmarchingevolution,newfastmarchingevolution,modifiedfastmarchingevolution and so on.

<2015-11-11>
1,The MethodMC.txt  from the code before didn't output the right endpoint world position.

<2015-11-17>
1,Both of Resolution and DMP increase windowing threshold from50 to 200.
2,Simulation sets reference points intensity based on only one corresponding position.
3,Updating intensity based on mean value of a variable range.
<2015-11-18>
1,I have add the model correction without moving the model to the code, and the results of the segment length are longer than previous one.
2,In addition I have change the length calculation approach in FindminPath() which has skipped some kinds of points. 
<2015-12-02>
1,I have divided the m_nModelPos into m_nOneSegModelVectorStartPos and m_nOneSegModelVectorEndPos
<2015-12-03>
1, If the distance of Point C and endpoint of the model smaller than 1mm, the endpoint will be Point C （CorrectModelPointBasedSegLengthdontMoveModel）. Updating segment-point will find Point D backward（UpdateModelSgmtPoint）. Next model vector used for next segmentation (also last segmentation) will be Point D->Point C.
<2015-12-14>
1,Change the names of 3 apps folders.
   1)zxhcaeDMPMultiModel->zxhcaeDMPMToNewImg//map points to new image.
   2)zxhcaeDMPMultiResolution->zxhcaeDMPResoImage
   3)zxhcaeDMPSimulation->zxhcaeDMPSimuModel
   4)zxhcaeDMP->zxhcaeDMP_1.0
2,Change the deformation of the last point(5 points for deforming)in Simu-model as "foward deforming".
3, submit the zxhcaeDMP1.0
<2016-01-14>
1,add badpoints detection into zxhcaeDMP_1.0
2,zxhcaeDMP_1.0_WM is the original one that extrction didn't reach the end of reference.
<2016-01-14>

1,zxhcaeDMP_1.0_WM is the new one that extrction fastmarchingevolution has the badpointdetection and linemask.
<2016-01-18>
1,ModifiedFastMarchingEvolutiondontMoveModelwithMask()has a longger centerline results ,however, ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCwithMask() has some problems to solve;
<2016-01-23>
1,To add some function to get the average status of segmentpoint ,forwardsegmentpoint ,and segmentvector;
2,To improve the function of finding C Points which combines the tangentplane distance and nearest distance insdead of single tangentplane.
<2016-01-23>
1,After testing the zxhcaeDMP_1.0_WM to get the centerline with mask, I modified the zxhcaeDMP_1.0_WM into zxhcaeDMP_1.1 as a new one, which extracts the low resolution images doesn't use mask.
2,I have use the two kinds of methods(zxhcaeDMP_1.1,zxhcaeDMP_1.0_WM)to tesing 20 simumodel from deforming referece. However, the current code also has few bit error results.
************
1,用zxhcaeDMP_1.1 对low resolution的图像进行初步提取。结果为MCLineX.vtk
2，用命令行定位路径Paht:
F:\Coronary_0\code\AboutMinimalPathForExtractArtery\coronary_extraction\package_arteryextraction\package_arteryextraction\taosvncopy_build\bin\Debug
使用以下代码找出最长的MCLineX.vtk->MCLine.vtk
miiFindMaxPath.exe F:/Coronary_0/miiMinPathModelApp_results/mod02Simu_to_unseen02_results_physnoimgres/ori_vessel0/SimuasModel/Simu19/MCLine 20
3,用zxhcaeDMPMToNewImg 把MCLine.vtk投射到high resolution图像上->MCLine.nii.gz
4,用以下代码对当前路径下的MCLine.nii.gz进行dialate->MCLine_Mask.nii.gz
zxhimageop -int MCLine.nii.gz -o MCLine_Mask -di 15
5,用zxhcaeDMP_1.0_WM 对higt resolution图像加入MCLine_Mask进行分割->MCLineX.vtk
6,再用第2步提取最长的MCLine.vtk为最终结果
miiFindMaxPath.exe F:/Coronary_0/miiMinPathModelApp_results/mol02Reso_to_unseen02_results/Simu19/MCLine 20
可视化代码:
1,可视化第2步
zxhCardCoronaryRender F:\Coronary_0\Coronary_Niessen\ProcessByLL\training\dataset02\image02_res1mm.nii.gz -centerline 2 F:\Coronary_0\Coronary_Niessen\ProcessByLL\training\dataset02\vessel0\reference.vtk F:\Coronary_0\miiMinPathModelApp_results\mod02Simu_to_unseen02_results_physnoimgres\ori_vessel0\SimuasModel\Simu3\MCLine.vtk
zxhCardCoronaryRender F:\Coronary_0\Coronary_Niessen\ProcessByLL\training\dataset02\image02_res1mm.nii.gz -centerline 3 F:\Coronary_0\Coronary_Niessen\ProcessByLL\training\dataset02\vessel0\reference.vtk F:\Coronary_0\miiMinPathModelApp_results\mod02Simu_to_unseen02_results_physnoimgres\ori_vessel0\SimuasModel\Simu19\MCLine.vtk F:\Coronary_0\Coronary_Niessen\ProcessByJDQ\training\dataset02\vessel0\Simu_reference\SetintensiyByOnePoint\Simulation19.vtk
2，可视化第6步
zxhCardCoronaryRender F:\Coronary_0\Coronary_Niessen\ProcessByLL\training\dataset02\image02_res1mm.nii.gz -centerline 3 F:\Coronary_0\Coronary_Niessen\ProcessByLL\training\dataset02\vessel0\reference.vtk F:\Coronary_0\miiMinPathModelApp_results\mol02Reso_to_unseen02_results\Simu19\MCLine.vtk F:\Coronary_0\Coronary_Niessen\ProcessByJDQ\training\dataset02\vessel0\Simu_reference\SetintensiyByOnePoint\Simulation19.vtk
************
<2016-01-28>
1，Add trainingdata extraction into the zxhcaeDMP_1.1 and zxhcaeDMP_1.0_WM
<2016-02-02>
1，I Have tested training dataset01 ,the results of vessel0 and vessel1 are good enough ,but the vessel2 and vessel3 have a bit errors;
2,the zxhcaeDMPModelIntenGen.cpp is used to get meancenterline intensity from modelline and original image;
<2016-02-14>
1,I have used the zxhcaeDMP_1.1 to test dataset01-dataset04, and it has not so bad results.
2,zxhcaeDMP_1.1_MLL was a copy of zxhcaeDMP_1.1 and I will change and test it.
<2016-03-10>
1,final code before experiments
<2016-04-22>
1, Adjust terminating conditions about end point of the model especially point "D" selection;
2, Add the semi-automatic extraction code as "*SM";
3, *MLL and *MLLCOWM are the low-reso and high-reso code for script testing based on training data;
4, *MLLT and *WMMLLT are the low-reso and high-reso code for script testing based on testing data;
5, zxhcaeDMPSelectBestBr4T is for testing;
6, zxhcaeDMPSmooth is for smooth the result centerlines of the extraction;
7, zxhcaeDMPCMBORIIMG is used for fusion original image and upsampled results of low-resolution extraction
<2016-04-22>
1，the code is used for testing "decresing the weight of the direction and intensity", also for testing "intensity threhold for incorporating low-resolution extraction results"
2,zxhcaetransImagetool is used to transfor the world coordinate into image coordinate and output as .txt file.
<2016-04-22>
1,Daily Update the code.
<2016-04-22>
1, zxhcaeDMPProjection.cpp is created for projecting the 3D image into the first slice(2D).
2, function ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCwithMaskMLL_D and ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCwithMaskMLL_DeaF_D are added for extraction with new mask(intersection of old dilated mask and new threhold vesselness)
<2016-05-25>
the code correct the bad results than before, but it need to be correct terminate condition.
<2016-05-26>
terminate condition of zxhcaeDMP_1.1_MLLM has been corrected and tested on bad results cases.
<2016-05-30>
1,zxhcaeDMP_1.1_MLLM has been change for ‘L’ ‘H’ and 'DH' mode corresponding Low-resolution extraction, High-resolution extraction and High-resoltuion extraction directly.
2,3 kinds of terminating conditions and conditions for segmentpoint selecting have been added into evolution.
<2016-05-30>
zxhcaeDMPMToNewImg is changed into -R(map to input image ) and -N(map to totally new image) modes 
<2016-07-12>
zxhcaeDMP_1.1_MLLMN is added in order to use u1 vector to calculate direction parameter
<2016-07-22>
...NB is add into the code for chosing the maximum length by backtracking.
<2016-07-23>
zxhcaeDMP_1.1_MLLMN is updated after meeting
<2016-08-04>
zxhcaeDMP_1.1_MLLMN+SBP is added before testing.
<2016-08-08>
old random function in zxhcaeDMP_1.1_MLLMN+SBP is updated before testing.
<2016-08-09>
random function is modfied for generating not repeating numbers.
<2016-08-24>
convergence parameter is added into the code
<2016-09-11>
zxhcaeSoanch.cpp is added into the code before testing.
<2016-09-17>
zxhcaeDMP_1.1_MLLMN+SBP has been tested before running by matlab
zxhcaeDMPBesselCurve is added into the code
<2016-09-18>
bugs in if sentences are corrected inModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaFM_DN
<2016-09-24>
SelectBPasFinalSegPoint_RandomNoOv15WithoutSSP_WORec has added and parameter adjusted into the code.
<2016-09-27>
DMP debug is adjusted for synthesis data；
<2016-10-11>
JDQ re-generated vsls image and u1 image for training data;
and the code is used for testing.
<2016-10-16>
Selecting best branch code is modified for whole cases (Matlab);
<2016-11-21>
ModifiedFastMarchingEvolutiondontMoveModelwithFindMinPathPointCNewLL_DeaFM_DN_VVSel is added into the code for selecting vessel vector
<2016-12-12>
badpoint correction is modified and added  in ModifiedFindImgMidPoint_DistanceThresholdPlusVectorAngledontMoveSetPointBMLL_JDQN_VSPSC 
jdqimageop_recipr is added into the code.
<2017-1-12>
The code is used for mia paper
<2017-1-25>
multi-stecil fast marching is tested in this code.
<2017-6-9>
prepare for journal paper
<2017-06-25>
zxhcaeSoanchFL is an new and all range of techniqus for selecting  best branch.
It includes and improve the old technique from zxhcaeSoanch before and added new methods.
It still need testing on all results from trainning and test data.
<2017-07-07>
this type of the code is prepared for jounal paper.
the results by this code is submited today.
<2017-10-03> Gen_CeL_From_MS: generate centerlines from the manual segment points: start coding.
<2017-10-10>Gen_CeL_From_MS can identify the key points of the branch.
<2017-10-11>Gen_CeL_From_MS can identify the leaf points; order the points within each branch; identify the key points of the branch.
<2017-10-12>Gen_CeL_From_MS has been tested by dataset00 and its manual key points;
<2017-10-19> JDQ_Curve_Fusion is used for fusing all the curves from zxhcae.
<2017-10-22>Curfusion_jdq2017 is used for fusiong all curves from zxhcae. But it borrow the idea of Theo::cat2008
<2017-12-19>the framework of curve fusion is constructed and modified.
<2017-12-20> Planefitting is modified.
<2017-12-25> curve fitting is consistent with the plane fitting technique.
<2017-12-29> rotatation part is modified and tested right.
<2018-01-03> bifurcation point is coding
<2018-01-11> I have implemented the results fusion code successfully, although meeting many bugs. 
But the current code only uses two candidate curves. Next I will use all the results we have obtained before as soon as possible.
<2018-01-16> Read curves from a folder
<2018-01-17> Begin Collec2_ALL
<2018-01-19> The code results fusion method are tested right.
<2018-01-25>BifurDec is built for bifurcation point detection
<2018-01-30>bwlabel is used to detect connected region.
<2018-02-05>Backlabeled is added into the bwlabel for labeling the connected region in the original matrix;
<2018-02-06>BWLABEL3 is built to detect connected region in 3D case. One test has taken.
<2018-02-07>FIRSTPASS and replaceSameLabel is added into the BWLABEL3 and have tested.
<2018-02-07>BWLABEL3 is tested in BifurDec.
<2018-02-07>BWLABEL3 has a little bug to correct.
<2018-02-08>Point_select is used to detect the bifurcation points in a window.
<2018-02-08>Point_select is modified for detection bifurcation points.
<2018-02-09>I find some errors in detecting connected region.
<2018-02-11> Point_select is tested right for bifurcation points.
<2018-02-11> Point_select is tested right for root points.
<2018-02-13> this is the test for github
<2018-02-13> Point_link is first added into the code