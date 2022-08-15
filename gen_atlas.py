import SimpleITK as sitk
import nipype.interfaces.fsl as fsl
import os

#DATA = "./data"
DATA="/mnt/rstor/CSE_BME_AXM788/home/soumya/DATA/CCF"
RESULTS = "./results"
RESULTS2 = "/mnt/rstor/CSE_BME_AXM788/home/rakesh/shapeAnalysisProstate/results50"

OR = f"{DATA}/Or_ref.nii.gz"

TEMPLATEID=["Pat54"]
TEMPLATERECID=["Pat1"]
REF=f"{RESULTS}/{TEMPLATEID[0]}_MaskedImg_BSpline_Iso.nii.gz"
REF_REC=f"{RESULTS}/{TEMPLATERECID[0]}_MaskedImg_BSpline_Iso.nii.gz"

#cases = ["Pat28", "Pat29"]
#cases_rec = ["Pat1", "Pat2"]

cases=["Pat28", "Pat29", "Pat50", "Pat51", "Pat53", "Pat54", "Pat56", "Pat64", "Pat65", "Pat66", "Pat67", "Pat68", "Pat69", "Pat70", "Pat71", "Pat73", "Pat75", "Pat76", "Pat77", "Pat79", "Pat80", "Pat82", "Pat83", "Pat85", "Pat86"]

cases_rec=["Pat1", "Pat2", "Pat3", "Pat7", "Pat8", "Pat10", "Pat11", "Pat12", "Pat13", "Pat14", "Pat16", "Pat17", "Pat18", "Pat19", "Pat20", "Pat21", "Pat22", "Pat23", "Pat24", "Pat25", "Pat26", "Pat27", "Pat31", "Pat32", "Pat33"]

for case in cases_rec:
    os.system(f"mkdir {RESULTS}/{case}/cross_template")
    os.system(f"mkdir {RESULTS}/{case}/cross_template/MaskedImg")
    os.system(f"mkdir {RESULTS}/{case}/cross_template/Prostate")
    #os.system(f"elastix -f {REF} -m {RESULTS}/{case}_MaskedImg_Aff_BSpline_Iso.nii.gz -p rigid.txt -out {RESULTS}/{case}/cross_template -fMask {RESULTS}/{TEMPLATEID[0]}_Prostate_Aff_BSpline_Iso.nii.gz -mMask {RESULTS}/{case}_Prostate_Aff_BSpline_Iso.nii.gz")
    os.system(f"elastix -f {REF} -m {RESULTS}/{case}_MaskedImg_BSpline_Iso.nii.gz -p rigid.txt -out {RESULTS}/{case}/cross_template")
    os.system(f"transformix -in {RESULTS}/{case}_MaskedImg_BSpline_Iso.nii.gz -out {RESULTS}/{case}/cross_template/MaskedImg -tp {RESULTS}/{case}/cross_template/TransformParameters.0.txt") 
    #change interpolation order for mask:
    with open(f"{RESULTS}/{case}/cross_template/TransformParameters.0.txt", "r") as tp, open(f"{RESULTS}/{case}/cross_template/TransformParameters.0.mask.txt", "w") as tpm:
        for line in tp:
            if(line == "(FinalBSplineInterpolationOrder 3)\n"):
                tpm.write("(FinalBSplineInterpolationOrder 0)\n")
            else:
                tpm.write(line)
    os.system(f"transformix -in {RESULTS}/{case}_Prostate_BSpline_Iso.nii.gz -out {RESULTS}/{case}/cross_template/Prostate -tp {RESULTS}/{case}/cross_template/TransformParameters.0.mask.txt") 

smoother = sitk.SmoothingRecursiveGaussianImageFilter()
smoother.SetSigma(2)

thresh = sitk.BinaryThresholdImageFilter()
thresh.SetLowerThreshold(0.1)
thresh.SetUpperThreshold(2)
thresh.SetInsideValue(1)
thresh.SetOutsideValue(0)

sdf_filt = sitk.SignedDanielssonDistanceMapImageFilter()
for case in cases_rec:
    #masked = sitk.ReadImage(f"{RESULTS}/{case}/cross_template/MaskedImg/result.mha")
    mask = sitk.ReadImage(f"{RESULTS}/{case}/cross_template/Prostate/result.mha")
    
    smooth = smoother.Execute(mask)
    threshed = thresh.Execute(smooth)
    sdf = sdf_filt.Execute(threshed)

    sitk.WriteImage(mask, f"{RESULTS}/{case}_cross_Prostate_Smoothed_ArgMax_SDF.nii.gz")

Img4d = [f"{RESULTS}/{case}_Prostate_Aff_BSpline_Iso_Smoothed_ArgMax_SDF.nii.gz" for case in cases] + [f"{RESULTS}/{case_rec}_cross_Prostate_Smoothed_ArgMax_SDF.nii.gz" for case_rec in cases_rec]

merger = fsl.Merge()
merger.inputs.in_files = Img4d
merger.inputs.merged_file = f"{RESULTS}/ProstPopMerge_New.nii.gz"
merger.inputs.dimension = 'a'
merger.inputs.output_type = "NIFTI_GZ"
print(merger.cmdline)
merged = merger.run()

#construct matrix and contrast files:
with open(f"{RESULTS}/Design_New.mat", "w") as f:
    f.write(f"/NumWaves 2\n/NumPoints {len(cases)+len(cases_rec)}\n/Matrix\n")
    f.writelines("1 0\n" for i in cases)
    f.writelines("0 1\n" for i in cases_rec)
with open(f"{RESULTS}/Design_New.con", "w") as f:
    f.write(f"/NumWaves 2\n/NumPoints 1\n/Matrix\n")
    f.write("1 -1")

#T-test

rand = fsl.Randomise(in_file=f"{RESULTS}/ProstPopMerge_New.nii.gz", base_name=f"{RESULTS}/ProstPopNormal_New", design_mat=f"{RESULTS}/Design_New.mat", tcon=f"{RESULTS}/Design_New.con", mask=f"{RESULTS}/{TEMPLATEID[0]}_Prostate_Aff_BSpline_Iso_Smoothed_ArgMax.nii.gz", num_perm=500, tfce=True)
print(rand.cmdline)
res = rand.run()

tstat = sitk.ReadImage(f"{RESULTS}/ProstPopNormal_New_tfce_corrp_tstat1.nii")
msk_img = sitk.ReadImage(f"{RESULTS}/{TEMPLATEID[0]}_Prostate_Aff_BSpline_Iso_Smoothed_ArgMax.nii.gz")

tstat.SetOrigin(msk_img.GetOrigin())
tstat.SetDirection(msk_img.GetDirection())
tstat.SetSpacing(msk_img.GetSpacing())

writer = sitk.ImageFileWriter()
writer.SetFileName(f"{RESULTS}/ProstPopNormal_New_tfce_corrp_tstat1_Rsmp.nii.gz")
writer.Execute(tstat)

