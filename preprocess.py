import numpy as np
import nibabel as nib
from nibabel.processing import conform, resample_from_to
import SimpleITK as sitk
import os
import re
import sys
import threading

DATA="./data"
RESULTS="./results"

rec = True
if rec:
    TEMPLATERECID=("Pat5",)
    CASES=("Pat1", "Pat2", "Pat3", "Pat7", "Pat8", "Pat10", "Pat11", "Pat12", "Pat13", "Pat14", "Pat16", "Pat17", "Pat18", "Pat19", "Pat20", "Pat21", "Pat22", "Pat23", "Pat24", "Pat25", "Pat26", "Pat27", "Pat31", "Pat32", "Pat33")
else:
    TEMPLATERECID=("Pat54",)
    CASES=("Pat28", "Pat29", "Pat50", "Pat51", "Pat53", "Pat54", "Pat56", "Pat64", "Pat65", "Pat66", "Pat67", "Pat68", "Pat69", "Pat70", "Pat71", "Pat73", "Pat75", "Pat76", "Pat77", "Pat79", "Pat80", "Pat82", "Pat83", "Pat85", "Pat86")
REF=f"{DATA}/{TEMPLATERECID[0]}_T2_N4_Or_Masked.nii.gz"

OR=sitk.ReadImage(f"{DATA}/B006_LFOV_N4_Norm_Origin0.nii.gz")

N4_temp = sitk.ReadImage(f"{DATA}/{TEMPLATERECID[0]}_T2_N4.nii.gz")
mask_temp = sitk.ReadImage(f"{DATA}/{TEMPLATERECID[0]}_T2_N4_Prostate.nii.gz")

N4_temp.SetOrigin(OR.GetOrigin())
mask_temp.SetOrigin(OR.GetOrigin())

sitk.WriteImage(N4_temp, f"{DATA}/{TEMPLATERECID[0]}_T2_N4_Or.nii.gz")
sitk.WriteImage(mask_temp, f"{DATA}/{TEMPLATERECID[0]}_T2_N4_Prostate_Or.nii.gz")

N4_temp = nib.load(f"{DATA}/{TEMPLATERECID[0]}_T2_N4_Or.nii.gz")
mask_temp = nib.load(f"{DATA}/{TEMPLATERECID[0]}_T2_N4_Prostate_Or.nii.gz")
temp_iso = conform(N4_temp)
mask_iso = conform(mask_temp)

masked_temp = nib.Nifti1Image(np.array(temp_iso.dataobj) * np.array(mask_iso.dataobj), temp_iso.affine)

writer = sitk.ImageFileWriter()

nib.save(masked_temp, f"{DATA}/{TEMPLATERECID[0]}_T2_N4_Or_Masked.nii.gz")
def process(case):
    N4 = sitk.ReadImage(f"{DATA}/{case}_T2_N4.nii.gz")
    mask = sitk.ReadImage(f"{DATA}/{case}_T2_N4_Prostate.nii.gz")
    N4.SetOrigin(OR.GetOrigin())
    mask.SetOrigin(OR.GetOrigin())
    sitk.WriteImage(N4, f"{DATA}/{case}_T2_N4_Or.nii.gz")
    sitk.WriteImage(mask, f"{DATA}/{case}_T2_N4_Prostate_Or.nii.gz")

    N4 = nib.load(f"{DATA}/{case}_T2_N4_Or.nii.gz")
    mask = nib.load(f"{DATA}/{case}_T2_N4_Prostate_Or.nii.gz")

    N4 = conform(N4)
    mask = conform(mask)
    
    nib.save(mask, f"{DATA}/{case}_T2_N4_Prostate_Or_Rsmp.nii.gz")
    
    masked = nib.Nifti1Image(np.array(N4.dataobj) * np.array(mask.dataobj), N4.affine)
    nib.save(masked, f"{DATA}/{case}_T2_N4_Or_Masked.nii.gz")
    
    
    os.system(f"mkdir {RESULTS}/{case}")
    os.system(f"mkdir {RESULTS}/{case}/Rigid")
    os.system(f"mkdir {RESULTS}/{case}/Rigid/MaskedImg")
    os.system(f"mkdir {RESULTS}/{case}/Rigid/prostate")
    os.system(f"mkdir {RESULTS}/{case}/Deform")
    os.system(f"mkdir {RESULTS}/{case}/Deform/MaskedImg")
    os.system(f"mkdir {RESULTS}/{case}/Deform/prostate")
    os.system(f"mkdir {RESULTS}/registration_files")
    os.system(f"mkdir {RESULTS}/registration_files/{case}")
    os.system(f"mkdir {RESULTS}/registration_files/{case}/Rigid")
    os.system(f"mkdir {RESULTS}/registration_files/{case}/Deform")

    #print(f"elastix -f {REF} -m {DATA}/{case}_T2_N4_Or_Masked.nii.gz -p rigid.txt -out {RESULTS}/registration_files/{case}/Rigid")
    #sys.exit(0)
    
    os.system(f"elastix -f {REF} -m {DATA}/{case}_T2_N4_Or_Masked.nii.gz -p rigid.txt -out {RESULTS}/registration_files/{case}/Rigid")
    os.system(f"transformix -in {DATA}/{case}_T2_N4_Or_Masked.nii.gz -tp {RESULTS}/registration_files/{case}/Rigid/TransformParameters.0.txt -out {RESULTS}/{case}/Rigid/MaskedImg")
    os.system(f"elastix -f {REF} -m {RESULTS}/{case}/Rigid/MaskedImg/result.mha -p deform.txt -out {RESULTS}/registration_files/{case}/Deform")
    os.system(f"transformix -in {RESULTS}/{case}/Rigid/MaskedImg/result.mha -tp {RESULTS}/registration_files/{case}/Deform/TransformParameters.0.txt -out {RESULTS}/{case}/Deform/MaskedImg")
    
    #change interpolation order for mask:
    with open(f"{RESULTS}/registration_files/{case}/Rigid/TransformParameters.0.txt", "r") as tp, open(f"{RESULTS}/registration_files/{case}/Rigid/TransformParameters.0.mask.txt", "w") as tpm:
        for line in tp:
            if(line == "(FinalBSplineInterpolationOrder 3)\n"):
                tpm.write("(FinalBSplineInterpolationOrder 0)\n")
            else:
                tpm.write(line)

    with open(f"{RESULTS}/registration_files/{case}/Deform/TransformParameters.0.txt", "r") as tp, open(f"{RESULTS}/registration_files/{case}/Deform/TransformParameters.0.mask.txt", "w") as tpm:
        for line in tp:
            if(line == "(FinalBSplineInterpolationOrder 3)\n"):
                tpm.write("(FinalBSplineInterpolationOrder 0)\n")
            else:
                tpm.write(line)

    os.system(f"transformix -in {DATA}/{case}_T2_N4_Prostate_Or_Rsmp.nii.gz -tp {RESULTS}/registration_files/{case}/Rigid/TransformParameters.0.mask.txt -out {RESULTS}/{case}/Rigid/prostate")
    os.system(f"transformix -in {RESULTS}/{case}/Rigid/prostate/result.mha -tp {RESULTS}/registration_files/{case}/Deform/TransformParameters.0.mask.txt -out {RESULTS}/{case}/Deform/prostate")

    #niftyreg version:
    #os.system(f"reg_aladin -ref {REF} -flo {DATA}/{case}_T2_N4_Or_Masked.nii.gz -aff {RESULTS}/registration_files/{case}_MaskedImg_Aff.txt")
    #os.system(f"reg_f3d -ref {REF} -flo {DATA}/{case}_T2_N4_Or_Masked.nii.gz -aff {RESULTS}/registration_files/{case}_MaskedImg_Aff.txt -cpp {RESULTS}/registration_files/{case}_MaskedImg_cpp.nii.gz")

    #os.system(f"reg_resample -ref {REF} -flo {DATA}/{case}_T2_N4_Prostate_Or_Rsmp.nii.gz -cpp {RESULTS}/registration_files/{case}_MaskedImg_cpp.nii.gz -result {RESULTS}/{case}_Prostate_aff_BSpline.nii.gz -NN")
    #os.system(f"reg_resample -ref {REF} -flo {DATA}/{case}_T2_N4_Or_Masked.nii.gz -cpp {RESULTS}/registration_files/{case}_MaskedImg_cpp.nii.gz -result {RESULTS}/{case}_MaskedImg_Aff_BSpline.nii.gz -NN")


    #mask_iso = sitk.ReadImage(f"{RESULTS}/{case}_Prostate_BSpline_Iso.nii.gz")
    #masked_iso = sitk.ReadImage(f"{RESULTS}/{case}_MaskedImg_BSpline_Iso.nii.gz")
    
    mask_iso = sitk.ReadImage(f"{RESULTS}/{case}/Deform/prostate/result.mha")
    
    #smooth mask and argmax
    smoother = sitk.SmoothingRecursiveGaussianImageFilter()
    smoother.SetSigma(2)
    mask_smooth = smoother.Execute(mask_iso)

    thresh = sitk.BinaryThresholdImageFilter()
    thresh.SetLowerThreshold(0.1)
    thresh.SetUpperThreshold(2)
    thresh.SetInsideValue(1)
    thresh.SetOutsideValue(0)
    mask_thresh = thresh.Execute(mask_smooth)
    writer.SetFileName(f"{RESULTS}/{case}_Prostate_Aff_BSpline_Iso_Smoothed_ArgMax.nii.gz")
    writer.Execute(mask_thresh)

    #Convert to SDF
    sdf_filt = sitk.SignedDanielssonDistanceMapImageFilter()
    mask_sdf = sdf_filt.Execute(mask_thresh)
    writer.SetFileName(f"{RESULTS}/{case}_Prostate_Aff_BSpline_Iso_Smoothed_ArgMax_SDF.nii.gz")
    writer.Execute(mask_sdf)

process(CASES[0])
sys.exit(0)
threads = []

for case in CASES:
    threads.append(threading.Thread(target=process, args=(case,)))

while len(threads) > 0:
    print(f"{len(threads)} cases left to run")
    n = 10
    to_run = threads[:n]
    [i.start() for i in to_run]
    [i.join() for i in to_run]    
    threads = threads[n:]
