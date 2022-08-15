import SimpleITK as sitk
import numpy as np
import os
import threading
import sys
import nibabel as nib
from nibabel.processing import conform

#DATA="./data"
DATA="/mnt/rstor/CSE_BME_AXM788/home/soumya/DATA/CCF"
RESULTS="./results"

rec = False #True to generate recurrance template, False for non-recurrance
if rec:
    REG_TEMP="Pat5"
    #CASES=("Pat1", "Pat2", "Pat3", "Pat7", "Pat8", "Pat10", "Pat11", "Pat12", "Pat13", "Pat14", "Pat16", "Pat17", "Pat18", "Pat19", "Pat20", "Pat21", "Pat22", "Pat23", "Pat24", "Pat25", "Pat26", "Pat27", "Pat31", "Pat32", "Pat33")
    CASES=("Pat1", "Pat2", "Pat8", "Pat10", "Pat11", "Pat12", "Pat13", "Pat14", "Pat16", "Pat17", "Pat18", "Pat19", "Pat20", "Pat21", "Pat22", "Pat23", "Pat26", "Pat27")
else:
    REG_TEMP="Pat54"
    #CASES=("Pat28", "Pat29", "Pat50", "Pat51", "Pat53", "Pat54", "Pat56", "Pat64", "Pat65", "Pat66", "Pat67", "Pat68", "Pat69", "Pat70", "Pat71", "Pat73", "Pat75", "Pat76", "Pat77", "Pat79", "Pat80", "Pat82", "Pat83", "Pat85", "Pat86")
    CASES=("Pat50", "Pat53", "Pat54", "Pat56", "Pat64", "Pat65", "Pat66", "Pat67", "Pat68", "Pat69", "Pat70", "Pat71", "Pat73", "Pat76", "Pat77", "Pat79", "Pat80", "Pat82", "Pat83", "Pat85", "Pat86")

REF=f"{DATA}/{REG_TEMP}_T2_N4.nii.gz"

"""
OR = sitk.ReadImage(f"/mnt/rstor/CSE_BME_AXM788/home/soumya/DATA/CCF/B006_LFOV_N4_Norm_Origin0.nii.gz")
REF_img = sitk.ReadImage(REF)
REF_img.SetOrigin(OR.GetOrigin())
REF_img.SetDirection(OR.GetDirection())
REF_img.SetSpacing((1, 1, 1))
sitk.WriteImage(REF_img, f"{RESULTS}/template/{REG_TEMP}_Rsmp.nii.gz")
"""


REF_img = nib.load(REF)
REF_mask = nib.load(f"{DATA}/{REG_TEMP}_T2_N4_Prostate.nii.gz")
shape = (256, 256, 256)
REF_img = conform(REF_img, out_shape=shape, voxel_size=(1.0, 1.0, 1.0))
REF_mask = conform(REF_mask, out_shape=shape, voxel_size=(1.0, 1.0, 1.0))
nib.save(REF_img, f"{RESULTS}/template/{REG_TEMP}_Rsmp.nii.gz")
nib.save(REF_mask, f"{RESULTS}/template/{REG_TEMP}_Prostate_Rsmp.nii.gz")


REF = f"{RESULTS}/template/{REG_TEMP}_Rsmp.nii.gz"
REFMask = f"{RESULTS}/template/{REG_TEMP}_Prostate_Rsmp.nii.gz"

os.system(f"mkdir {RESULTS}/template")
#os.system(f"mkdir {RESULTS}/template/{REG_TEMP}")

REF_img = sitk.ReadImage(REF)
REF_mask = sitk.ReadImage(REFMask)
images = [sitk.GetArrayFromImage(REF_img)]
masks = [sitk.GetArrayFromImage(REF_mask)]
#images = [np.array(REF_img)]

def process(case):
    os.system(f"mkdir {RESULTS}/template/{case}")
    os.system(f"mkdir {RESULTS}/template/{case}/reg")
    os.system(f"mkdir {RESULTS}/template/{case}/reg/T2")
    os.system(f"mkdir {RESULTS}/template/{case}/reg/prostate")
    #img = sitk.ReadImage(f"{DATA}/{case}_T2_N4.nii.gz")

    #img.SetOrigin(REF_img.GetOrigin())
    #img.SetDirection(REF_img.GetDirection())
    #img.SetSpacing(REF_img.GetSpacing())
    """
    resampler = sitk.ResampleImageFilter()
    resampler.SetTransform(sitk.Transform(3, sitk.sitkIdentity))
    resampler.SetInterpolator(sitk.sitkLinear)
    resampler.SetReferenceImage(REF_img)

    img = resampler.Execute(img)
    """

    #sitk.WriteImage(img, f"{RESULTS}/template/{case}/Rsmp.nii.gz")
    
    img = nib.load(f"{DATA}/{case}_T2_N4.nii.gz")
    mask = nib.load(f"{DATA}/{case}_T2_N4_Prostate.nii.gz")
    img = conform(img, out_shape=shape, voxel_size=(1.0, 1.0, 1.0))
    mask = conform(mask, out_shape=shape, voxel_size=(1.0, 1.0, 1.0))    
    nib.save(img, f"{RESULTS}/template/{case}/Rsmp.nii.gz") 
    nib.save(mask, f"{RESULTS}/template/{case}/Rsmp_mask.nii.gz")

    os.system(f"elastix -f {REF} -m {RESULTS}/template/{case}/Rsmp.nii.gz -p deform.txt -out {RESULTS}/template/{case}/reg")
    os.system(f"transformix -in {RESULTS}/template/{case}/Rsmp.nii.gz -tp {RESULTS}/template/{case}/reg/TransformParameters.0.txt -out {RESULTS}/template/{case}/reg/T2")

    with open(f"{RESULTS}/template/{case}/reg/TransformParameters.0.txt", "r") as tp, open(f"{RESULTS}/template/{case}/reg/TransformParameters.0.mask.txt", "w") as tpm:
        for line in tp:
            if(line == "(FinalBSplineInterpolationOrder 3)\n"):
                tpm.write("(FinalBSplineInterpolationOrder 0)\n")
            else:
                tpm.write(line)

    os.system(f"transformix -in {RESULTS}/template/{case}/Rsmp_mask.nii.gz -tp {RESULTS}/template/{case}/reg/TransformParameters.0.mask.txt -out {RESULTS}/template/{case}/reg/prostate")
    #os.system(f"reg_aladin -ref {REF} -flo {RESULTS}/template/{case}/Rsmp.nii.gz -aff {RESULTS}/template/{case}/reg/aff.txt")
    #os.system(f"reg_f3d -ref {REF} -flo {RESULTS}/template/{case}/Rsmp.nii.gz -aff {RESULTS}/template/{case}/reg/aff.txt -cpp {RESULTS}/template/{case}/reg/CPP.nii.gz -res {RESULTS}/template/{case}/reg/nreg_result.nii.gz")

    images.append(sitk.GetArrayFromImage(sitk.ReadImage(f"{RESULTS}/template/{case}/reg/T2/result.mha")))
    masks.append(sitk.GetArrayFromImage(sitk.ReadImage(f"{RESULTS}/template/{case}/reg/prostate/result.mha")))
    

threads = []
for case in CASES:
    threads.append(threading.Thread(target=process, args=(case,)))

while len(threads) > 0:
    print(f"{len(threads)} cases left to run")
    n = len(CASES)
    to_run = threads[:n]
    [i.start() for i in to_run]
    [i.join() for i in to_run]
    threads = threads[n:]

print(f"{len(CASES)} original images, {len(images)} registered")

avg = np.median(images, axis=0)
avgMask = np.median(masks, axis = 0)

template = sitk.GetImageFromArray(avg)
template.SetOrigin(REF_img.GetOrigin())
template.SetSpacing(REF_img.GetSpacing())
template.SetDirection(REF_img.GetDirection())

templateMask = sitk.GetImageFromArray(avgMask)
templateMask.SetOrigin(REF_img.GetOrigin())
templateMask.SetSpacing(REF_img.GetSpacing())
templateMask.SetDirection(REF_img.GetDirection())
if rec:
    sitk.WriteImage(template, f"{RESULTS}/template/rec_template.nii.gz")
    sitk.WriteImage(templateMask, f"{RESULTS}/template/rec_template_mask.nii.gz")
else:
    sitk.WriteImage(template, f"{RESULTS}/template/norec_template.nii.gz")
    sitk.WriteImage(templateMask, f"{RESULTS}/template/norec_template_mask.nii.gz")
