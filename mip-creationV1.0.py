# To do 
# Faster registraion
# Skull stripping.. Why I am not using erosion... filter....
# More Modular approach.....
# Licensing.....

import SimpleITK as sitk
import numpy as np
import os,fnmatch
import nibabel as nb
import imageio
from dipy.align.reslice import reslice
from skimage.morphology import closing
from skimage import util
from skimage.morphology import disk,dilation

cc=os.getcwd()
selem = disk(6)
delem=disk(3)

for each in os.listdir('./'):
    print('===================================================================')
    print("Working on:"+each)
    Image=sitk.ReadImage(each)
    Mask=sitk.ReadImage('/Users/gaganjsharma/Downloads/20181217T154339.430000/mask.nii.gz')
    Image_Array=sitk.GetArrayFromImage(Image)
    Image_Mask=sitk.GetArrayFromImage(Mask)
    newMask=np.zeros([Image_Mask.shape[0],Image_Mask.shape[1],Image_Mask.shape[2]])
    for ieach in range(Image_Mask.shape[0]):
        newMask[ieach,:,:]=closing(Image_Mask[ieach,:,:],selem)
        
    newMaskImage=sitk.GetImageFromArray(newMask)
    newMaskImage.CopyInformation(Mask)
    #sitk.WriteImage(newMaskImage,'new_dbgmaskout.mha')
    new_ImageA=Image_Array
    new_ImageA[newMask==0]=0
    Axial_MIP=np.zeros([1,Image_Array.shape[1],Image_Array.shape[2]])
    Axial_MIP[0,:,:] = np.amax(new_ImageA[:,:,:],0)
    Axial_MIP=util.invert(Axial_MIP)
    Axial_MIP_Image=sitk.GetImageFromArray(Axial_MIP)
    Axial_MIP_Image.SetSpacing(Image.GetSpacing())
    Axial_MIP_Image.SetDirection(Image.GetDirection())
    Axial_MIP_Image.SetOrigin(Image.GetOrigin())
    sitk.WriteImage(Axial_MIP_Image,'Axial_MIP_'+str(each))


    Coronal_MIP=np.zeros([Image_Array.shape[0],1,Image_Array.shape[2]])
    Coronal_MIP[:,0,:] = np.amax(new_ImageA[:,:,:],1)
    Coronal_MIP=util.invert(Coronal_MIP)
    Coronal_MIP=np.reshape(Coronal_MIP,(Coronal_MIP.shape[1],Coronal_MIP.shape[0],Coronal_MIP.shape[2]))
    Coronal_MIP_Image=sitk.GetImageFromArray(Coronal_MIP)
    Coronal_MIP_Image.SetDirection(Image.GetDirection())
    Coronal_MIP_Image.SetOrigin(Image.GetOrigin())
    Coronal_MIP_Image.SetSpacing([Image.GetSpacing()[0],Image.GetSpacing()[2],Image.GetSpacing()[1]])
    sitk.WriteImage(Coronal_MIP_Image,'Coronal_MIP_'+each)
    #====================================================================================
    corimage=nb.load('Coronal_MIP_'+each)
    zoom_image=nb.load(each)
    cordata = corimage.get_fdata()
    zooms = corimage.header.get_zooms()[:3]
    print(zooms)
    new_zooms=[zoom_image.header.get_zooms()[0],zoom_image.header.get_zooms()[1],zoom_image.header.get_zooms()[0]]
    print(new_zooms)
    data2, affine2 = reslice(cordata, corimage.affine, zooms, new_zooms)
    array_img=nb.Nifti1Image(data2,affine2)
    nb.save(array_img,'Resliced_Coronal_MIP_'+each)
    #====================================================================================
    Sagittal_MIP_one=np.zeros([Image_Array.shape[0],Image_Array.shape[1],1])
    Sagittal_MIP_two=np.zeros([Image_Array.shape[0],Image_Array.shape[1],1])

    Sagittal_MIP_one[:,:,0]=np.amax(new_ImageA[:,:,0:int(new_ImageA.shape[2]/2-1)],2)
    Sagittal_MIP_one=util.invert(Sagittal_MIP_one)
    Sagittal_MIP_two[:,:,0]=np.amax(new_ImageA[:,:,int(new_ImageA.shape[2]/2):new_ImageA.shape[2]],2)
    Sagittal_MIP_two=util.invert(Sagittal_MIP_two)

    Sagittal_MIP_one_Image=sitk.GetImageFromArray(Sagittal_MIP_one)
    Sagittal_MIP_two_Image=sitk.GetImageFromArray(Sagittal_MIP_two)


    Sagittal_MIP_one_Image.SetSpacing(Image.GetSpacing())
    Sagittal_MIP_one_Image.SetDirection(Image.GetDirection())
    Sagittal_MIP_one_Image.SetOrigin(Image.GetOrigin())
    sitk.WriteImage(Sagittal_MIP_one_Image,'Sagittal_MIP1_'+each)

    Sagittal_MIP_two_Image.SetSpacing(Image.GetSpacing())
    Sagittal_MIP_two_Image.SetDirection(Image.GetDirection())
    Sagittal_MIP_two_Image.SetOrigin(Image.GetOrigin())
    sitk.WriteImage(Sagittal_MIP_two_Image,'Sagittal_MIP2_'+each)
    #======================================================================================
    sag1image=nb.load('Sagittal_MIP1_'+each)
    sag1data=sag1image.get_fdata()
    s1zooms = sag1image.header.get_zooms()[:3]
    data2s1, affine2s1 = reslice(sag1data, sag1image.affine,s1zooms, new_zooms)
    array_img=nb.Nifti1Image(data2s1,affine2s1)
    nb.save(array_img,'Resliced_SAG1_MIP1_'+each)
    #=======================================================================================
    sag2image=nb.load('Sagittal_MIP2_'+each)
    sag2data=sag2image.get_fdata()
    s2zooms = sag2image.header.get_zooms()[:3]
    data2s2, affine2s2 = reslice(sag2data, sag2image.affine,s1zooms, new_zooms)
    array_img=nb.Nifti1Image(data2s2,affine2s2)
    nb.save(array_img,'Resliced_SAG2_MIP2_'+each)

    #======================================================================================

    #Lets create Montage Nifti
    A_array=sitk.GetArrayFromImage(sitk.ReadImage('Axial_MIP_'+each))
    C_array=sitk.GetArrayFromImage(sitk.ReadImage('Resliced_Coronal_MIP_'+each))
    S1_array=sitk.GetArrayFromImage(sitk.ReadImage('Resliced_SAG1_MIP1_'+each))
    S2_array=sitk.GetArrayFromImage(sitk.ReadImage('Resliced_SAG2_MIP2_'+each))


    new_array=np.zeros([np.squeeze(A_array).shape[0]+np.squeeze(C_array).shape[0]+20,np.squeeze(A_array).shape[1]+np.squeeze(C_array).shape[1]+20])
    new_array[0:np.squeeze(A_array).shape[0],0:np.squeeze(A_array).shape[1]]=np.squeeze(A_array)
    new_array[np.squeeze(A_array).shape[0]+20:new_array.shape[0],0:np.squeeze(A_array).shape[1]]=np.squeeze(C_array)

    new_array[50:np.squeeze(S1_array).shape[0]+50,np.squeeze(A_array).shape[1]:new_array.shape[1]-20]=np.squeeze(S1_array)
    new_array[np.squeeze(A_array).shape[0]+20:new_array.shape[0],np.squeeze(S2_array).shape[1]:new_array.shape[1]-20]=np.squeeze(S2_array)
    new_array=np.reshape(new_array,[1,new_array.shape[0],new_array.shape[1]])
    new_array_Image=sitk.GetImageFromArray(new_array)
    new_array_Image.SetSpacing(Image.GetSpacing())
    new_array_Image.SetDirection(Image.GetDirection())
    new_array_Image.SetOrigin(Image.GetOrigin())
    sitk.WriteImage(new_array_Image,"new_MIP_consolidated"+each)
    gifarray=np.squeeze(new_array)
    gifarray=np.flipud(gifarray)
    imageio.imwrite("new_MIP_consolidated"+os.path.splitext(each)[0]+".png",gifarray)
    #imageio.imwrite("new_MIP_consolidated"+os.path.splitext(image_name)[0]+".png",gifarray)
    os.chdir(cc)
