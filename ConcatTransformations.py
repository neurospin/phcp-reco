#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 10:26:51 2024

@author: la272118
"""

import optparse, json, os
import ants
import nibabel as ni
import numpy as np
import SimpleITK as sitk
from scipy.ndimage import laplace

def print_message(message):
    print("==================================")
    print(message)
    print("==================================") 

"""Create Reference Init Space Files"""
def create_ReferenceInitSpace(filename, OutputDirectory):
    meta = ni.load(filename)
    xdim, ydim, zdim = meta.shape
    newarr = np.meshgrid(np.arange(xdim), np.arange(ydim), np.arange(zdim), indexing='ij')
    ni.save(ni.Nifti1Image(newarr[0], meta.affine, meta.header), os.path.join(OutputDirectory,'ReferenceInitSpace_x.nii.gz'))
    ni.save(ni.Nifti1Image(newarr[1], meta.affine, meta.header), os.path.join(OutputDirectory,'ReferenceInitSpace_y.nii.gz'))
    ni.save(ni.Nifti1Image(newarr[2], meta.affine, meta.header), os.path.join(OutputDirectory,'ReferenceInitSpace_z.nii.gz'))
    #Ajout 18/11/24
    ni.save(ni.Nifti1Image(np.ones(meta.shape), meta.affine, meta.header), os.path.join(OutputDirectory,'mask.nii.gz'))
    #
    return None

"""Load & Create Transform Files"""
def Create_transformationFilesList_from_JsonFilename(JsonFilename):
    with open(JsonFilename, 'r') as f :
        res = json.load(f)
    return res[list(res.keys())[0]]

#ADD 13/12/24

def compose_linear_transformations(MatrixTransformsFilenameList, nbr_trf, OutputDirectory):
    res = np.eye(4)
    restrictedList = MatrixTransformsFilenameList[:nbr_trf]
    # restrictedList.reverse()
    # print(type(restrictedList))
    center = np.asarray((0,0,0))
    for trans in restrictedList:
        information = sitk.ReadTransform(trans)
        
        if trans.split('.')[-1]=='mat':
            matrix = information.GetMatrix()
            trans = information.GetTranslation()
            if information.GetCenter() != list((0,0,0)):
                center = np.asarray(information.GetCenter())
            Euler =  create_matrix_matFile(matrix, trans)
        else :
            parameters = information.GetParameters()
            Euler = create_matrix_txtFile(parameters)
        res = res @ Euler
    OuputFilename = os.path.join(OutputDirectory, 'AffineTransform_nbr'+str(nbr_trf)+'.txt')
    with open(OuputFilename, 'w') as f:
        f.write('#Insight Transform File V1.0\n')
        f.write('#Transform 0\n')
        f.write('Transform: MatrixOffsetTransformBase_double_3_3\n')
        f.write('Parameters: ')
        f.write(" ".join(map(str, np.concatenate((np.asarray(res[:3,:3]).flatten(), np.asarray(res[:3,3].reshape((3,))) ) ) ) ))
        f.write("\n")
        f.write("FixedParameters: "+" ".join(map(str, center)))
    return OuputFilename

def identify_warpFiles(TransformationFilesList):
    warpfiles_indices = []
    for indice in range(len(TransformationFilesList)):
        if TransformationFilesList[indice].split(".")[0][-4:] == 'Warp' :
            warpfiles_indices.append(indice)
    return warpfiles_indices

def identify_NbrOfTransformations(warpfiles_indices):
    nbrOfLinearTransform = []
    for el in warpfiles_indices :
        nbrOfLinearTransform.append(el - warpfiles_indices.index(el))
    return nbrOfLinearTransform

def run_linear_transform_composer(TransformationFilesList, MatrixTransformsFilenameList, OutputDirectory):
    warpfiles_indices = identify_warpFiles(TransformationFilesList)
    nbrOfLinearTransform = identify_NbrOfTransformations(warpfiles_indices)
    ListLinearTrfComposed = []
    for nbr_trf in nbrOfLinearTransform:
        LinearTrfComposed = compose_linear_transformations(MatrixTransformsFilenameList, nbr_trf, OutputDirectory)
        ListLinearTrfComposed.append(LinearTrfComposed)
    return ListLinearTrfComposed

def make_complete_transformationlist(ListLinearTrfComposed, TransformationFilesList, warpfiles_indices):
    res = []
    for i in range(len(ListLinearTrfComposed)):
        res.append(ListLinearTrfComposed[i])
        res.append(TransformationFilesList[warpfiles_indices[i]])
        res.append('[%s, 1]' % ListLinearTrfComposed[i])
    res.reverse()
    return res



#END ADD 13/12/24


# def make_complete_transformationlist(TransformationFilesList):
#     res = []
#     for el in TransformationFilesList:
#         # print(el)
#         if el.split('.')[1] != 'nii':
#             # print('[%s, 1]' % el)
#             res.append('[%s, 1]' % el)
#     TransformationFilesList.reverse()
#     for el in TransformationFilesList:
#         res.append(el)
#     TransformationFilesList.reverse()
#     return res

"""Apply Transform Files"""
def antsApplyTransforms_fromtransformationList(fixedFilename, movingFilename, TransformationFilesList, interpolator, OutputFilename, verbose):
    command = "antsApplyTransforms -d 3"
    command += " -i %s" % movingFilename
    command += " -r %s" % fixedFilename    
    command += " -o %s" % OutputFilename
    command += " -n %s" % interpolator
    for el in TransformationFilesList:
        command += " -t %s" % el
    if verbose:
        command+= " -v 1"
    os.system(command)
    return None

def apply_transformationlist_To_ReferenceInitSpace_Files(InputFilename, OutputDirectory, TransformationFilesList):
    interpolator = "BSpline[5]"
    # interpolator = "Linear"
    movingFilename = os.path.join(OutputDirectory,'ReferenceInitSpace_x.nii.gz')
    OutputRefWarpedxFilename = os.path.join(OutputDirectory,'ReferenceInitSpaceWarped_x.nii.gz')
    antsApplyTransforms_fromtransformationList(InputFilename, movingFilename, TransformationFilesList, interpolator, OutputRefWarpedxFilename, verbose = True)
    
    movingFilename = os.path.join(OutputDirectory,'ReferenceInitSpace_y.nii.gz')
    OutputRefWarpedxFilename = os.path.join(OutputDirectory,'ReferenceInitSpaceWarped_y.nii.gz')
    antsApplyTransforms_fromtransformationList(InputFilename, movingFilename, TransformationFilesList, interpolator, OutputRefWarpedxFilename, verbose = True)

    movingFilename = os.path.join(OutputDirectory,'ReferenceInitSpace_z.nii.gz')
    OutputRefWarpedxFilename = os.path.join(OutputDirectory,'ReferenceInitSpaceWarped_z.nii.gz')
    antsApplyTransforms_fromtransformationList(InputFilename, movingFilename, TransformationFilesList, interpolator, OutputRefWarpedxFilename, verbose = True)
    
    #Ajout 18/11/24
    movingFilename = os.path.join(OutputDirectory,'mask.nii.gz')
    OutputRefWarpedxFilename = os.path.join(OutputDirectory,'MaskWarped.nii.gz')
    antsApplyTransforms_fromtransformationList(InputFilename, movingFilename, TransformationFilesList, "NearestNeighbor", OutputRefWarpedxFilename, verbose = True)
    #

    
    return None

"""Create Deformation Field"""
def change_header_to_deformation_fieldheader(ni_im):
    ni_im.header.set_intent(1007)
    ni_im.header['regular']=b'r'
    ni_im.header.set_xyzt_units(2)
    sform = ni_im.header.get_sform()
    ni_im.header.set_qform(sform, code = 1)
    ni_im.header.set_sform(sform, code = 1)
    return ni_im
    
def create_deformation_field_rapid(Filename, OutputDirectory):
    meta = ni.load(Filename)
    xdim, ydim, zdim = meta.shape
    resx, resy, resz = meta.header['pixdim'][1:4]
    metax = ni.load(OutputDirectory+'ReferenceInitSpaceWarped_x.nii.gz')
    metay = ni.load(OutputDirectory+'ReferenceInitSpaceWarped_y.nii.gz')
    metaz = ni.load(OutputDirectory+'ReferenceInitSpaceWarped_z.nii.gz')
    meta_mask = ni.load(OutputDirectory+'MaskWarped.nii.gz')

    newarr = np.meshgrid(np.arange(xdim), np.arange(ydim), np.arange(zdim), indexing='ij')
    arr_mask = meta_mask.get_fdata()

    # arr_Ref = np.concatenate( ( -resx*np.reshape(arr_mask*newarr[0], meta.shape+(1,1)), resy* np.reshape(arr_mask*newarr[1], meta.shape+(1,1)), -resz*np.reshape(arr_mask*newarr[2], meta.shape+(1,1))) , -1)
    # arr_SyN = np.concatenate( (-resx*np.reshape(arr_mask*metax.get_fdata(),(metax.shape+(1,1))), resy*np.reshape(arr_mask*metay.get_fdata(),(metay.shape+(1,1))), -resz*np.reshape(arr_mask*metaz.get_fdata(),(metaz.shape+(1,1)))) ,-1)
    
    # Fonctionne pour GIS DWI
    # arr_Ref = np.concatenate( ( resx*np.reshape(arr_mask*newarr[0], meta.shape+(1,1)), resy* np.reshape(arr_mask*newarr[1], meta.shape+(1,1)) ,-resz*np.reshape(arr_mask*newarr[2], meta.shape+(1,1)) ) , -1)
    # arr_SyN = np.concatenate( (resx*np.reshape(arr_mask*metax.get_fdata(),(metax.shape+(1,1))), resy*np.reshape(arr_mask*metay.get_fdata(),(metay.shape+(1,1))) , -resz*np.reshape(arr_mask*metaz.get_fdata(),(metaz.shape+(1,1)))) ,-1)

    
    
    # Fonctionne pour headerfliprawdata
    arr_Ref = np.concatenate( ( -resx*np.reshape(arr_mask*newarr[0], meta.shape+(1,1)),resz*np.reshape(arr_mask*newarr[2], meta.shape+(1,1)), resy* np.reshape(arr_mask*newarr[1], meta.shape+(1,1)) ) , -1)
    arr_SyN = np.concatenate( (-resx*np.reshape(arr_mask*metax.get_fdata(),(metax.shape+(1,1))), resz*np.reshape(arr_mask*metaz.get_fdata(),(metaz.shape+(1,1))), resy*np.reshape(arr_mask*metay.get_fdata(),(metay.shape+(1,1)))) ,-1)

    res = (arr_SyN-arr_Ref)
    
    ni_im = ni.Nifti1Image(res , metax.affine, metax.header) 
    ni_im = change_header_to_deformation_fieldheader(ni_im)    
    ni.save(ni_im, os.path.join(OutputDirectory,'Deformation_field_SyN.nii.gz'))
    return None

# def create_deformation_field_rapid(Filename, OutputDirectory):
#     meta = ni.load(Filename)
#     xdim, ydim, zdim = meta.shape
#     resx, resy, resz = meta.header['pixdim'][1:4]
#     metax = ni.load(OutputDirectory+'ReferenceInitSpaceWarped_x.nii.gz')
#     metay = ni.load(OutputDirectory+'ReferenceInitSpaceWarped_y.nii.gz')
#     metaz = ni.load(OutputDirectory+'ReferenceInitSpaceWarped_z.nii.gz')
#     meta_mask = ni.load(OutputDirectory+'MaskWarped.nii.gz')

#     newarr = np.meshgrid(np.arange(xdim), np.arange(ydim), np.arange(zdim), indexing='ij')
#     arr_mask = meta_mask.get_fdata()

#     arr_Ref = np.concatenate( ( resx*np.reshape(arr_mask*newarr[0], meta.shape+(1,1)), resy* np.reshape(arr_mask*newarr[1], meta.shape+(1,1)), -resz*np.reshape(arr_mask*newarr[2], meta.shape+(1,1))) , -1)
#     arr_SyN = np.concatenate( (resx*np.reshape(arr_mask*metax.get_fdata(),(metax.shape+(1,1))), resy*np.reshape(arr_mask*metay.get_fdata(),(metay.shape+(1,1))), -resz*np.reshape(arr_mask*metaz.get_fdata(),(metaz.shape+(1,1)))) ,-1)

#     ni_im = ni.Nifti1Image( (arr_SyN-arr_Ref), metax.affine, metax.header) 
#     ni_im = change_header_to_deformation_fieldheader(ni_im)    
#     ni.save(ni_im, os.path.join(OutputDirectory,'Deformation_field_SyN.nii.gz'))
#     return None

""" Create Jacobian files """

def create_jacobian_files(Filename, OutputDirectory):
    ants_ref = ants.image_read(Filename)
    jacobian = ants.create_jacobian_determinant_image(ants_ref,
                                                      os.path.join(OutputDirectory,'Deformation_field_SyN.nii.gz'),
                                                      do_log = False)
    ants.image_write(jacobian, os.path.join(OutputDirectory,'Jacobian_Deformation_field_SyN.nii.gz'))
    jacobian = ants.create_jacobian_determinant_image(ants_ref,
                                                      os.path.join(OutputDirectory,'Deformation_field_SyN.nii.gz'),
                                                      do_log = True)
    ants.image_write(jacobian, os.path.join(OutputDirectory,'JacobianLog_Deformation_field_SyN.nii.gz'))
    return None

"""Create TotalAffineTransform File"""
def ExtractMatrixTransforms(TransformationFilesList):
    MatrixTransformsFilenameList = []
    for el in TransformationFilesList:
        if el.split('.')[1] != 'nii':
            MatrixTransformsFilenameList.append(el)
    return MatrixTransformsFilenameList

def create_matrix_matFile(matrix, trans):
    matrix, trans = np.asarray(matrix).reshape((3,3)), np.asarray(trans).reshape((3,1))
    ligne = np.asarray([0,0,0,1]).reshape((1,4))
    hst = np.hstack((matrix, trans.reshape((3,1))))
    return np.vstack((hst, ligne))

def create_matrix_txtFile(parameters):
    matrix, trans = np.asarray(parameters)[:9].reshape((3,3)), np.asarray(parameters)[9:].reshape((3,1))
    ligne = np.asarray([0,0,0,1]).reshape((1,4))
    hst = np.hstack((matrix, trans.reshape((3,1))))
    return np.vstack((hst, ligne))

def compose_transformations(MatrixTransformsFilenameList, OutputDirectory):
    res = np.eye(4)
    center = np.asarray((0,0,0))
    for trans in MatrixTransformsFilenameList:
        information = sitk.ReadTransform(trans)

        if trans.split('.')[-1]=='mat':
            matrix = information.GetMatrix()
            trans = information.GetTranslation()
            if information.GetCenter() != list((0,0,0)):
                center = np.asarray(information.GetCenter())
            Euler = create_matrix_matFile(matrix, trans)
        else :
            parameters = information.GetParameters()
            Euler = create_matrix_txtFile(parameters)

        res = res@Euler
    T_center, T_neg_center = np.eye(4), np.eye(4)  #Need to take into account inside the affine matrix the Center of rotation for the next script (ConvertAntsToGkgFormat.py) 
    T_center[:3, 3], T_neg_center[:3, 3] = center, -center
    res = T_center @ res @ T_neg_center
    OuputFilename = os.path.join(OutputDirectory, 'TotalAffineTransform.txt')
    with open(OuputFilename, 'w') as f:
        f.write('#Insight Transform File V1.0\n')
        f.write('#Transform 0\n')
        f.write('Transform: MatrixOffsetTransformBase_double_3_3\n')
        f.write('Parameters: ')
        f.write(" ".join(map(str, np.concatenate((np.asarray(res[:3,:3]).flatten(), np.asarray(res[:3,3].reshape((3,))) ) ) ) ))
        f.write("\n")
        f.write("FixedParameters: 0 0 0") #+" ".join(map(str, center)))
        
    return None


### Ajout 18/11/2024 début

def regularize_laplacian_smoothing(deformationField_array, iteration, Lambda, component):
    arr = deformationField_array[:,:,:,0,component]
    for _ in range(iteration):
        arr += Lambda* laplace(arr) 
    return arr, laplace(arr) 

def apply_laplacian_smoothing(OutputDirectory):
    meta_deformationField = ni.load(os.path.join(OutputDirectory,'Deformation_field_SyN.nii.gz'))
    arr_deformationField = meta_deformationField.get_fdata()
    
    arrsmoothed_deformationfield_component_x, Laplacian_deformationfield_component_x = regularize_laplacian_smoothing(arr_deformationField, 3, 0.2, 0)
    ni_im = ni.Nifti1Image(Laplacian_deformationfield_component_x, meta_deformationField.affine)
    ni.save(ni_im, os.path.join(OutputDirectory,'Laplacian_component_x.nii.gz'))

    arrsmoothed_deformationfield_component_y, Laplacian_deformationfield_component_y = regularize_laplacian_smoothing(arr_deformationField, 3, 0.2, 1)
    ni_im = ni.Nifti1Image(Laplacian_deformationfield_component_y, meta_deformationField.affine)
    ni.save(ni_im, os.path.join(OutputDirectory,'Laplacian_component_y.nii.gz'))

    arrsmoothed_deformationfield_component_z, Laplacian_deformationfield_component_z = regularize_laplacian_smoothing(arr_deformationField, 3, 0.2, 2)
    ni_im = ni.Nifti1Image(Laplacian_deformationfield_component_z, meta_deformationField.affine)
    ni.save(ni_im, os.path.join(OutputDirectory,'Laplacian_component_z.nii.gz'))
    
    
    SmoothedDeformationField = np.concatenate( (np.reshape(arrsmoothed_deformationfield_component_x, meta_deformationField.shape[:4] + (1,)), 
                                                np.reshape(arrsmoothed_deformationfield_component_y, meta_deformationField.shape[:4] + (1,)),
                                                np.reshape(arrsmoothed_deformationfield_component_z, meta_deformationField.shape[:4] + (1,)) ), axis = -1 )
    
    ni_im = ni.Nifti1Image(SmoothedDeformationField, meta_deformationField.affine, meta_deformationField.header)
    ni.save(ni_im, os.path.join(OutputDirectory,'Deformation_field_SyN_Smoothed.nii.gz'))
    
    return None

### fin Ajout 18/11/2024. 


        
def RunPipeline(InputFilename, JsonFilename, OutputDirectory):
    print_message("Create Reference Init Space Files")
    create_ReferenceInitSpace(InputFilename, OutputDirectory)
    
    print_message("Load & Create Transform Files")
    TransformationFilesList = Create_transformationFilesList_from_JsonFilename(JsonFilename)
    MatrixTransformsFilenameList = ExtractMatrixTransforms(TransformationFilesList) 
    ListLinearTrfComposed = run_linear_transform_composer(TransformationFilesList, MatrixTransformsFilenameList, OutputDirectory)
    warpfiles_indices =identify_warpFiles(TransformationFilesList)
    NewTransformationFilesList = make_complete_transformationlist(ListLinearTrfComposed, TransformationFilesList, warpfiles_indices)
    
    print_message("Create TotalAffineTransform File")
    compose_transformations(MatrixTransformsFilenameList, OutputDirectory)
    
    # NewTransformationFilesList = make_complete_transformationlist(TransformationFilesList)
    
    print_message("Apply Transform Files")
    apply_transformationlist_To_ReferenceInitSpace_Files(InputFilename, OutputDirectory, NewTransformationFilesList)
    
    print_message("Create Deformation Field")
    create_deformation_field_rapid(InputFilename, OutputDirectory)
    
    print_message("Create Jacobian Determinant files")
    create_jacobian_files(InputFilename, OutputDirectory)
    

    print_message("Deformation Field - Laplacian Smoothing")
    apply_laplacian_smoothing(OutputDirectory)
    return None



parser = optparse.OptionParser()
parser.add_option( '-i', '--input',
                   dest = 'InputFilename',
                   help = 'input filename' )
parser.add_option( '-j', '--json',
                   dest = 'JsonFilename',
                   help = 'Json filename containing transformation from the init to the final space' )
parser.add_option( '-o', '--output',
                   dest = 'OutputDirectory',
                   help = 'output directory' )

( options, args ) = parser.parse_args()


################################################################################
# 1) Value Extractor 
################################################################################

RunPipeline( options.InputFilename,
        options.JsonFilename,
        options.OutputDirectory)



# def compose_transformations2(JSON_filename, Output):
#     Transforms = dict()
#     with open(JSON_filename, 'r') as f:
#         Transforms = json.load(f)
    
#     res = np.eye(4)
#     for trans in Transforms['TransformationFiles']:
#         print(trans)
#         information = sitk.ReadTransform(trans)

#         if trans.split('.')[-1]=='mat':
#             matrix = information.GetMatrix()
#             trans = information.GetTranslation()
#             Euler = create_matrix_matFile(matrix, trans)
#         else :
#             parameters = information.GetParameters()
#             Euler = create_matrix_txtFile(parameters)

#         res = res@Euler
#     with open(Output, 'w') as f:
#         f.write('#Insight Transform File V1.0\n')
#         f.write('#Transform 0\n')
#         f.write('Transform: MatrixOffsetTransformBase_double_3_3\n')
#         f.write('Parameters: ')
#         f.write(" ".join(map(str, np.concatenate((np.asarray(res[:3,:3]).flatten(), np.asarray(res[:3,3].reshape((3,))) ) ) ) ))
#         f.write("\n")
#         f.write("FixedParameters: 0 0 0")
        
#     # data = {'AffineTransform_float_3_3': res[:3,:4].reshape((12,1)), 'fixed':np.asarray([0,0,0]).reshape((3,1))}
#     # scipy.io.savemat(Output, data)
#     return None