#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 15:59:26 2025

@author: la272118
"""

import subprocess, json, optparse, os
import numpy as np
import nibabel as ni
import nibabel.processing as proc
import glob
import scipy.ndimage as nd

QMAP_VALUES = ['QT2', 'QT2star']
# QMAP_VALUES = ['QT1', 'QT2', 'QT2star', 'proton-density']
# CM_QMAP = ['CM_T1', 'CM_T2', 'CM_T2star', 'CM_T2star']
CM_QMAP = ['CM_T2', 'CM_T2star']


""" Miscellaneous algorithms """
def print_message(message):
    print("==================================")
    print(message)
    print("==================================") 


def ConversionGisToNifti(InputGisFilename, OutputNiftiFilename):
    subprocess.run(["singularity","exec", "--bind", "/neurospin:/neurospin:rw",
                    "/neurospin/phcp/code/gkg/2022-12-20_gkg/2022-12-20_gkg.sif",
                    "GkgExecuteCommand", "Gis2NiftiConverter", "-i", InputGisFilename, 
                    "-o", OutputNiftiFilename, "-verbose"])

      
def GetMaskGkg(InputGisFilename, OutputGisFilename):
    subprocess.run(["singularity","exec", "--bind", "/neurospin:/neurospin:rw",
                    "/neurospin/phcp/code/gkg/2022-12-20_gkg/2022-12-20_gkg.sif",
                    "GkgExecuteCommand", "GetMask", "-i", InputGisFilename, 
                    "-o", OutputGisFilename, "-a", "1", "-format", "Gis", "-verbose"])

# Utilisation de cette fonction sur GIS file geomean pour avoir un masque propre
def create_ArrayFromGisFilename_16bit(GisFilename):
    DimFilename = GisFilename[:-3]+'dim'
    TxtInDimFile = open(DimFilename, 'r').read()
    TxtFirstLine = TxtInDimFile.split('\n')[0]
    shape = tuple(np.int16(TxtFirstLine.split(' ')))[:-1]
    arr = np.fromfile(GisFilename, np.float16())
    arr = np.reshape(arr, shape, order = 'F')
    return arr

def create_ArrayFromGisFilename_32bit(GisFilename):
    DimFilename = GisFilename[:-3]+'dim'
    TxtInDimFile = open(DimFilename, 'r').read()
    TxtFirstLine = TxtInDimFile.split('\n')[0]
    shape = tuple(np.int16(TxtFirstLine.split(' ')))[:-1]
    arr = np.fromfile(GisFilename, np.float32())
    arr = np.reshape(arr, shape, order = 'F')
    return arr

def arrange_ArrayToFlipHeaderFormat(arr):
    newarr = np.swapaxes(arr, 1, 2)
    newarr = np.flip(newarr, 1)
    newarr = np.flip(newarr, 0)
    return newarr

def load_jsonfile(JsonFilename):
    with open(JsonFilename, 'r') as f :
        dictionary = json.load(f)
    return dictionary

def send_to_RefSpace(InputFilename, OutputFilename, JsonFilename, Interpolator):
    dictionary = load_jsonfile(JsonFilename)
    RefSpaceFilename = dictionary['RefSpaceFilename']
    tparams = dictionary['tparams']
    
    str_command = ["antsApplyTransforms",
                    "-d", "3",
                    "-i", InputFilename,
                    "-r", RefSpaceFilename,
                    "-o", OutputFilename,
                    "-n", Interpolator]
    for i in range(len(tparams)):
        tparams.insert(2*i, "-t")
    str_command = str_command + tparams + ["-v", "1"]
    subprocess.run(str_command)
    return None

def send_to_RefSpace_Mask(InputFilename, OutputFilename, JsonFilename, Interpolator):
    dictionary = load_jsonfile(JsonFilename)
    RefSpaceFilename = dictionary['RefSpaceFilename']
    tparams = dictionary['tparams']
    InitSpaceFilename= dictionary['InitSpaceFilename']
    IntermediateFilename = InputFilename.split('.')[0] +'_inter.nii.gz'
    
    str_command = ["antsApplyTransforms",
                    "-d", "3",
                    "-i", InputFilename,
                    "-r", InitSpaceFilename,
                    "-o", IntermediateFilename,
                    "-n", Interpolator]
    
    tparams_1 = tparams[-3:]
    for i in range(len(tparams_1)):
        tparams_1.insert(2*i, "-t")
    str_command_1 = str_command + tparams_1 + ["-v", "1"]
    subprocess.run(str_command_1)
    
    str_command = ["antsApplyTransforms",
                    "-d", "3",
                    "-i", IntermediateFilename,
                    "-r", RefSpaceFilename,
                    "-o", OutputFilename,
                    "-n", Interpolator]
    
    tparams_2 = tparams[:-3]
    for i in range(len(tparams_2)):
        tparams_2.insert(2*i, "-t")
    str_command_2 = str_command + tparams_2 + ["-v", "1"]
    subprocess.run(str_command_2)
    return None

def create_mask_FOV(inputfilename, outputfilename):
    meta_input = ni.load(inputfilename)
    res = np.ones(meta_input.shape)
    ni_im = ni.Nifti1Image(res, meta_input.affine, meta_input.header)
    ni.save(ni_im, outputfilename)

def sigmoid(x,k, x0):
    return 1/(1+np.exp(-k*(x-x0)))

def create_mask_from_geomean(ArchitecturePath, InputGisFilename, OutputMaskGisFilename, OutputFilename):
    
    if not(os.path.exists(OutputMaskGisFilename)):
        GetMaskGkg(InputGisFilename, OutputMaskGisFilename)
    OutputMaskNiftiFilename = OutputMaskGisFilename[:-3]+ 'nii.gz'
    
    if not(os.path.exists(OutputMaskNiftiFilename)):
        ConversionGisToNifti(OutputMaskGisFilename, OutputMaskNiftiFilename)
    
    JsonFilename = os.path.join(ArchitecturePath,
                                    "SendToRefSpace_DWI_" +\
                                        InputGisFilename.split('.')[0].split('_')[-1]+\
                                            ".json")
    if not(os.path.exists(OutputFilename)):
        
        send_to_RefSpace_Mask(OutputMaskNiftiFilename, OutputFilename, JsonFilename, 'NearestNeighbor')
    
    Mask_FOV_Filename = OutputMaskNiftiFilename.split('.')[0]+'_FOV.nii.gz'
    Mask_FOV_refSpace_Filename = OutputFilename.split('Mask')[0]+ Mask_FOV_Filename.split('/')[-1]
    weight_filename = Mask_FOV_refSpace_Filename.split('.')[0] + "_weight.nii.gz"

    if not(os.path.exists(weight_filename)):
        create_mask_FOV(OutputMaskNiftiFilename, Mask_FOV_Filename)
        send_to_RefSpace_Mask(Mask_FOV_Filename, Mask_FOV_refSpace_Filename, JsonFilename, 'NearestNeighbor')
        meta = ni.load(Mask_FOV_refSpace_Filename)
        arr = meta.get_fdata()
        newarr = nd.distance_transform_edt(arr, meta.header['pixdim'][1])
        newarr = sigmoid(newarr, 1, 6.5)
        ni_im = ni.Nifti1Image(newarr, meta.affine, meta.header)
        ni.save(ni_im, weight_filename)
    
    return None


def create_weight(CM_filename1, CM_filename2):
    CM_1 = ni.load(CM_filename1)
    CM_2 = ni.load(CM_filename2)
    
    CM_arr1 = proc.smooth_image(CM_1, proc.sigma2fwhm(0.2)).get_fdata()
    CM_arr2 = proc.smooth_image(CM_2, proc.sigma2fwhm(0.2)).get_fdata()

    return sigmoid(np.log(CM_arr1/CM_arr2), 4,0)

def create_mask_overlap_block(Mask_FilenameList):
    overlap_mask = np.zeros(ni.load(Mask_FilenameList[0]).shape, dtype = np.int8())
    for i in range(len(Mask_FilenameList)-1):
        overlap_mask |= (np.int8(ni.load(Mask_FilenameList[i]).get_fdata()) &\
                np.int8(ni.load(Mask_FilenameList[i+1]).get_fdata()))
    return overlap_mask

def fill_Qvalues_nonOverlapping(Mask_FilenameList, Qmap_FilenameList, overlap_mask):
    Qvalues_nonOverlapping = np.zeros(ni.load(Mask_FilenameList[0]).shape)
    for mask, Qmap in zip(Mask_FilenameList, Qmap_FilenameList):
        Qvalues_nonOverlapping += (ni.load(mask).get_fdata() * ni.load(Qmap).get_fdata())*np.logical_not(overlap_mask)
    return Qvalues_nonOverlapping

def fill_Qvalues_Overlapping(CM_FilenameList, Qmap_FilenameList,overlap_mask):
    Qvalues_Overlapping = np.zeros(ni.load(CM_FilenameList[0]).shape)
    for i in range(len(CM_FilenameList)-1):
        weight_2 = create_weight(CM_FilenameList[i], CM_FilenameList[i+1])
        weight_1 = 1-weight_2
        Qvalues_Overlapping += (weight_1 * ni.load(Qmap_FilenameList[i]).get_fdata() +\
                                weight_2 * ni.load(Qmap_FilenameList[i+1]).get_fdata())*overlap_mask
    return Qvalues_Overlapping

def reconstructionv2(CM_FilenameList, Qmap_FilenameList, Mask_FilenameList):
    res = (ni.load(Mask_FilenameList[0]).get_fdata() * ni.load(Qmap_FilenameList[0]).get_fdata())
    for i in range(len(CM_FilenameList)-1):
    
        meta_mask = ni.load(Mask_FilenameList[i+1])
        arr_mask = meta_mask.get_fdata()
        meta_mask_i = ni.load(Mask_FilenameList[i])
        arr_mask_i = meta_mask_i.get_fdata()
            
        overlap_mask = arr_mask*arr_mask_i
        res = (ni.load(Mask_FilenameList[i+1]).get_fdata() * ni.load(Qmap_FilenameList[i+1]).get_fdata() + res)*np.logical_not(overlap_mask) #fill_Qvalues_nonOverlapping2(Mask_FilenameList, Qmap_FilenameList, res, overlap_mask, i)
        res += np.nan_to_num(fill_Qvalues_Overlapping2(CM_FilenameList, Qmap_FilenameList, Mask_FilenameList,overlap_mask, i))
    return res
        
    
def fill_Qvalues_Overlapping2(CM_FilenameList, Qmap_FilenameList, Mask_FilenameList,overlap_mask, i):
    Qvalues_Overlapping = np.zeros(ni.load(CM_FilenameList[0]).shape)
    mask_weight_filename = Mask_FilenameList[i+1].split('.')[0] + '_FOV_weight.nii.gz'
    meta_mask_weight = ni.load(mask_weight_filename)
    arr_mask_weight = meta_mask_weight.get_fdata()
    mask_weight_filename_i = Mask_FilenameList[i].split('.')[0] + '_FOV_weight.nii.gz'
    meta_mask_weight_i = ni.load(mask_weight_filename_i)
    arr_mask_weight_i = meta_mask_weight_i.get_fdata()

    
    weight_2 = create_weight(CM_FilenameList[i], CM_FilenameList[i+1]) * arr_mask_weight
    weight_1 = (1-weight_2)*arr_mask_weight_i
    weight_2 = 1- weight_1
    Qvalues_Overlapping += (weight_1 * ni.load(Qmap_FilenameList[i]).get_fdata() +\
                            weight_2 * ni.load(Qmap_FilenameList[i+1]).get_fdata())*overlap_mask
    return Qvalues_Overlapping
    
def fill_Qvalues_Overlappingv2(CM_FilenameList, Qmap_FilenameList, Mask_FilenameList, overlap_mask):
    Qvalues_Overlapping = np.zeros(ni.load(CM_FilenameList[0]).shape)
    for i in range(len(CM_FilenameList)-1):
        mask_weight_filename = Mask_FilenameList[i+1].split('.')[0] + '_FOV_weight.nii.gz'
        meta_mask_weight = ni.load(mask_weight_filename)
        arr_mask_weight = meta_mask_weight.get_fdata()
        mask_weight_filename_i = Mask_FilenameList[i].split('.')[0] + '_FOV_weight.nii.gz'
        meta_mask_weight_i = ni.load(mask_weight_filename_i)
        arr_mask_weight_i = meta_mask_weight_i.get_fdata()
        
        weight_2 = create_weight(CM_FilenameList[i], CM_FilenameList[i+1]) * arr_mask_weight
        weight_1 = (1-weight_2)*arr_mask_weight_i
        weight_2 = 1- weight_1
        Qvalues_Overlapping += (weight_1 * ni.load(Qmap_FilenameList[i]).get_fdata() +\
                                weight_2 * ni.load(Qmap_FilenameList[i+1]).get_fdata())*overlap_mask
    return Qvalues_Overlapping

def merge_QMRI_FOVs(ArchitecturePath, JsonFilename, QmapToMerge):
    FOVmaterialFilename = load_jsonfile(JsonFilename)
        
    Mask_FilenameList = FOVmaterialFilename['mask']
    Qmap_FilenameList = FOVmaterialFilename[QmapToMerge]
    CM_FilenameList = FOVmaterialFilename[CM_QMAP[QMAP_VALUES.index(QmapToMerge)]]
    
    meta = ni.load(Qmap_FilenameList[0])
    
    Qmap_bloc = reconstructionv2(CM_FilenameList, Qmap_FilenameList, Mask_FilenameList)
    Qmap_bloc = np.nan_to_num(Qmap_bloc)
    
    ni_im = ni.Nifti1Image(Qmap_bloc, meta.affine, meta.header)
    return ni_im

def merge_QMRI_blocks(ArchitecturePath, value, nbrBlocks):
    BlocksPath = os.path.join(ArchitecturePath, "04-Blocks")
    for i in range(nbrBlocks):
        blockFilename =  os.path.join(BlocksPath, 'bloc_'+str(i+1)+'_'+value+'.nii.gz' )
        if i < 1:
            meta_first_block = ni.load(blockFilename)
            res = meta_first_block.get_fdata()
        else :
            res = np.maximum(res, ni.load(blockFilename).get_fdata())
    return ni.Nifti1Image(res, meta_first_block.affine, meta_first_block.header)

def correct_proton_density_intensities(FOVmaterialFilename):
    Mask_FilenameList = FOVmaterialFilename['mask']
    ProtonDensity_FilenameList = FOVmaterialFilename['proton-density']

    for i in range(len(Mask_FilenameList)-1):
        ProtonDensity_corrected_filename = ProtonDensity_FilenameList[i+1].split('.')[0]+'_corr.nii.gz'
        if not(os.path.exists(ProtonDensity_corrected_filename)):
            Overlap = ni.load(Mask_FilenameList[i]).get_fdata() * ni.load(Mask_FilenameList[i+1]).get_fdata()
            ProtonDensityTocorrect_meta = ni.load(ProtonDensity_FilenameList[i+1])
            ProtonDensityTocorrect_arr = ProtonDensityTocorrect_meta.get_fdata()
            ProtonDensityREF_meta = ni.load(ProtonDensity_FilenameList[i])
            ProtonDensityREF_arr = ProtonDensityREF_meta.get_fdata()
            Ratio_correction = Overlap*(ProtonDensityREF_arr/ProtonDensityTocorrect_arr)
            Ratio_correction = Ratio_correction[~(np.isinf(Ratio_correction)|np.isnan(Ratio_correction)|(Ratio_correction==0))]
            Ratio_correction = np.mean(Ratio_correction)
            print(Ratio_correction)
            ProtonDensity_corrected_arr = Ratio_correction*ProtonDensityTocorrect_arr
            ProtonDensity_corrected_meta = ni.Nifti1Image(ProtonDensity_corrected_arr,
                                                          ProtonDensityTocorrect_meta.affine,
                                                          ProtonDensityTocorrect_meta.header)
            ni.save(ProtonDensity_corrected_meta, ProtonDensity_corrected_filename)
    return None 

def run_DWI_all_blocks(ArchitecturePath, nbrBlocks):
    
    for i in range(nbrBlocks):
        JsonFilename = os.path.join(ArchitecturePath, 'bloc_'+str(i+1)+'_intensity.json' )
        FOVmaterialFilename = load_jsonfile(JsonFilename)
        correct_DWI_intensities(FOVmaterialFilename)
    return None

def correct_DWI_intensities(FOVmaterialFilename):
    Mask_FilenameList = FOVmaterialFilename['mask']
    ProtonDensity_FilenameList = FOVmaterialFilename['DWI_1500']

    for i in range(len(Mask_FilenameList)-1):
        ProtonDensity_corrected_filename = ProtonDensity_FilenameList[i+1].split('.')[0]+'_corr.nii.gz'
        if not(os.path.exists(ProtonDensity_corrected_filename)):
            Overlap = ni.load(Mask_FilenameList[i]).get_fdata() * ni.load(Mask_FilenameList[i+1]).get_fdata()
            ProtonDensityTocorrect_meta = ni.load(ProtonDensity_FilenameList[i+1])
            ProtonDensityTocorrect_arr = ProtonDensityTocorrect_meta.get_fdata()
            ProtonDensityREF_meta = ni.load(ProtonDensity_FilenameList[i])
            ProtonDensityREF_arr = ProtonDensityREF_meta.get_fdata()
            Ratio_correction = Overlap*(ProtonDensityREF_arr/ProtonDensityTocorrect_arr)
            Ratio_correction = Ratio_correction[~(np.isinf(Ratio_correction)|np.isnan(Ratio_correction)|(Ratio_correction==0))]
            Ratio_correction = np.mean(Ratio_correction)
            print(Ratio_correction)
            ProtonDensity_corrected_arr = Ratio_correction*ProtonDensityTocorrect_arr
            ProtonDensity_corrected_meta = ni.Nifti1Image(ProtonDensity_corrected_arr,
                                                          ProtonDensityTocorrect_meta.affine,
                                                          ProtonDensityTocorrect_meta.header)
            ni.save(ProtonDensity_corrected_meta, ProtonDensity_corrected_filename)
    return None 

def run_proton_density_all_blocks(ArchitecturePath, nbrBlocks):
    
    for i in range(nbrBlocks):
        JsonFilename = os.path.join(ArchitecturePath, 'bloc_'+str(i+1)+'_intensity.json' )
        FOVmaterialFilename = load_jsonfile(JsonFilename)
        correct_proton_density_intensities(FOVmaterialFilename)
    return None

def run_Fusion_QMRI(ArchitecturePath, nbrBlocks):
    BlocksPath = os.path.join(ArchitecturePath, "04-Blocks")
    WholeHemispherePath = os.path.join(ArchitecturePath, "05-Output")
    for value in QMAP_VALUES:
        
        for i in range(nbrBlocks):
            print_message("Merge "+value+ ' FOVs from block n°' + str(i+1))
            
            JsonFilename = os.path.join(ArchitecturePath, 'bloc_'+str(i+1)+'.json' )
            FOVmaterialFilename = load_jsonfile(JsonFilename)
            outputFilename =  os.path.join(BlocksPath, 'bloc_'+str(i+1)+'_'+value+'.nii.gz' )
                        
            if len(FOVmaterialFilename[value])>0 :
                if not(os.path.exists(outputFilename)):
                    fus1 = merge_QMRI_FOVs(ArchitecturePath, JsonFilename, value)
                    ni.save(fus1, outputFilename)
                else : 
                    print('Reconstructed block already exist')
            else :
                print('No files found for '+ value + ' reconstruction.')
        
        print_message("Merge "+ value + ' block ' + str(i+1))
        ReconstructedHemisphereFilename = os.path.join(WholeHemispherePath, "Reconstructed_"+value+".nii.gz")
        if not(os.path.exists(ReconstructedHemisphereFilename)) and len(FOVmaterialFilename[value])>0:
            ReconstructedWH_nifti = merge_QMRI_blocks(ArchitecturePath, value, nbrBlocks)  
            ni.save(ReconstructedWH_nifti, ReconstructedHemisphereFilename)
        else :
            print("Whole Hemisphere reconstruction already exist")
    return None
    
def RunPipeline(ArchitecturePath, nbrBlocks,  runMerger):
    if not(runMerger) :
        GeoMeanPath = os.path.join(ArchitecturePath, "01-GeoMean")
        InitSpacePath = os.path.join(ArchitecturePath, "02-InitSpace")
        RefSpacePath = os.path.join(ArchitecturePath, "03-RefSpace")
        
        ImaGenerator = glob.iglob(os.path.join(GeoMeanPath, '*.ima'))
        print_message("CREATE MASK IN REFSPACE FROM GEOMEAN FILES")
        for GeoMeanfilename in ImaGenerator:
            suffixe = GeoMeanfilename.split('.')[0].split('_')[-1]
            OutputGeoMeanMask_gis = os.path.join(GeoMeanPath, "Mask_"+suffixe+'.ima')
            OutputGeoMeanMask_RefSpace_nifti = os.path.join(RefSpacePath, "Mask_"+suffixe+'.nii.gz')
            create_mask_from_geomean(ArchitecturePath,
                                     GeoMeanfilename,
                                     OutputGeoMeanMask_gis,
                                     OutputGeoMeanMask_RefSpace_nifti)
        
        print_message("CONVERT GIS FILES IN INITSPACE")
        ImaGenerator = glob.iglob(os.path.join(InitSpacePath, '*.ima'))
        for GISFilenameInInitSpaceFolder in ImaGenerator:
            OutputNiftiFilename = GISFilenameInInitSpaceFolder[:-3]+'nii.gz'
            if not(os.path.exists(OutputNiftiFilename)):
                ConversionGisToNifti(GISFilenameInInitSpaceFolder, OutputNiftiFilename)
            
        print_message("SEND NIFTI FILE IN INITSPACE TO REFSPACE")
        ImaGenerator = glob.iglob(os.path.join(InitSpacePath, '*.nii.gz'))
        for NiftiFilenameInInitSpaceFolder in ImaGenerator:
            prefix_values = NiftiFilenameInInitSpaceFolder.split('/')[-1].split('.')[0].split('_')
            # print(prefix_values, NiftiFilenameInInitSpaceFolder)
            JsonFilename = "SendToRefSpace_" + prefix_values[0] + "_" + prefix_values[1] + '.json'
            Output_RefSpace_nifti = os.path.join(RefSpacePath,
                                                 'RefSpace_' + NiftiFilenameInInitSpaceFolder.split('/')[-1])
            if not(os.path.exists(Output_RefSpace_nifti)):
                send_to_RefSpace(NiftiFilenameInInitSpaceFolder,
                             Output_RefSpace_nifti,
                             os.path.join(ArchitecturePath, JsonFilename),
                             "Linear")
        
        print_message("CORRECT PROTON-DENSITY INTENSITIES")
        run_proton_density_all_blocks(ArchitecturePath, int(nbrBlocks))
        
        # print_message("CORRECT DWI INTENSITIES")
        # run_DWI_all_blocks(ArchitecturePath, int(nbrBlocks))
        
    else :
        print_message("MERGE QMRI")
        run_Fusion_QMRI(ArchitecturePath, int(nbrBlocks))
    
    return None


parser = optparse.OptionParser()
parser.add_option( '-p', '--path',
                    dest = 'ArchitecturePath',
                    help = 'Architectural file path' )
parser.add_option( '-n', '--nbrBlocks',
                    dest = 'nbrBlocks',
                    help = 'Number of blocks' )
parser.add_option( '-r', '--run',
                    dest = 'runMerger',
                    action = "store_true",
                    help = 'Run Merger' )


( options, args ) = parser.parse_args()


################################################################################
# 1) Value Extractor 
################################################################################

RunPipeline(options.ArchitecturePath,
            options.nbrBlocks,
            options.runMerger) 















# def send_to_RefSpace(InputFilename, OutputFilename, RefSpaceFilename, Interpolator, tparams):
#     str_command = ["antsApplyTransforms",
#                     "-d", "3",
#                     "-i", InputFilename,
#                     "-r", RefSpaceFilename,
#                     "-o", OutputFilename,
#                     "-n", Interpolator]
#     for i in range(len(tparams)):
#         tparams.insert(2*i, "-t")
#     str_command = str_command + tparams + ["-v", "1"]
#     print(str_command)
#     subprocess.run(str_command)
#     return None