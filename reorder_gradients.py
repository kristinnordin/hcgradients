"""
Reorder Gradient Maps Script

Description:
This script is designed to rearrange subject-level gradients according to the group-level gradients. 
The script rearranges the gradient based on the strength of the correlation between the subject-level and group-level gradients. 
The strongest correlation is used to determine the first gradients followed by the second and third strongest correlations. 
Additionally, based on the sign of the correlation, the script flips the gradients if necessary.
The outputs are reordered maps alongside a CSV file that summarizing the changes.

Execution:
Run the script via the command line with the following arguments:
    --data_dir        Path to the data directory containing the required data.
    --subjects_list   Path to a text file listing the subject IDs to process (one per line).

Expected Inputs:
1. Data Directory (--data_dir):
   - Contains subject-specific folders named by their IDs.
   - Includes a 'masks' folder with:
       * ROI mask files (e.g., 'rh.hippocampus.nii.gz', 'lh.hippocampus.nii.gz').
       * GM mask file (e.g., 'cortex.nii.gz').
   - Includes a 'congrads_group_all' folder with:
       * Group-level connectopic maps for each ROI (e.g., 'rh.hippocampus.cmaps.nii.gz').

2. Subjects List (--subjects_list):
   - A text file with each line containing a subject ID to process.

Outputs:
1. Rearranged connectopic maps.
   - Saved within the respective subject folder in the data directory.
2. A CSV file summarizing the results.
   - Saved in the data directory as 'results.csv'.

Usage Example:
python script_name.py --data_dir /path/to/data_directory --subjects_list /path/to/subjects_list.txt

Created on: Dec 6 2021
last modified: 10 Jan 2025 
@farshad.falahati
"""


import os
import numpy as np
import nibabel as nib
import scipy.stats
import pandas as pd
import argparse


def swap_vol(vol, idx1, idx2):
    vol[:, [idx1, idx2]] = vol[:, [idx2, idx1]]
    return vol
    
    
def main(data_dir, subjects_list_file):
    masks_dir = os.path.join(data_dir, 'masks')
    group_dir = os.path.join(data_dir, 'congrads_group_all')


    # Read subjects list
    with open(subjects_list_file, 'r') as tf:
        subjects = [line.strip() for line in tf if line.strip()]

    rois = ['rh.hippocampus', 'lh.hippocampus', 'rh.putamen', 'lh.putamen', 'rh.cud_acc', 'lh.cud_acc', 'lh.striatum', 'rh.striatum']

    # non-normalized or normalized cmaps
    norm_flag=True

    # Regions of interest 
    lst_dict = []
    df = pd.DataFrame(columns=['Subject', 'ROI', 
                            'CMAP1xGL1', 'CMAP1xGL2', 'CMAP1xGL3', 'G1_corr_rank', 'G1_rank_assigned',
                            'CMAP2xGL1', 'CMAP2xGL2', 'CMAP2xGL3', 'G2_corr_rank', 'G2_rank_assigned', 
                            'CMAP3xGL1', 'CMAP3xGL2', 'CMAP3xGL3', 'G3_corr_rank', 'G3_rank_assigned',
                            'Comments'])

    for roi in rois: 
        print('ROI: %s' % roi)
        if norm_flag:
            roi_dir = roi + '.norm'
        else:
            roi_dir = roi
            
        # load rio masks 
        # one 3D image containes the binary mask of ROI per hemisphere
        mask_file = os.path.join(masks_dir, roi+'.nii.gz')
        print('Load ROI mask from: %s' % mask_file)
        mask_img = nib.load(mask_file)
        mask_data = mask_img.get_data()
        mask_dim = mask_data.shape
        mask_vect = np.reshape(mask_data,(np.prod(mask_dim[0:3])))
        maskIndices = np.where(mask_vect==1)[0]
        
        # load GM mask 
        # one 3D image containes the binary mask of ROI per hemisphere
        gm_mask_file = os.path.join(masks_dir, 'cortex.nii.gz')
        print('Load GM mask from: %s' % gm_mask_file)
        gm_mask_img = nib.load(gm_mask_file)
        gm_mask_data = gm_mask_img.get_data()
        
        # Group-level maps
        gl_cmap_file = os.path.join(group_dir, roi_dir, roi+'.cmaps.nii.gz')
        print('Load group-level cmap from: %s' % gl_cmap_file)
        gl_cmap_img = nib.load(gl_cmap_file)
        gl_cmap_data = gl_cmap_img.get_data()
        gl_cmap_dim = gl_cmap_data.shape
        if len(gl_cmap_dim) <= 3:
            gl_cmap_dim = gl_cmap_dim + (1,)
        gl_cmap_vect = np.reshape(gl_cmap_data,(np.prod(gl_cmap_dim[0:3]),gl_cmap_dim[3]))
        gl_cmap_masked = gl_cmap_vect[maskIndices,:]

        for subject in subjects:
            print('ROI: %s - Subject # %s' % (roi, subject))
            
            # load cmap data
            # one 4D image per subject, including 3 gradient (3 3D images)
            cmap_file = os.path.join(data_dir, subject, roi_dir, roi+'.cmaps.nii.gz')
            print('Load subject cmap file from: %s' % cmap_file)
            cmap_img = nib.load(cmap_file)
            cmap_data = cmap_img.get_data()
            cmap_dim = cmap_data.shape
            if len(cmap_dim) <= 3:
                cmap_dim = cmap_dim + (1,)
            cmap_vect = np.reshape(cmap_data,(np.prod(cmap_dim[0:3]),cmap_dim[3]))
            cmap_masked = cmap_vect[maskIndices,:]

            # load pmap data
            pmap_file = os.path.join(data_dir, subject, roi_dir, roi+'.pmaps.nii.gz')
            print('Load subject pmaps file from: %s' % pmap_file)
            pmap_img = nib.load(pmap_file)
            pmap_data = pmap_img.get_data()

            #for grd in range(0,cmap_masked.shape[1]): 
            #    cmap_masked[:,grd] =  np.multiply(cmap_masked[:,grd], np.sign(np.corrcoef(cmap_masked[:,grd], gl_cmap_masked[:,grd])[0,1])) 
            
            r1o, p1 = scipy.stats.pearsonr(cmap_masked[:,0], gl_cmap_masked[:,0]) 
            r2o, p2 = scipy.stats.pearsonr(cmap_masked[:,0], gl_cmap_masked[:,1]) 
            r3o, p3 = scipy.stats.pearsonr(cmap_masked[:,0], gl_cmap_masked[:,2]) 
            g1o_r = [abs(r1o), abs(r2o), abs(r3o)]
            g1o = g1o_r.index(max(g1o_r))+1

            r4o, p4 = scipy.stats.pearsonr(cmap_masked[:,1], gl_cmap_masked[:,0]) 
            r5o, p5 = scipy.stats.pearsonr(cmap_masked[:,1], gl_cmap_masked[:,1]) 
            r6o, p6 = scipy.stats.pearsonr(cmap_masked[:,1], gl_cmap_masked[:,2]) 
            g2o_r = [abs(r4o), abs(r5o), abs(r6o)]
            g2o = g2o_r.index(max(g2o_r))+1

            r7o, p7 = scipy.stats.pearsonr(cmap_masked[:,2], gl_cmap_masked[:,0]) 
            r8o, p8 = scipy.stats.pearsonr(cmap_masked[:,2], gl_cmap_masked[:,1]) 
            r9o, p9 = scipy.stats.pearsonr(cmap_masked[:,2], gl_cmap_masked[:,2]) 
            g3o_r = [abs(r7o), abs(r8o), abs(r9o)]
            g3o = g3o_r.index(max(g3o_r))+1

            r11, r12, r13 = g1o_r
            r21, r22, r23 = g2o_r
            r31, r32, r33 = g3o_r
            #g1, g2, g3 = g1o, g2o ,g3o 

            if ((g1o == g2o) | (g1o == g3o) | (g2o == g3o)):
                comment = 'FLAGGED'
            elif ((g1o == 1) & (g2o == 2) & (g3o == 3)):
                comment ='IN-ORDER'        
            elif ((g1o != g2o) & (g1o != g3o) & (g2o != g3o)):
                comment = 'RE-ARRANGE'
            
            if (r11 >= r21) & (r11 >= r31):
                g1 = 1
                re1 = 0
                if (r22 >= r32):
                    g2, g3 = (2, 3)
                    re2, re3 = (1, 2)
                else: 
                    g2, g3 = (3, 2)
                    re2, re3 = (2, 1)
            if (r21 >= r11) & (r21 >= r31):
                g2 = 1 
                re2 = 0
                if (r12 > r32): 
                    g1, g3 = (2, 3)
                    re1, re3 = (1, 2)
                else: 
                    g1, g3 = (3, 2)
                    re1, re3 = (2, 1)
            if (r31 >= r11) & (r31 >= r21):
                g3 = 1 
                re3 = 0
                if (r12 > r22): 
                    g1, g2 = (2, 3)
                    re1, re2 = (1, 2)
                else: 
                    g1, g2 = (3, 2)
                    re1, re2 = (2, 1)
                
            re_cmap_data = zeros(shape(cmap_data))
            re_cmap_data[:,:,:,re1] = cmap_data[:,:,:,0]
            re_cmap_data[:,:,:,re2] = cmap_data[:,:,:,1]
            re_cmap_data[:,:,:,re3] = cmap_data[:,:,:,2]
            
            # rearrange corresponding pmaps 
            re_pmap_data = zeros(shape(pmap_data))
            re_pmap_data[:,:,:,re1] = pmap_data[:,:,:,0]
            re_pmap_data[:,:,:,re2] = pmap_data[:,:,:,1]
            re_pmap_data[:,:,:,re3] = pmap_data[:,:,:,2]
            
            re_cmap_vect = np.reshape(re_cmap_data,(np.prod(cmap_dim[0:3]),cmap_dim[3]))
            re_cmap_masked = re_cmap_vect[maskIndices,:]
            for grad in range(0,re_cmap_masked.shape[1]):
                if np.sign(np.corrcoef(re_cmap_masked[:,grad], gl_cmap_masked[:,grad])[0,1]) < 0:
                    if norm_flag:
                        re_cmap_data[:,:,:,grad] = np.multiply(1-re_cmap_data[:,:,:,grad], mask_data)
                        re_pmap_data[:,:,:,grad] = np.multiply(1-re_pmap_data[:,:,:,grad], gm_mask_data)
                    else: 
                        re_cmap_data[:,:,:,grad] = -re_cmap_data[:,:,:,grad]
                        re_pmap_data[:,:,:,grad] = -re_pmap_data[:,:,:,grad]
            
            re_cmap_img = nib.Nifti1Image(re_cmap_data, cmap_img.get_affine(), cmap_img.get_header())
            cmap_outfile = os.path.join(data_dir, subject, roi_dir, roi+'.cmaps.rearranged.nii.gz')
            nib.save(re_cmap_img, cmap_outfile)

            re_pmap_img = nib.Nifti1Image(re_pmap_data, pmap_img.get_affine(), pmap_img.get_header())
            pmap_outfile = os.path.join(data_dir, subject, roi_dir, roi+'.pmaps.rearranged.nii.gz')
            nib.save(re_pmap_img, pmap_outfile)

            lst_dict.append({'Subject':subject, 'ROI':roi, 
                            'CMAP1xGL1':r1o, 'CMAP1xGL2':r2o, 'CMAP1xGL3':r3o, 'G1_corr_rank':g1o, 'G1_rank_assigned':g1, 
                            'CMAP2xGL1':r4o, 'CMAP2xGL2':r5o, 'CMAP2xGL3':r6o, 'G2_corr_rank':g2o, 'G2_rank_assigned':g2,
                            'CMAP3xGL1':r7o, 'CMAP3xGL2':r8o, 'CMAP3xGL3':r9o, 'G3_corr_rank':g3o, 'G3_rank_assigned':g3,
                            'Comments':comment})

    df1 = df.append(lst_dict, ignore_index=True)
    output_file = os.path.join(data_dir, 'results.csv')
    df1.to_csv(output_file, index=False)
    print(f'Results saved to {output_file}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Reorder Gradient Maps")
    parser.add_argument("--data_dir", type=str, required=True, help="Path to the data directory.")
    parser.add_argument("--subjects_list", type=str, required=True, help="Path to the subjects list file.")
    args = parser.parse_args()

    main(args.data_dir, args.subjects_list)

