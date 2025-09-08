# A method for analysing tissue motion and deformation during mammalian organogenesis

This repository contains the MATLAB 2023b implementation of a workflow designed to analyze and validate myocardial motion and deformation during heart tissue (HT) morphogenesis. :mouse: :anatomical_heart:

**Paper:** *A method for analysing tissue motion and deformation during mammalian cardiogenesis* (Raiola M et al., 2025b)

The order of the code follows the order of Figure 1 in the paper.

## Quick Overview

<details>
  <summary><strong>Requirements</strong></summary>
  <ul>
    <li><a href="#software">Software</a></li>
    <li><a href="#matlab-toolboxes">MATLAB Toolboxes</a></li>
  </ul>
</details>

<details>
  <summary><strong>Workflow Overview</strong></summary>
  <ul>
    <li><a href="#1-estimating-individual-live-image-motion">1. Estimating Individual Live Image Motion</a>
      <ul>
        <li><a href="#11-image-preprocessing">1.1 Image Preprocessing</a></li>
        <li><a href="#12-mirt-algorithm">1.2 MIRT Algorithm</a></li>
        <li><a href="#13-ht-segmentation">1.3 HT Segmentation</a></li>
        <li><a href="#14-continuous-description-of-ht-morphogenesis">1.4 Continuous Description of HT Morphogenesis</a></li>
        <li><a href="#15-validating-motion-estimation">1.5 Validating Motion Estimation</a></li>
        <li><a href="#16-validating-ht-morphogenesis-description">1.6 Validating HT Morphogenesis Description</a></li>
      </ul>
    </li>
    <li><a href="#2-integrating-multiple-live-images-into-a-consensus-temporal-reference">2. Integrating Multiple Live Images into a Consensus Temporal Reference</a>
      <ul>
        <li><a href="#21-staging-system">2.1 Staging System</a>
          <ul>
            <li><a href="#211-morphometric-feature-definition">2.1.1 Morphometric Feature Definition</a></li>
            <li><a href="#212-staging-system-modeling">2.1.2 Staging System Modeling</a></li>
          </ul>
        </li>
        <li><a href="#22-spatial-mapping">2.2 Spatial Mapping</a>
          <ul>
            <li><a href="#221-rigid-registration-and-masking">2.2.1 Rigid Registration and Masking</a></li>
            <li><a href="#222-validating-spatial-correspondences-between-atlas-and-live-shape">2.2.2 Validating Spatial Correspondences between ATLAS and Live-Shape</a></li>
          </ul>
        </li>
      </ul>
    </li>
  </ul>
</details>

## Requirements

### Software
- [Meshlab](https://www.meshlab.net/)
- [Fiji](https://imagej.net/software/fiji/downloads)
- [ITK-SNAP](http://www.itksnap.org/pmwiki/pmwiki.php)

### MATLAB Toolboxes
- [Mij Toolbox](https://es.mathworks.com/matlabcentral/fileexchange/47545-mij-running-imagej-and-fiji-within-matlab)
- MIRT Toolbox
- [Saveastiff Toolbox](https://es.mathworks.com/matlabcentral/fileexchange/35684-multipage-tiff-stack)
- [Iso2mesh Toolbox](https://iso2mesh.sourceforge.net/cgi-bin/index.cgi)
- [Toolbox_graph](https://github.com/gpeyre/matlab-toolboxes/tree/master/toolbox_graph)
- [Geom3d Toolbox](https://es.mathworks.com/matlabcentral/fileexchange/24484-geom3d)
- [Groupwise Registration Toolbox](https://es.mathworks.com/matlabcentral/fileexchange/63693-robust-group-wise-registration-of-point-sets-using-multi-resolution-t-mixture-model)
- [Meshes3d](https://github.com/mattools/matGeom)
- [Geodesic MATLAB Master](https://es.mathworks.com/matlabcentral/fileexchange/18168-exact-geodesic-for-triangular-meshes)
- [ObjWriter](https://github.com/JBKacerovsky/objWriter?tab=readme-ov-file)
- 
### DATA
The dataset used in this project is publicly available on Mendeley Data @ doi: 10.17632/54gbvnsgnp.1

## Workflow Overview

### 1. Estimating Individual Live Image Motion

#### 1.1 Image Preprocessing
**Script:** [Preprocessing.m](./1.EstimatingIndividualLiveImageMotion/1.1Preprocessing.m)

This step involves preprocessing 3D+t images:
- Loading images via MIJ (MATLAB-ImageJ).
- Applying Gaussian and Median filtering (if required).
- Cropping and resizing images (down to 25%).
- Reslicing volumes starting from the left side using Fiji.

Embryo = '... \Embryo1\Data\Embryo1.tif'; % Raw Data 
Output = '... \Embryo1\Data'; % Resampled Data at 25%

#### 1.2 MIRT Algorithm
**Script:** [MIRT.m](./1.EstimatingIndividualLiveImageMotion/1.2MIRT.m)

In this step, live images are registered using the MIRT3D algorithm:
- Registration begins from the midpoint (N/2) of the time series (N = number of frames).
- Outputs transformation sets to warp images sequentially.

#### 1.3 HT Segmentation
- Segment heart tissues using ITK-SNAP and convert processed .tif images to .mha format using the 3D IO ImageJ plugin.
  - **Myocardium:** Label = 1
  - **Splanchnic Mesoderm:** Label = 2 (linked in the middle)
  - **Endoderm:** Label = 3
- Export segmented images as NIFTI (.nii.gz) and reslice from the left using ImageJ, then save as .tif.

Output : *\\1. Estimating Individual Live Image Motion\\Embryo1\\ITK_Data*

#### 1.4 Continuous Description of HT Morphogenesis
**Script:** [ContinuousHTMorphogenesis.m](./1.EstimatingIndividualLiveImageMotion/1.4ContinousHTMorphogenesis.m)

This script creates a continuous description of HT morphogenesis:
- Generate triangular meshes from the segmented heart tissue using the Iso2mesh toolbox.
- Interpolate the mesh at time point T(N/2), both backward and forward in time, to produce a 3D+t mesh sequence describing continuous morphogenesis.

Output : *\\1. Estimating Individual Live Image Motion\\Embryo1\\Shape*

#### 1.5 Validating Motion Estimation
**Script:** [Error1.m](./1.EstimatingIndividualLiveImageMotion/1.5Error1.m)

This script evaluates the accuracy of motion estimation:
- Compare manual tracking (ground truth) with the propagated tracking (punctual and sequential test sets).
- Calculate the error in micrometers (µm) as the Eulerian distance.

#### 1.6 Validating HT Morphogenesis Description
**Script:** [Error2.m](./1.EstimatingIndividualLiveImageMotion/1.6Error2.m)

This script validates the accuracy of tracking cell division during morphogenesis:
- Segment eight dividing cells during the morphogenesis of embryo e02 using ITK-SNAP.
- Convert the segmented cells into meshes, and reshape the segmentation as in step 4.
- Compare the direction of cell division by fitting lines to the vertices of the 3D mesh and using the cosine similarity function.

### 2. Integrating Multiple Live Images into a Consensus Temporal Reference

#### 2.1 Staging System

##### 2.1.1 Morphometric Feature Definition
**Script:** [FeatureExtraction.m](./2.IntegratingMultipleLiveImagesIntoAConsensusTemporalReference/StagingSystem/2.1.1FeatureExtraction.m)

- Manually select the points pt1-pt2-pt3-pt4 following the rules reported in the manuscript.
- Compute the Eulerian distances between the points, followed by the height/width (h/w) computation.

##### 2.1.2 Staging System Modeling
**Script:** [StagingSystem.m](./2.IntegratingMultipleLiveImagesIntoAConsensusTemporalReference/StagingSystem/2.1.2StagingSystem.m)

- Build a Gaussian Naive Bayesian classifier on h/w Atlas features.
- Predict the staging grade (Gr) for each frame of each live image.
- Among different equal-staged frames, consider only the one with the highest probability.

#### 2.2 Spatial Mapping
**Script:** [SpatialMapping.m](./2.IntegratingMultipleLiveImagesIntoAConsensusTemporalReference/SpatialMapping/2.2SpatialMapping.m)

1. Perform rigid registration of Atlas to live-shape using the TGMM algorithm.
2. Manually cut Atlas missing parts in MeshLab ((*)_Cut.ply) and fill gaps in MeshLab ((*)_Cut.ply).
3. Mask Atlas by creating the ATLAS_Cut mask with the surf2volz function, recreating a binary image with the same size as the live-shape mask.
4. Import the segmentation into Fiji, manually fill gaps, and close holes to transition from surface segmentation to whole segmentation (Fill_(*).tif).
5. Mask the live-shape, creating it with continuous HT description (image).
6. Select only staged frames from the live-shape mask.
7. Perform non-rigid registration with MIRT, transforming the live-shape mask into the ATLAS_Cut mask. The output is the deformed live-shape mask and the transformation T.
8. Morph the live-shape into the Atlas, aligning the live-shape point cloud onto the Atlas mask edge. Output SurfaceMap.
9. Perform face-to-face matching between live-shape and ATLAS_Cut. Return IdxCUT.mat (closest point in Atlas for each point in ATLAS_Cut) and IdxMatch.mat (closest point in Mapped for each point in Atlas (IdxCut)).

##### 2.2.1 Validating Spatial Correspondences between ATLAS and Live-Shape
**Script:** [Validation.m](./2.IntegratingMultipleLiveImagesIntoAConsensusTemporalReference/SpatialMapping/2.2.1Validation.m)

- Validate spatial correspondences between Atlas and live-shape.
- Compute the live-shape area mesh and map face-to-face values into the Atlas.

### 3. Extracting Tissue Deformation

#### 3.1 Individual Tissue Deformation
**Script:** [ExtractingTissueDeformation.m](./3.QuantifyingTissueDeformation/3.1ExtractingTissueDeformation.m)

Extract the mesh deformation between the rest and deformed shapes Gr-Gr+1.

#### 3.2 Mapping Individual Deformation In Atlas
**Script:** [MappingIndividualDeformationInATLAS.m](./3.QuantifyingTissueDeformation/3.2MappingIndividualDeformationInATLAS.m)

Plotting Deformation onto the ATLAS Shape using Face-to-Face Matching.

#### 3.3 HT Cumulative Deformation
**Script:** [HTCumulativeDeformation.m](./3.QuantifyingTissueDeformation/3.3HTCumulativeDeformation.m)

### 4. In Silico Fate Map
**Script:** [ATLASMotionProfile](./4.InSilicoFateMap/ATLASMotionProfile.m)

Compute the dynamic ATLAS by combining the motion profiles of individual SurfaceMaps.
