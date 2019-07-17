## Identifying Lesions

Assume we have example data in OFAMS00001/ples_lpa_mrFLAIR.nii which is floating point with some probability values per voxel. The value 0 represents the background. Values larger than 0 represent voxel that are part of a lesion.

The call
```
docker run --rm -it connectedcomponents OFAMS00001/ples_lpa_mrFLAIR.nii /tmp/
```
will write a file (representing a binary mask) for each detected connected component in the input to the /tmp/ directory. Also, the call will save a .json file in the same directory with voxel information (volume, roundness, elongation, principal moments, ...). In order to create a spreadsheet from this JSON file use a call like:
```
jq -r '.lesions | map(.filename), map(.id), map(.num_voxel), map(.flatness), map(.roundness), map(.elongation) | @csv' /tmp//ples_lpa_mrFLAIR_label.json > /tmp//ples_lpa_mrFLAIR_label.csv
```

An additional file "/tmp/ples_lpa_mrFLAIR_label.nii" contains all lesions identified by increasing label values.

### Options

You may change the threshold (-t 0.00001) used to binarize the input image. The program is using second maximum threshold value of 100 that cannot be changed. Voxel between the lower and upper threshold are treated as potential lesions.

You may filter out lesions that are too small by specifying a threshold (-m 1) for the minimum number of voxel in a single lesion that
should appear in the output.

### Create the dockerized version of ConnectedComponents

Checkout this repository and build the containerized version:
```
> git clone https://github.com/mmiv-center/LesionProject.git
> cd LesionProject
> docker build -t connectedcomponents -f Dockerfile .
> docker run --rm -it connectedcomponents
Option infile is required but not defined
Option outdir is required but not defined
 Command tags: 
   [ -t < threshold > ]
      = Specify the threshold applied to the input to create a mask (0.00001).
   [ -m < minPixel > ]
      = Specify the minimum number of voxel in a lesion (1).
   [ -v ]
      = Print more verbose output
 Command fields: 
   < infile > 
      = Input mask
   < outdir > 
      = Output masks directory
```
