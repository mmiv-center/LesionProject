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

You may change the threshold (-t 0.00001) used by the program to initially binarize the input image. The program is using a second maximum threshold value of 100 that cannot be changed. Only voxel values between the lower and upper threshold are treated as potential lesions.

You may filter out lesions that are too small by specifying a threshold (-m 1) for the minimum number of voxel in a single lesion that
should appear in the output.

### Create the dockerized version of ConnectedComponents

Checkout this repository and build the containerized version:
```
> git clone https://github.com/mmiv-center/LesionProject.git
> cd LesionProject
> docker build -t connectedcomponents -f Dockerfile .
> docker run --rm -it connectedcomponents -h
Usage : /ConnectedComponents/ConnectedComponents
 System tags: 
   [ -v ] or [ -h ]
      = List options in short format
   [ -V ] or [ -H ]
      = List options in long format
   [ -vxml ] or [ -hxml ] or [ -exportXML ]
      = List options in xml format for BatchMake
   [ --xml ]
      = List options in xml format for Slicer
   [ -vgad ] or [ -hgad ] or [ -exportGAD ]
      = List options in Grid Application Description format
   [ -version ]
      = return the version number
   [ -date ]
      = return the cvs checkout date
 Command tags: 
   [ -t < threshold > ]
      = Specify the threshold applied to the input to create a mask (0).
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

## InPainting

After detecting lesions in-painting can be used to synthetically create a version of the raw data were lesions
are masked with assumed intensity values similar to the neighboring voxel in the data. Such synthetic data can be
used for algorithms that are sensitive to the lesions otherwise - such as FreeSurfer.

The provided algorithm performs a region growing of initially 2 voxel to create a lesion border. This border might
be affected by partial volume effect. Afterwards another 2 voxel morphological grow operation defines a region of
background voxel used for the interpolation of the lesion and the lesion border voxel intensities.