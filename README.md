## Identifying Lesions

Assume we have example data in OFAMS00001/ples_lpa_mrFLAIR.nii which is floating point with some probability value per voxel. The value 0 represents the background. Values larger than 0 represent voxel that are part of a lesion.

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

Using the above ConnectedComponents example we an perform a T1.nii InPainting using the output
summary lesion field (each lesion is encoded using a whole number > 0):
```
./InPainting OFAMS00001/T1.nii /tmp/ples_lpa_mrFLAIR_label.nii /tmp/
```
Or using the docker container from above
```
docker run --rm -it --entrypoint "/ConnectedComponents/InPainting" connectedcomponents OFAMS00001/T1.nii /tmp/ples_lpa_mrFLAIR_label.nii /tmp/
```

For lesions that are close to the border of white matter the interpolation might not be correct as it would blurr intensities from different tissue types across the lesion volume. Instead it might be more appropriate to limit the intensities for interpolation to the voxel of a single material. For these pusposes you can provide an additional mask argument - limiting the sample points for the interpolation to the white matter material only.

## Lesion distance measures

In order to quantify the location of a lesion relative to the cortical surface an approach can be used that calculates curvilinear (geodesic) distances between two labels, the ventricles located in the center of the brain and the cortical gray to white matter surface. This approach mimics the general direction of the path neurons travel during cortex development.

Please visit the https://github.com/mmiv-center/HeatEquation project that implements such a method.

## Lesion matching over time

To match two lesion segmented volumes the "MatchPair" function can be used:
```
./MatchPair abc/20161111/location_labeled_lesions.nii.gz \
            abc/20171104/location_labeled_lesions.nii.gz abc/output/ \
            -i abc
```
The program will calculate the connected components in both volumes and register the point clouds of found lesion centroids with each other. This is a crude way of performing a registration of the two time points and works only if there are sufficient number of lesions found. A table (csv) is created with all shape measures and classes for "new lesion", "deleted lesion", and "mapped lesion" as well as ratio of volume change between matching lesions. 

The '-i' option is used to add a column PatientID to the output CSV so each file can be trivially merged across patients.
