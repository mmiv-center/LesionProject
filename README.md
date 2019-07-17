### Identifying lessions

Assume we have example data in OFAMS00001/ples_lpa_mrFLAIR.nii which is floating point with some probability values per voxel. The value 0 represents the background.

```
docker run --rm -it connectedcomponents OFAMS00001/ples_lpa_mrFLAIR.nii /tmp/
```
will write a file representing a binary mask) for each detected connected component to the /tmp/ directory. Also the call will save a .json file in the same directory with voxel information (volume, roundness, elongation, principal moments, ...). Create a spreadsheet from this file using
```
jq -r '.lesions | map(.filename), map(.id), map(.num_voxel), map(.flatness), map(.roundness), map(.elongation) | @csv' /tmp//ples_lpa_mrFLAIR_label.json > /tmp//ples_lpa_mrFLAIR_label.csv
```

### Create the dockerized version of ConnectedComponents

```
docker build -t connectedcomponents -f Dockerfile .
docker run --rm -it connectedcomponents
```
