// The idea is simple
// 1) do a connected components on the lesion input
// 2) do a 3D region growing based on each lesion for 3+maskBorder iterations
// 3) use the volume ring of the last 3 iterations and sample from all input image
//    voxel that have a mask image > 0, use these intensities and locations as the
//    starting voxel for the in-painting operation
// 4) set the in-painting voxel intensities into a copy of the input image and export that

//#include "itkGDCMImageIO.h"
//#include "itkGDCMSeriesFileNames.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
//#include "itkImageSeriesReader.h"
//#include "itkMetaDataObject.h"
//#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
//#include "itkBinaryFillholeImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
//#include "itkExtractImageFilter.h"
//#include "itkPasteImageFilter.h"
//#include "itkDiscreteGaussianImageFilter.h"
//#include "itkHessianRecursiveGaussianImageFilter.h"
//#include "itkImageAdaptor.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkLabelObject.h"
#include "itkLabelShapeKeepNObjectsImageFilter.h"

#include "itkLabelImageToLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelSelectionLabelMapFilter.h"

//#include "itkMinimumMaximumImageCalculator.h"
//#include "itkRGBPixel.h"
//#include "itkSliceBySliceImageFilter.h"
//#include "itkSymmetricEigenAnalysisImageFilter.h"
//#include "itkSymmetricSecondRankTensor.h"
//#include <itkPixelAccessor.h>
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"

//#include "itkBSplineInterpolateImageFunction.h"
//#include "itkExtractImageFilter.h"
//#include "itkResampleImageFilter.h"
//#include "itkScalarImageToHistogramGenerator.h"
//#include "itkWindowedSincInterpolateImageFunction.h"

//#include "gdcmAnonymizer.h"
//#include "gdcmAttribute.h"
//#include "gdcmDataSetHelper.h"
//#include "gdcmFileDerivation.h"
//#include "gdcmFileExplicitFilter.h"
//#include "gdcmGlobal.h"
//#include "gdcmImageApplyLookupTable.h"
//#include "gdcmImageChangePlanarConfiguration.h"
//#include "gdcmImageChangeTransferSyntax.h"
//#include "gdcmImageHelper.h"
//#include "gdcmImageReader.h"
//#include "gdcmImageWriter.h"
//#include "gdcmMediaStorage.h"
//#include "gdcmRescaler.h"
//#include "gdcmStringFilter.h"
//#include "gdcmUIDGenerator.h"
//#include "itkConstantPadImageFilter.h"
//#include "itkShrinkImageFilter.h"

//#include "itkGDCMImageIO.h"

float pointValue(float x, float y, float z, float power, float smoothing, std::vector<int> xv, std::vector<int> yv, std::vector<int> zv, std::vector<float> values) {
  float nominator = 0.0f;
  float denominator = 0.0f;
  for (int i = 0; i < values.size(); i++) {
    float dist = sqrt( (x-xv[i])*(x-xv[i]) + (y-yv[i])*(y-yv[i]) + (z-zv[i])*(z-zv[i]) + smoothing*smoothing);
    if (dist<0.0000000001)
      return values[i];  
    nominator = nominator+(values[i]/pow(dist,power));
    denominator = denominator+(1.0/pow(dist,power));
  }
  float value = 0.0f;
  if (denominator > 0)  
    value = nominator/denominator;  
  else
    value = -9999;
  return value;
}

/* 
from math import pow  
from math import sqrt  
import numpy as np  
import matplotlib.pyplot as plt  

def pointValue(x,y,power,smoothing,xv,yv,values):  
    nominator=0  
    denominator=0  
    for i in range(0,len(values)):  
        dist = sqrt((x-xv[i])*(x-xv[i])+(y-yv[i])*(y-yv[i])+smoothing*smoothing);  
        #If the point is really close to one of the data points, return the data point value to avoid singularities  
        if(dist<0.0000000001):  
            return values[i]  
        nominator=nominator+(values[i]/pow(dist,power))  
        denominator=denominator+(1/pow(dist,power))  
    #Return NODATA if the denominator is zero  
    if denominator > 0:  
        value = nominator/denominator  
    else:  
        value = -9999  
    return value  

def invDist(xv,yv,values,xsize=100,ysize=100,power=2,smoothing=0):  
    valuesGrid = np.zeros((ysize,xsize))  
    for x in range(0,xsize):  
        for y in range(0,ysize):  
            valuesGrid[y][x] = pointValue(x,y,power,smoothing,xv,yv,values)  
    return valuesGrid  


if __name__ == "__main__":  
    power=1  
    smoothing=20  

    #Creating some data, with each coodinate and the values stored in separated lists  
    xv = [10,60,40,70,10,50,20,70,30,60]  
    yv = [10,20,30,30,40,50,60,70,80,90]  
    values = [1,2,2,3,4,6,7,7,8,10]  

    #Creating the output grid (100x100, in the example)  
    ti = np.linspace(0, 100, 100)  
    XI, YI = np.meshgrid(ti, ti)  

    #Creating the interpolation function and populating the output matrix value  
    ZI = invDist(xv,yv,values,100,100,power,smoothing)  


    # Plotting the result  
    n = plt.normalize(0.0, 100.0)  
    plt.subplot(1, 1, 1)  
    plt.pcolor(XI, YI, ZI)  
    plt.scatter(xv, yv, 100, values)  
    plt.title('Inv dist interpolation - power: ' + str(power) + ' smoothing: ' + str(smoothing))  
    plt.xlim(0, 100)  
    plt.ylim(0, 100)  
    plt.colorbar()  

    plt.show() 
    */



//#include "itkMetaDataDictionary.h"
#include "json.hpp"
#include "metaCommand.h"
#include <boost/filesystem.hpp>
#include <map>

using json = nlohmann::json;
using namespace boost::filesystem;

json resultJSON;

int main(int argc, char *argv[]) {

  itk::MultiThreaderBase::SetGlobalMaximumNumberOfThreads(4);

  MetaCommand command;
  command.SetAuthor("Hauke Bartsch");
  command.SetDescription("InPainting an intensity image to fill in values in small islands defined in a binary mask. The allowed voxel from which intensities are sampled can be specified by another (white matter) mask image.");
  command.AddField("imagefile", "Input intensity volume", MetaCommand::STRING, true);
  command.AddField("lesionfile", "Input lesion volume", MetaCommand::STRING, true);
  command.AddField("outdir", "Output directory for in-painted volume", MetaCommand::STRING, true);


  command.SetOption("maskfile", "m", false, "Input mask volume for white matter");
  command.AddOptionField("maskfile", "maskfile", MetaCommand::STRING, true);

  command.SetOption("borderPixel", "b", false, "Specify a border in pixel around the lesion (2). No voxel from this perimeter will be sampled.");
  command.AddOptionField("borderPixel", "borderPixel", MetaCommand::INT, true);

  command.SetOption("Verbose", "V", false, "Print more verbose output");

  if (!command.Parse(argc, argv)) {
    return 1;
  }

  std::string image = command.GetValueAsString("imagefile");
  std::string lesions = command.GetValueAsString("lesionfile");
  // std::string mask = command.GetValueAsString("maskfile");
  std::string outdir = command.GetValueAsString("outdir");

  if (!boost::filesystem::exists(image)) {
    std::cout << "Could not find the input file..." << std::endl;
    exit(1);
  }
  if (!boost::filesystem::exists(lesions)) {
    std::cout << "Could not find the lesions file..." << std::endl;
    exit(1);
  }

  int borderPixel = 2;
  if (command.GetOptionWasSet("borderPixel"))
    borderPixel = command.GetValueAsInt("borderPixel", "borderPixel");

  std::string mask_filename = "";
  if (command.GetOptionWasSet("mask"))
    mask_filename = command.GetValueAsString("maskfile", "maskfile");

  bool verbose = false;
  if (command.GetOptionWasSet("Verbose"))
    verbose = true;

  // store information in the result json file
  resultJSON["command_line"] = json::array();
  for (int i = 0; i < argc; i++) {
    resultJSON["command_line"].push_back(std::string(argv[i]));
  }
  path p(image);
  std::string fn = p.filename().string();
  size_t lastdot = fn.find_last_of(".");
  std::string output_filename;
  if (lastdot == std::string::npos)
    output_filename = fn + "_inpainted.nii";
  else
    output_filename = fn.substr(0, lastdot) + "_inpainted.nii";

  resultJSON["output_volume"] = outdir + "/" + output_filename;

  constexpr unsigned int ImageDimension = 3;
  using PixelType = float;
  typedef itk::Image<float, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ImageReaderType;
  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName(image);
  imageReader->Update();

  // after importing the intensity image also import the lesion mask
  using MaskPixelType = unsigned short;
  using MaskImageType = itk::Image<MaskPixelType, ImageDimension>;
  typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
  MaskReaderType::Pointer maskReader = MaskReaderType::New();
  maskReader->SetFileName(lesions);
  maskReader->Update();

  MaskReaderType::Pointer whiteMatterReader = MaskReaderType::New();
  if (command.GetOptionWasSet("mask")) {
    if (verbose) {
      fprintf(stdout, "read the white matter mask...\n");
    }
      whiteMatterReader->SetFileName(mask_filename);
      whiteMatterReader->Update();
  }

  typedef itk::ConnectedComponentImageFilter<MaskImageType, MaskImageType> ConnectedComponentImageFilterType;
  ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
  connected->SetBackgroundValue(0);
  connected->SetInput(maskReader->GetOutput());
  connected->Update();

  MaskImageType::Pointer con = connected->GetOutput();
  con->SetOrigin(imageReader->GetOutput()->GetOrigin());
  con->SetSpacing(imageReader->GetOutput()->GetSpacing());
  con->SetDirection(imageReader->GetOutput()->GetDirection());

  // using LabelType = unsigned short;
  using ShapeLabelObjectType = itk::ShapeLabelObject<MaskPixelType, ImageDimension>;
  using LabelMapType = itk::LabelMap<ShapeLabelObjectType>;
  using LabelType = itk::LabelImageToShapeLabelMapFilter<MaskImageType, LabelMapType>;
  LabelType::Pointer label = LabelType::New();
  label->SetInput(connected->GetOutput());
  label->SetComputePerimeter(true);
  label->Update();

  LabelMapType *labelMap = label->GetOutput();
  if (labelMap->GetNumberOfLabelObjects() == 0) {
    // error case
    fprintf(stderr, "Error: Could not find any lesions in the lesion input.\n");
  }

  // do the inpainting in this volume
  ImageType::Pointer outimg = imageReader->GetOutput();

  resultJSON["voxel_size"] = json::array();
  resultJSON["voxel_size"].push_back(outimg->GetSpacing()[0]);
  resultJSON["voxel_size"].push_back(outimg->GetSpacing()[1]);
  resultJSON["voxel_size"].push_back(outimg->GetSpacing()[2]);

  resultJSON["lesions"] = json::array();
  int counter = 0;
  size_t totalVolume = 0;
  for (unsigned int n = 0; n < labelMap->GetNumberOfLabelObjects(); ++n) {
    if (verbose) {
      fprintf(stdout, "process %d of %d lesions...\n", n, labelMap->GetNumberOfLabelObjects());
    }
    ShapeLabelObjectType *labelObject = labelMap->GetNthLabelObject(n);
    //if (labelObject->GetNumberOfPixels() < minPixel)
    //  continue; // ignore this region
    // labelObject->GetNumberOfPixels()
    json lesion;
    lesion["id"] = counter;
    lesion["input_value"] = labelObject->GetLabel();
    lesion["num_voxel"] = labelObject->GetNumberOfPixels();
    lesion["physical_size"] = labelObject->GetPhysicalSize();
    lesion["flatness"] = labelObject->GetFlatness();
    lesion["roundness"] = labelObject->GetRoundness();
    lesion["perimeter"] = labelObject->GetPerimeter();
    lesion["elongation"] = labelObject->GetElongation();
    lesion["number_pixel_on_border"] = labelObject->GetNumberOfPixelsOnBorder();
    lesion["feret_diameter"] = labelObject->GetFeretDiameter();
    lesion["perimeter_on_border_ratio"] = labelObject->GetPerimeterOnBorderRatio();
    lesion["centroid"] = json::array();
    lesion["centroid"].push_back(labelObject->GetCentroid()[0]);
    lesion["centroid"].push_back(labelObject->GetCentroid()[1]);
    lesion["centroid"].push_back(labelObject->GetCentroid()[2]);
    lesion["principal_moments"] = json::array();
    lesion["principal_moments"].push_back(labelObject->GetPrincipalMoments()[0]);
    lesion["principal_moments"].push_back(labelObject->GetPrincipalMoments()[1]);
    lesion["principal_moments"].push_back(labelObject->GetPrincipalMoments()[2]);
    lesion["equivalent_spherical_radius"] = labelObject->GetEquivalentSphericalRadius();
    totalVolume += labelObject->GetNumberOfPixels();

    ImageType::RegionType region = connected->GetOutput()->GetLargestPossibleRegion();
    // now create a new volume just for this lesion

    MaskImageType::Pointer mask = MaskImageType::New();
    mask->SetRegions(region);
    mask->Allocate();
    mask->FillBuffer(0); // density for air
    mask->SetOrigin(imageReader->GetOutput()->GetOrigin());
    mask->SetSpacing(imageReader->GetOutput()->GetSpacing());
    mask->SetDirection(imageReader->GetOutput()->GetDirection());
    itk::ImageRegionIterator<MaskImageType> imageIterator(connected->GetOutput(), region);
    itk::ImageRegionIterator<MaskImageType> maskIterator(mask, region);
    while (!imageIterator.IsAtEnd() && !maskIterator.IsAtEnd()) {
      if (imageIterator.Get() == labelObject->GetLabel()) {
        maskIterator.Set(1);
      }
      ++imageIterator;
      ++maskIterator;
    }
    // ok, we have a single lesion now
    // start with step 1 of making the lesion larger by using a structuring element
    using StructuringElementType = itk::BinaryBallStructuringElement<MaskPixelType, ImageDimension>;
    using ErodeFilterType = itk::BinaryErodeImageFilter<MaskImageType, MaskImageType, StructuringElementType>;
    using DilateFilterType = itk::BinaryDilateImageFilter<MaskImageType, MaskImageType, StructuringElementType>;
    StructuringElementType structuringElement;
    structuringElement.SetRadius(1); // 3x3 structuring element
    structuringElement.CreateStructuringElement();
    // do dilate for the radius borderPixel to make the control points
    // independed of the lesion border
    MaskImageType::Pointer m = mask;
    for (int b = 0; b < borderPixel; b++) {
      DilateFilterType::Pointer binaryDilate = DilateFilterType::New(); // grows inside the tissue

      binaryDilate->SetKernel(structuringElement);
      binaryDilate->SetInput(m);
      binaryDilate->SetDilateValue(1);
      binaryDilate->Update();

      m = binaryDilate->GetOutput();
      m->DisconnectPipeline();
    }
    // ok, now we have a larger mask with the border in m, we want to keep that

    // now we can grow again by 2times maybe to get the voxel that are the seed points
    // for the in-painting
    MaskImageType::Pointer m2 = m;
    for (int b = 0; b < 2; b++) {
      DilateFilterType::Pointer binaryDilate = DilateFilterType::New(); // grows inside the tissue

      binaryDilate->SetKernel(structuringElement);
      binaryDilate->SetInput(m2);
      binaryDilate->SetDilateValue(1);
      binaryDilate->Update();

      m2 = binaryDilate->GetOutput();
      m2->DisconnectPipeline();
    }

    // now compute the difference between the two m2 - m
    /* MaskImageType::Pointer diffMask = MaskImageType::New();
    diffMask->SetRegions(region);
    diffMask->Allocate();
    diffMask->FillBuffer(0); // density for air
    diffMask->SetOrigin(imageReader->GetOutput()->GetOrigin());
    diffMask->SetSpacing(imageReader->GetOutput()->GetSpacing());
    diffMask->SetDirection(imageReader->GetOutput()->GetDirection()); */

    // If we have a mask here we should use it as well to make sure
    // that we don't use background voxel that don't belong to the tissue
    // type.
    itk::ImageRegionIterator<MaskImageType> whiteMatterMaskIterator;
    if (command.GetOptionWasSet("mask")) {
      whiteMatterMaskIterator = itk::ImageRegionIterator<MaskImageType> (whiteMatterReader->GetOutput(), region);
    }
    itk::ImageRegionIterator<MaskImageType> maskIterator1(m, region);
    itk::ImageRegionIterator<MaskImageType> maskIterator2(m2, region);
    itk::ImageRegionIterator<ImageType> imageIterator3(imageReader->GetOutput(), region);
    // ok, store the cooridnates of the pixel and their value
    std::vector<int> xv;
    std::vector<int> yv;
    std::vector<int> zv;
    std::vector<float> iv;
    using Index3DType = MaskImageType::IndexType;
    Index3DType index;
    while (!maskIterator1.IsAtEnd() && !maskIterator2.IsAtEnd() && !imageIterator3.IsAtEnd()) {
      bool bail = false;
      if (command.GetOptionWasSet("mask")) {
        if (whiteMatterMaskIterator.Get() == 0)
            // don't use this voxel 
            bail = true;
      }
      if (maskIterator2.Get() == 1 && maskIterator1.Get() == 0 && !bail) {
        // where are we?
        index = maskIterator1.GetIndex();
        xv.push_back(index[0]);
        yv.push_back(index[1]);
        zv.push_back(index[2]);
        iv.push_back(imageIterator3.Get());
      }
      ++maskIterator1;
      ++maskIterator2;
      ++imageIterator3;
      if (command.GetOptionWasSet("mask")) {
        ++whiteMatterMaskIterator;
      }
    }

    // now we should interpolate and write the values back to the input image
    // before saving that again as the only output
    itk::ImageRegionIterator<MaskImageType> maskIterator11(m, region);
    itk::ImageRegionIterator<ImageType> imageIterator33(outimg, region);
    while (!maskIterator11.IsAtEnd() && !imageIterator33.IsAtEnd()) {
      if (maskIterator11.Get() == 1) { // the inner mask
        index = maskIterator11.GetIndex();
        float val = pointValue(index[0], index[1], index[2], 5.0, 0.01, xv, yv, zv, iv);
        imageIterator33.Set(val);
      }
      ++maskIterator11;
      ++imageIterator33;
    }
    resultJSON["lesions"].push_back(lesion);
    counter++;
  }
  resultJSON["num_lesions"] = counter;
  resultJSON["total_lesion_size"] = totalVolume;

  if (1) { // save output image as nifti again
      typedef itk::ImageFileWriter<ImageType> WriterType;
      WriterType::Pointer writer = WriterType::New();
      // check if that directory exists, create before writing
      std::string fn = resultJSON["output_volume"];
      size_t lastdot = fn.find_last_of(".");
      std::string filename("");
      if (lastdot == std::string::npos)
        filename = fn + ".nii.gz";
      else
        filename = fn.substr(0, lastdot) + ".nii.gz";

      writer->SetFileName(filename);
      writer->SetInput(outimg);

      std::cout << "Writing output " << std::endl;
      std::cout << " to " << filename << std::endl;
      try {
        writer->Update();
      } catch (itk::ExceptionObject &ex) {
        std::cout << ex << std::endl;
        return EXIT_FAILURE;
      }
  }

  std::ostringstream o;
  std::string si(resultJSON["output_volume"]);
  si.erase(std::remove(si.begin(), si.end(), '\"'), si.end());
  lastdot = si.find_last_of(".");
  if (lastdot == std::string::npos)
    si = si + ".json";
  else
    si = si.substr(0, lastdot) + ".json";

  o << si;
  /* resultJSON["z_comment"] =
      std::string("jq -r '.lesions | map(.filename), map(.id), map(.num_voxel), map(.flatness), map(.roundness), map(.elongation) | @csv' ") + o.str(); */

  std::ofstream out(o.str());
  std::string res = resultJSON.dump(4) + "\n";
  out << res;
  out.close();

  fprintf(stdout, "%s", res.c_str());

  return EXIT_SUCCESS;
}
