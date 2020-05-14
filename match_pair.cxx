// just read in a volume and run connected components on it

//#include "itkGDCMImageIO.h"
//#include "itkGDCMSeriesFileNames.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
//#include "itkImageSeriesReader.h"
//#include "itkMetaDataObject.h"
//#include "itkSmoothingRecursiveGaussianImageFilter.h"

//#include "itkBinaryBallStructuringElement.h"
//#include "itkBinaryDilateImageFilter.h"
//#include "itkBinaryErodeImageFilter.h"
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
// point cloud registration
#include "itkEuclideanDistancePointMetric.h"
#include "itkLevenbergMarquardtOptimizer.h"
#include "itkPointSetToPointSetRegistrationMethod.h"
#include "itkTranslationTransform.h"

//#include "itkMetaDataDictionary.h"
#include "json.hpp"
#include "metaCommand.h"
#include <algorithm>
#include <boost/filesystem.hpp>
#include <map>
#include <math.h>
#include <set>

// for lesions from icometrix we have
// periventricular (2)
// juxtacortical (3)
// Infratentorial
// Deep white matter lesions

using json = nlohmann::json;
using namespace boost::filesystem;

// forward declaration
void CopyDictionary(itk::MetaDataDictionary &fromDict, itk::MetaDataDictionary &toDict);

template <typename TFilter> class CommandIterationUpdate : public itk::Command {
public:
  typedef CommandIterationUpdate Self;
  typedef itk::Command Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  itkNewMacro(Self);

protected:
  CommandIterationUpdate() {}

public:
  virtual void Execute(itk::Object *caller, const itk::EventObject &event) ITK_OVERRIDE { Execute((const itk::Object *)caller, event); }

  virtual void Execute(const itk::Object *object, const itk::EventObject &event) ITK_OVERRIDE {
    const TFilter *filter = dynamic_cast<const TFilter *>(object);

    if (typeid(event) != typeid(itk::IterationEvent)) {
      return;
    }
    if (filter->GetElapsedIterations() == 1) {
      std::cout << "Current level = " << filter->GetCurrentLevel() + 1 << std::endl;
    }
    std::cout << "  Iteration " << filter->GetElapsedIterations() << " (of " << filter->GetMaximumNumberOfIterations()[filter->GetCurrentLevel()] << ").  ";
    std::cout << " Current convergence value = " << filter->GetCurrentConvergenceMeasurement() << " (threshold = " << filter->GetConvergenceThreshold() << ")"
              << std::endl;
  }
};

template <typename TValue> TValue Convert(std::string optionString) {
  TValue value;
  std::istringstream iss(optionString);

  iss >> value;
  return value;
}

void replaceAll(std::string &s, const std::string &search, const std::string &replace) {
  for (size_t pos = 0;; pos += replace.length()) {
    // Locate the substring to replace
    pos = s.find(search, pos);
    if (pos == std::string::npos)
      break;
    // Replace by erasing and inserting
    s.erase(pos, search.length());
    s.insert(pos, replace);
  }
}

json robustAnalysis(std::map<int, std::vector<double>> change) {
  json res;
  std::map<int, std::vector<double>>::iterator it_map = change.begin();
  while (it_map != change.end()) {
    int key = it_map->first;
    double mean = 0.0f;
    for (int i = 0; i < it_map->second.size(); i++) {
      mean += it_map->second[i];
    }
    if (it_map->second.size() > 0)
      mean /= it_map->second.size();
    double std = 0.0f;
    for (int i = 0; i < it_map->second.size(); i++) {
      std += (it_map->second[i] - mean) * (it_map->second[i] - mean);
    }
    if (it_map->second.size() > 1)
      std = sqrt(std / (it_map->second.size() - 1));

    // we should do a robust mean here, sometimes we end up with very large values...
    // so we don't use entries that are larger than 3 std away from mean
    double robust_mean = 0.0f;
    int robust_n = 0;
    double robust_std = 0.0f;
    for (int i = 0; i < it_map->second.size(); i++) {
      if (fabs(it_map->second[i] - mean) <= 1.959964 * std) { // value from https://en.wikipedia.org/wiki/Standard_deviation for 1/20 95%
        robust_mean += it_map->second[i];
        robust_n++;
      }
    }
    if (robust_n > 0)
      robust_mean /= robust_n;
    for (int i = 0; i < it_map->second.size(); i++) {
      if (fabs(it_map->second[i] - mean) <= 3 * std) {
        robust_std += (it_map->second[i] - robust_mean) * (it_map->second[i] - robust_mean);
      }
    }
    if (robust_n > 1) {
      robust_std = sqrt(robust_std / (robust_n - 1));
    }

    json entry;
    entry["mean"] = mean;
    entry["n"] = it_map->second.size();
    entry["std"] = std;
    entry["robust_mean"] = robust_mean;
    entry["robust_n"] = robust_n;
    entry["robust_std"] = robust_std;
    res[std::to_string(key)] = entry;
    ++it_map;
  }
  return res;
}

std::map<int, std::vector<double>> groupBy(std::vector<std::map<std::string, std::string>> csv, std::string group, std::string measure) {
  std::map<int, std::vector<double>> repr2relativeSizeChange;
  for (int i = 0; i < csv.size(); i++) {
    std::string val = "";
    std::string val2 = "";
    std::map<std::string, std::string>::iterator it = csv[i].find(measure);
    if (it != csv[i].end()) {
      val = it->second;
    }
    std::map<std::string, std::string>::iterator it2 = csv[i].find(group);
    if (it2 != csv[i].end()) {
      val2 = it2->second;
      if (val2 == std::string(""))
        continue;
      int k = std::stoi(val2);
      // fprintf(stdout, "key is: %d for \"%s\"\n", k, val2.c_str());
      std::map<int, std::vector<double>>::iterator it3 = repr2relativeSizeChange.find(k);
      if (it3 == repr2relativeSizeChange.end()) {
        // fprintf(stdout, "Add a map key for : %d\n", k);
        repr2relativeSizeChange.insert({k, std::vector<double>()});
      }
      if (val != std::string("")) {
        // fprintf(stdout, "val as float is: %f\n", std::stof(val));
        // add the value to that entry
        std::map<int, std::vector<double>>::iterator it4 = repr2relativeSizeChange.find(k);
        if (it4 != repr2relativeSizeChange.end()) {
          std::vector<double> old = it4->second;
          old.push_back(std::stof(val));
          it4->second = old;
        }
      }
    }
  }
  return repr2relativeSizeChange;
}

json resultJSON;

int main(int argc, char *argv[]) {

  itk::MultiThreaderBase::SetGlobalMaximumNumberOfThreads(4);

  MetaCommand command;
  command.SetAuthor("Hauke Bartsch");
  command.SetDescription("Match pairs of lesions. This program requires already segmented lesions in two input volumes.");
  command.AddField("fixed", "Input lesion mask for fixed volume (t0)", MetaCommand::STRING, true);
  command.AddField("moving", "Input lesion mask moving volume (t1)", MetaCommand::STRING, true);
  command.AddField("outdir", "Output directory", MetaCommand::STRING, true);

  command.SetOption("Threshold", "t", false, "Specify the threshold applied to the input to create a mask (0.00001).");
  command.AddOptionField("Threshold", "threshold", MetaCommand::FLOAT, true);

  command.SetOption("minPixel", "m", false, "Specify the minimum number of voxel in a lesion (1).");
  command.AddOptionField("minPixel", "minPixel", MetaCommand::INT, true);

  command.SetOption("PatientID", "i", false, "Provide a patient id used in the result spreadsheet as ID.");
  command.AddOptionField("PatientID", "patientid", MetaCommand::STRING, true);

  command.SetOption("Verbose", "v", false, "Print more verbose output");

  command.SetOption("exportIndividualLables", "e", false, "Export individual label files for each lesion.");

  if (!command.Parse(argc, argv)) {
    return 1;
  }

  std::string fixed = command.GetValueAsString("fixed");
  std::string moving = command.GetValueAsString("moving");
  std::string outdir = command.GetValueAsString("outdir");

  if (!boost::filesystem::exists(fixed)) {
    std::cout << "Could not find the fixed file..." << std::endl;
    exit(1);
  }
  if (!boost::filesystem::exists(moving)) {
    std::cout << "Could not find the moving file..." << std::endl;
    exit(1);
  }

  std::string PatientID = "";
  if (command.GetOptionWasSet("PatientID")) {
    PatientID = command.GetValueAsString("PatientID", "patientid");
    // fprintf(stdout, "got a patient id %s\n", PatientID.c_str());
  }

  int minPixel = 1;
  if (command.GetOptionWasSet("minPixel"))
    minPixel = command.GetValueAsInt("minPixel", "minPixel");

  float threshold = 0.00001; // > 0
  if (command.GetOptionWasSet("Threshold"))
    threshold = command.GetValueAsFloat("Threshold", "threshold");

  bool verbose = false;
  if (command.GetOptionWasSet("Verbose"))
    verbose = true;

  bool exportIndividualLabels = false;
  if (command.GetOptionWasSet("exportIndividualLabels"))
    exportIndividualLabels = true;

  // store information in the result json file
  resultJSON["command_line"] = json::array();
  for (int i = 0; i < argc; i++) {
    resultJSON["command_line"].push_back(std::string(argv[i]));
  }
  path p(fixed);
  std::string fn = p.filename().string();

  replaceAll(fn, ".nii.gz", ".nii"); // pretend we have gunziped here already
  size_t lastdot = fn.find_last_of(".");
  std::string output_filename_fixed;
  if (lastdot == std::string::npos)
    output_filename_fixed = fn + "_label_fixed.nii";
  else
    output_filename_fixed = fn.substr(0, lastdot) + "_label_fixed.nii";
  std::string output_filename_moving;
  if (lastdot == std::string::npos)
    output_filename_moving = fn + "_label_moving.nii";
  else
    output_filename_moving = fn.substr(0, lastdot) + "_label_moving.nii";

  resultJSON["output_labels_fixed"] = outdir + "/" + output_filename_fixed;
  resultJSON["output_labels_moving"] = outdir + "/" + output_filename_moving;

  constexpr unsigned int ImageDimension = 3;
  using PixelType = float;
  typedef itk::Image<float, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ImageReaderType;
  ImageReaderType::Pointer imageReaderFixed = ImageReaderType::New();
  imageReaderFixed->SetFileName(fixed);
  imageReaderFixed->Update();

  ImageReaderType::Pointer imageReaderMoving = ImageReaderType::New();
  imageReaderMoving->SetFileName(moving);
  imageReaderMoving->Update();

  using OutputPixelType = unsigned short;
  using OutputImageType = itk::Image<OutputPixelType, ImageDimension>;
  using OutputMaskType = itk::Image<unsigned char, ImageDimension>;
  using FilterType = itk::BinaryThresholdImageFilter<ImageType, OutputImageType>;
  FilterType::Pointer filterFixed = FilterType::New();
  filterFixed->SetInput(imageReaderFixed->GetOutput());
  filterFixed->SetOutsideValue(0);
  filterFixed->SetInsideValue(1);
  filterFixed->SetLowerThreshold(threshold);
  filterFixed->SetUpperThreshold(255); // max value in volume
  filterFixed->Update();

  FilterType::Pointer filterMoving = FilterType::New();
  filterMoving->SetInput(imageReaderMoving->GetOutput());
  filterMoving->SetOutsideValue(0);
  filterMoving->SetInsideValue(1);
  filterMoving->SetLowerThreshold(threshold);
  filterMoving->SetUpperThreshold(255); // max value in volume
  filterMoving->Update();

  OutputImageType::Pointer f = filterFixed->GetOutput();
  f->SetOrigin(imageReaderFixed->GetOutput()->GetOrigin());
  f->SetSpacing(imageReaderFixed->GetOutput()->GetSpacing());
  f->SetDirection(imageReaderFixed->GetOutput()->GetDirection());

  typedef itk::ConnectedComponentImageFilter<OutputImageType, OutputImageType> ConnectedComponentImageFilterType;
  ConnectedComponentImageFilterType::Pointer connectedFixed = ConnectedComponentImageFilterType::New();
  connectedFixed->SetBackgroundValue(0);
  connectedFixed->SetInput(f);
  connectedFixed->Update();

  OutputImageType::Pointer conFixed = connectedFixed->GetOutput();
  conFixed->SetOrigin(imageReaderFixed->GetOutput()->GetOrigin());
  conFixed->SetSpacing(imageReaderFixed->GetOutput()->GetSpacing());
  conFixed->SetDirection(imageReaderFixed->GetOutput()->GetDirection());

  f = filterMoving->GetOutput();
  f->SetOrigin(imageReaderMoving->GetOutput()->GetOrigin());
  f->SetSpacing(imageReaderMoving->GetOutput()->GetSpacing());
  f->SetDirection(imageReaderMoving->GetOutput()->GetDirection());

  ConnectedComponentImageFilterType::Pointer connectedMoving = ConnectedComponentImageFilterType::New();
  connectedMoving->SetBackgroundValue(0);
  connectedMoving->SetInput(f);
  connectedMoving->Update();
  OutputImageType::Pointer conMoving = connectedMoving->GetOutput();
  conMoving->SetOrigin(imageReaderMoving->GetOutput()->GetOrigin());
  conMoving->SetSpacing(imageReaderMoving->GetOutput()->GetSpacing());
  conMoving->SetDirection(imageReaderMoving->GetOutput()->GetDirection());

  // create a data dictionary for all values
  std::map<std::string, std::string> data_dictionary;

  // create outdir if it does not exist already
  path p_outdir(outdir);
  if (!is_directory(p_outdir)) {
    fprintf(stdout, "Create output directory...");
    create_directories(p_outdir);
  }

  if (1) { // save the connected components image as a single volume
    typedef itk::ImageFileWriter<OutputImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    // check if that directory exists, create before writing
    writer->SetFileName(resultJSON["output_labels_fixed"]);
    writer->SetInput(conFixed);

    std::cout << "Writing all detected lesions as a single file " << std::endl;
    std::cout << resultJSON["output_labels_fixed"] << std::endl << std::endl;
    resultJSON["output_all_lesions_fixed"] = resultJSON["output_labels_fixed"];
    try {
      writer->Update();
    } catch (itk::ExceptionObject &ex) {
      std::cout << ex << std::endl;
      return EXIT_FAILURE;
    }
  }

  if (1) { // save the connected components image as a single volume
    typedef itk::ImageFileWriter<OutputImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    // check if that directory exists, create before writing
    writer->SetFileName(resultJSON["output_labels_moving"]);
    writer->SetInput(conMoving);

    std::cout << "Writing all detected lesions as a single file " << std::endl;
    std::cout << resultJSON["output_labels_moving"] << std::endl << std::endl;
    resultJSON["output_all_lesions_moving"] = resultJSON["output_labels_moving"];
    try {
      writer->Update();
    } catch (itk::ExceptionObject &ex) {
      std::cout << ex << std::endl;
      return EXIT_FAILURE;
    }
  }

  // using LabelType = unsigned short;
  using ShapeLabelObjectType = itk::ShapeLabelObject<OutputPixelType, ImageDimension>;
  using LabelMapType = itk::LabelMap<ShapeLabelObjectType>;
  using LabelType = itk::LabelImageToShapeLabelMapFilter<OutputImageType, LabelMapType>;
  LabelType::Pointer labelFixed = LabelType::New();
  labelFixed->SetInput(connectedFixed->GetOutput());
  labelFixed->SetComputePerimeter(true);
  labelFixed->SetComputeFeretDiameter(true);
  labelFixed->Update();

  LabelType::Pointer labelMoving = LabelType::New();
  labelMoving->SetInput(connectedMoving->GetOutput());
  labelMoving->SetComputePerimeter(true);
  labelMoving->SetComputeFeretDiameter(true);
  labelMoving->Update();

  LabelMapType *labelMapFixed = labelFixed->GetOutput();
  if (labelMapFixed->GetNumberOfLabelObjects() == 0) {
    // error case
    fprintf(stderr, "Error: Could not find any lesions using the current set of thresholds\n");
  }
  labelMapFixed->SetOrigin(imageReaderFixed->GetOutput()->GetOrigin());
  labelMapFixed->SetSpacing(imageReaderFixed->GetOutput()->GetSpacing());
  labelMapFixed->SetDirection(imageReaderFixed->GetOutput()->GetDirection());

  LabelMapType *labelMapMoving = labelMoving->GetOutput();
  if (labelMapMoving->GetNumberOfLabelObjects() == 0) {
    // error case
    fprintf(stderr, "Error: Could not find any lesions using the current set of thresholds\n");
  }
  labelMapMoving->SetOrigin(imageReaderMoving->GetOutput()->GetOrigin());
  labelMapMoving->SetSpacing(imageReaderMoving->GetOutput()->GetSpacing());
  labelMapMoving->SetDirection(imageReaderMoving->GetOutput()->GetDirection());

  resultJSON["voxel_size_fixed"] = json::array();
  resultJSON["voxel_size_fixed"].push_back(imageReaderFixed->GetOutput()->GetSpacing()[0]);
  resultJSON["voxel_size_fixed"].push_back(imageReaderFixed->GetOutput()->GetSpacing()[1]);
  resultJSON["voxel_size_fixed"].push_back(imageReaderFixed->GetOutput()->GetSpacing()[2]);

  resultJSON["voxel_size_moving"] = json::array();
  resultJSON["voxel_size_moving"].push_back(imageReaderMoving->GetOutput()->GetSpacing()[0]);
  resultJSON["voxel_size_moving"].push_back(imageReaderMoving->GetOutput()->GetSpacing()[1]);
  resultJSON["voxel_size_moving"].push_back(imageReaderMoving->GetOutput()->GetSpacing()[2]);

  resultJSON["lesions_fixed"] = json::array();
  resultJSON["lesions_moving"] = json::array();

  // store the points in a container
  using PointSetType = itk::PointSet<float, 3>;
  PointSetType::Pointer fixedPointSet = PointSetType::New();
  PointSetType::Pointer movingPointSet = PointSetType::New();
  using PointType = PointSetType::PointType;
  using PointsContainer = PointSetType::PointsContainer;
  PointsContainer::Pointer fixedPointContainer = PointsContainer::New();
  PointsContainer::Pointer movingPointContainer = PointsContainer::New();

  int counter = 0;
  size_t totalVolume = 0;
  for (unsigned int n = 0; n < labelMapFixed->GetNumberOfLabelObjects(); ++n) {
    ShapeLabelObjectType *labelObject = labelMapFixed->GetNthLabelObject(n);
    if (labelObject->GetNumberOfPixels() < minPixel)
      continue; // ignore this region
    // labelObject->GetNumberOfPixels()
    json lesion;
    lesion["id"] = counter;
    std::string prefix = "lesion_";
    data_dictionary.insert({prefix + "id", "Index of the lesion object. Starts counting with 0."});
    lesion["connected_component_value"] = labelObject->GetLabel();
    data_dictionary.insert({prefix + "connected_component_value", "Each detected connected component has an index. Counting starts with 1."});
    lesion["num_voxel"] = labelObject->GetNumberOfPixels();
    data_dictionary.insert({prefix + "num_voxel", "Number of voxel that are part of this lesion."});
    lesion["physical_size"] = labelObject->GetPhysicalSize();
    data_dictionary.insert(
        {prefix + "physical_size",
         "PhysicalSize is the size of the object in physical units. It is equal to the NumberOfPixels multiplied by the physical pixel size."});
    lesion["flatness"] = labelObject->GetFlatness();
    data_dictionary.insert({prefix + "flatness", "A volumetric shape attribute computed by itk."});
    lesion["equivalent_radius"] = labelObject->GetEquivalentSphericalRadius();
    data_dictionary.insert(
        {prefix + "equivalent_radius",
         "EquivalentRadius is the equivalent radius of the hypersphere of the same size than the label object. The value depends on the image spacing."});
    lesion["roundness"] = labelObject->GetRoundness();
    data_dictionary.insert({prefix + "roundness", "A volumetric shape attribute computed by itk."});
    lesion["perimeter"] = labelObject->GetPerimeter();
    data_dictionary.insert({prefix + "perimeter", "A volumetric shape attribute computed by itk."});
    lesion["elongation"] = labelObject->GetElongation();
    data_dictionary.insert(
        {prefix + "elongation",
         "Elongation is the  ratio of the largest principal moment to the second largest principal moment. Its value is greater or equal to 1."});
    lesion["number_pixel_on_border"] = labelObject->GetNumberOfPixelsOnBorder();
    data_dictionary.insert({prefix + "number_pixel_on_border",
                            "NumberOfPixelsOnBorder is the number of pixels in the objects which are on the border of the image. A pixel on several borders (a "
                            "pixel in a corner) is counted only one time, so the size on border can't be greater than the size of the object. This attribute "
                            "is particularly useful to remove the objects which are touching too much the border."});
    lesion["feret_diameter"] = labelObject->GetFeretDiameter();
    data_dictionary.insert({prefix + "feret_diameter", "FeretDiameter is the diameter in physical units of the sphere which includes all the objects voxel."});
    lesion["perimeter_on_border_ratio"] = labelObject->GetPerimeterOnBorderRatio();
    data_dictionary.insert({prefix + "perimeter_on_border_ratio", "A volumetric shape attribute computed by itk."});
    lesion["centroid"] = json::array();
    lesion["centroid"].push_back(labelObject->GetCentroid()[0]);
    lesion["centroid"].push_back(labelObject->GetCentroid()[1]);
    lesion["centroid"].push_back(labelObject->GetCentroid()[2]);
    data_dictionary.insert({prefix + "centroid", "Location of center of mass for this lesion in bounding box coordinates. It is not constrained to be in the "
                                                 "object, and thus can be outside if the object is not convex."});
    PointType fPoint;
    fPoint[0] = labelObject->GetCentroid()[0];
    fPoint[1] = labelObject->GetCentroid()[1];
    fPoint[2] = labelObject->GetCentroid()[2];
    fixedPointContainer->InsertElement(n, fPoint);

    lesion["principal_moments"] = json::array();
    lesion["principal_moments"].push_back(labelObject->GetPrincipalMoments()[2]);
    lesion["principal_moments"].push_back(labelObject->GetPrincipalMoments()[1]);
    lesion["principal_moments"].push_back(labelObject->GetPrincipalMoments()[0]);
    data_dictionary.insert({prefix + "principal_moments", "The sorted principal moments of the shape."});
    lesion["equivalent_spherical_radius"] = labelObject->GetEquivalentSphericalRadius();
    data_dictionary.insert(
        {prefix + "equivalent_spherical_radius",
         "EquivalentRadius is the equivalent radius of the hypersphere of the same size than the label object. The value depends on the image spacing."});
    lesion["equivalent_spherical_perimeter"] = labelObject->GetEquivalentSphericalPerimeter();
    data_dictionary.insert(
        {prefix + "equivalent_spherical_perimeter",
         "EquivalentPerimeter is the equivalent perimeter of the hypersphere of the same size than the label object. The value depends on the image spacing."});
    totalVolume += labelObject->GetNumberOfPixels();

    ImageType::RegionType region = connectedFixed->GetOutput()->GetLargestPossibleRegion();
    // now create a new volume just for this lesion

    OutputMaskType::Pointer mask = OutputMaskType::New();
    mask->SetRegions(region);
    mask->Allocate();
    mask->FillBuffer(0); // density for air
    mask->SetOrigin(imageReaderFixed->GetOutput()->GetOrigin());
    mask->SetSpacing(imageReaderFixed->GetOutput()->GetSpacing());
    mask->SetDirection(imageReaderFixed->GetOutput()->GetDirection());
    itk::ImageRegionIterator<OutputImageType> imageIterator(connectedFixed->GetOutput(), region);
    itk::ImageRegionIterator<OutputMaskType> maskIterator(mask, region);
    itk::ImageRegionIterator<ImageType> inputIterator(imageReaderFixed->GetOutput(), region);
    float representative_value = 0;
    while (!imageIterator.IsAtEnd() && !maskIterator.IsAtEnd()) {
      if (imageIterator.Get() == labelObject->GetLabel()) {
        maskIterator.Set(1);
        // remember a representative value from the input for a pixel inside that lesion
        // this would be the type of lesion if we use a label field as input
        representative_value = inputIterator.Get();
      }
      ++imageIterator;
      ++maskIterator;
      ++inputIterator;
    }
    lesion["representative_input_value"] = representative_value;
    data_dictionary.insert(
        {prefix + "representative_input_value",
         "Contains a value for this lesion from the input file. If the input file contains types of lesions this value will correspond to those types."});
    // and safe that volume now
    if (exportIndividualLabels) { // save the connected components image as a single volume
      typedef itk::ImageFileWriter<OutputMaskType> WriterType;
      WriterType::Pointer writer = WriterType::New();
      // check if that directory exists, create before writing
      std::string fn = resultJSON["output_labels_fixed"];
      size_t lastdot = fn.find_last_of(".");
      std::string filename("");
      char numb[1024];
      sprintf(numb, "%04d", counter);
      if (lastdot == std::string::npos)
        filename = fn + "_ID" + numb + ".nii.gz";
      else
        filename = fn.substr(0, lastdot) + "_ID" + numb + ".nii.gz";

      writer->SetFileName(filename);
      writer->SetInput(mask);

      std::cout << "Writing fixed lesion id " << counter << std::endl;
      std::cout << " to " << filename << std::endl;
      lesion["filename"] = filename;
      try {
        writer->Update();
      } catch (itk::ExceptionObject &ex) {
        std::cout << ex << std::endl;
        return EXIT_FAILURE;
      }
    }
    resultJSON["lesions_fixed"].push_back(lesion);
    counter++;
  }
  resultJSON["num_lesions_fixed"] = counter;
  resultJSON["total_lesion_size_fixed"] = totalVolume;

  counter = 0;
  totalVolume = 0;
  for (unsigned int n = 0; n < labelMapMoving->GetNumberOfLabelObjects(); ++n) {
    ShapeLabelObjectType *labelObject = labelMapMoving->GetNthLabelObject(n);
    if (labelObject->GetNumberOfPixels() < minPixel)
      continue; // ignore this region
    // labelObject->GetNumberOfPixels()
    json lesion;
    lesion["id"] = counter;
    lesion["connected_component_value"] = labelObject->GetLabel();
    lesion["num_voxel"] = labelObject->GetNumberOfPixels();
    lesion["physical_size"] = labelObject->GetPhysicalSize();
    lesion["flatness"] = labelObject->GetFlatness();
    lesion["roundness"] = labelObject->GetRoundness();
    lesion["perimeter"] = labelObject->GetPerimeter();
    lesion["elongation"] = labelObject->GetElongation();
    lesion["number_pixel_on_border"] = labelObject->GetNumberOfPixelsOnBorder();
    lesion["feret_diameter"] = labelObject->GetFeretDiameter();
    lesion["perimeter_on_border_ratio"] = labelObject->GetPerimeterOnBorderRatio();
    lesion["equivalent_radius"] = labelObject->GetEquivalentSphericalRadius();
    lesion["centroid"] = json::array();
    lesion["centroid"].push_back(labelObject->GetCentroid()[0]);
    lesion["centroid"].push_back(labelObject->GetCentroid()[1]);
    lesion["centroid"].push_back(labelObject->GetCentroid()[2]);
    PointType fPoint;
    fPoint[0] = labelObject->GetCentroid()[0];
    fPoint[1] = labelObject->GetCentroid()[1];
    fPoint[2] = labelObject->GetCentroid()[2];
    movingPointContainer->InsertElement(n, fPoint);

    lesion["principal_moments"] = json::array();
    lesion["principal_moments"].push_back(labelObject->GetPrincipalMoments()[2]);
    lesion["principal_moments"].push_back(labelObject->GetPrincipalMoments()[1]);
    lesion["principal_moments"].push_back(labelObject->GetPrincipalMoments()[0]);
    lesion["equivalent_spherical_radius"] = labelObject->GetEquivalentSphericalRadius();
    lesion["equivalent_spherical_perimeter"] = labelObject->GetEquivalentSphericalPerimeter();
    totalVolume += labelObject->GetNumberOfPixels();

    ImageType::RegionType region = connectedMoving->GetOutput()->GetLargestPossibleRegion();
    // now create a new volume just for this lesion

    OutputMaskType::Pointer mask = OutputMaskType::New();
    mask->SetRegions(region);
    mask->Allocate();
    mask->FillBuffer(0); // density for air
    mask->SetOrigin(imageReaderMoving->GetOutput()->GetOrigin());
    mask->SetSpacing(imageReaderMoving->GetOutput()->GetSpacing());
    mask->SetDirection(imageReaderMoving->GetOutput()->GetDirection());
    itk::ImageRegionIterator<OutputImageType> imageIterator(connectedMoving->GetOutput(), region);
    itk::ImageRegionIterator<OutputMaskType> maskIterator(mask, region);
    itk::ImageRegionIterator<ImageType> inputIterator(imageReaderMoving->GetOutput(), region);
    float representative_value = 0;
    while (!imageIterator.IsAtEnd() && !maskIterator.IsAtEnd()) {
      if (imageIterator.Get() == labelObject->GetLabel()) {
        maskIterator.Set(1);
        representative_value = inputIterator.Get();
      }
      ++imageIterator;
      ++maskIterator;
      ++inputIterator;
    }
    lesion["representative_input_value"] = representative_value;

    // and safe that volume now
    if (exportIndividualLabels) { // save the connected components image as a single volume
      typedef itk::ImageFileWriter<OutputMaskType> WriterType;
      WriterType::Pointer writer = WriterType::New();
      // check if that directory exists, create before writing
      std::string fn = resultJSON["output_labels_moving"];
      size_t lastdot = fn.find_last_of(".");
      std::string filename("");
      char numb[1024];
      sprintf(numb, "%04d", counter);
      if (lastdot == std::string::npos)
        filename = fn + "_ID" + numb + ".nii.gz";
      else
        filename = fn.substr(0, lastdot) + "_ID" + numb + ".nii.gz";

      writer->SetFileName(filename);
      writer->SetInput(mask);

      std::cout << "Writing moving lesion id " << counter << std::endl;
      std::cout << " to " << filename << std::endl;
      lesion["filename"] = filename;
      try {
        writer->Update();
      } catch (itk::ExceptionObject &ex) {
        std::cout << ex << std::endl;
        return EXIT_FAILURE;
      }
    }
    resultJSON["lesions_moving"].push_back(lesion);
    counter++;
  }
  resultJSON["num_lesions_moving"] = counter;
  resultJSON["total_lesion_size_moving"] = totalVolume;

  // we have now 2 point clouds with many features that we want to
  // a) align
  // b) copy labels from one to the other to mark related objects (with scaling, removal, creation, merging)
  fixedPointSet->SetPoints(fixedPointContainer);
  movingPointSet->SetPoints(movingPointContainer);

  using MetricType = itk::EuclideanDistancePointMetric<PointSetType, PointSetType>;
  MetricType::Pointer metric = MetricType::New();

  using TransformType = itk::AffineTransform<double, 3>; // itk::TranslationTransform<double, 3>;
  TransformType::Pointer transform = TransformType::New();

  using OptimizerType = itk::LevenbergMarquardtOptimizer;
  OptimizerType::Pointer optimizer = OptimizerType::New();
  optimizer->SetUseCostFunctionGradient(false);
  using RegistrationType = itk::PointSetToPointSetRegistrationMethod<PointSetType, PointSetType>;
  RegistrationType::Pointer registration = RegistrationType::New();
  // Scale the translation components of the Transform in the Optimizer
  OptimizerType::ScalesType scales(transform->GetNumberOfParameters());
  constexpr double translationScale = 1000.0; // dynamic range of translations
  constexpr double rotationScale = 1.0;       // dynamic range of rotations
  scales.Fill(0.01);
  scales[0] = 1.0 / rotationScale;
  scales[1] = 1.0 / rotationScale;
  scales[2] = 1.0 / rotationScale;
  scales[3] = 1.0 / translationScale;
  scales[4] = 1.0 / translationScale;
  scales[5] = 1.0 / translationScale;
  unsigned long numberOfIterations = 2000;
  double gradientTolerance = 1e-4; // convergence criterion
  double valueTolerance = 1e-4;    // convergence criterion
  double epsilonFunction = 1e-5;   // convergence criterion
  optimizer->SetScales(scales);
  optimizer->SetNumberOfIterations(numberOfIterations);
  optimizer->SetValueTolerance(valueTolerance);
  optimizer->SetGradientTolerance(gradientTolerance);
  optimizer->SetEpsilonFunction(epsilonFunction);
  // Start from an Identity transform (in a normal case, the user
  // can probably provide a better guess than the identity...
  transform->SetIdentity();
  registration->SetInitialTransformParameters(transform->GetParameters());

  registration->SetMetric(metric);
  registration->SetOptimizer(optimizer);
  registration->SetTransform(transform);
  fprintf(stdout, "we have %lu in fixed, %lu points moving\n", fixedPointSet->GetNumberOfPoints(), movingPointSet->GetNumberOfPoints());
  registration->SetFixedPointSet(fixedPointSet);
  registration->SetMovingPointSet(movingPointSet);
  try {
    registration->Update();
  } catch (itk::ExceptionObject &e) {
    std::cout << e << std::endl;
    return EXIT_FAILURE;
  }
  // std::cout << "Solution = " << registration->GetTransform()->GetParameters() << std::endl;
  std::cout << "PointCloud transform = " << transform->GetParameters() << std::endl;

  // ok, if the two point clouds are sufficiently close now we can try to match closest pairs
  // this should be done to get a globally consistent solution - if two points are matched
  // to the same target the target should be bigger (merged points). What is a common growth rate?
  // after this assignment?

  //           pointsA1 pointsA2 ...
  // pointsB1      1        0  ...
  // pointsB2      0        1  ...
  // mapping points in moving to space of fixed using the transform
  if (1) { // debugging, safe two landmark sets for Amira after transformation
    std::string out_file = outdir + "/LandmarkSet_fixed.am";
    FILE *fp = fopen(out_file.c_str(), "w");
    if (fp != NULL) {
      fprintf(fp, "# HyperMesh 3D ASCII 1.0\n\ndefine Markers %lu\n\nParameters {\n   ContentType \"LandmarkSet\",\n   NumSets 1 }\n\n",
              fixedPointSet->GetNumberOfPoints());
      fprintf(fp, "Markers { float[3] Coordinates } @1\n\n");
      fprintf(fp, "# Data section follows\n\n@1\n");
      for (int i = 0; i < fixedPointSet->GetNumberOfPoints(); i++) {
        PointType fPoint;
        fixedPointSet->GetPoint(i, &fPoint);
        // PointType oPoint = transform->TransformPoint(fPoint);
        fprintf(fp, "%f %f %f\n", fPoint[0], fPoint[1], fPoint[2]);
      }
      fclose(fp);
    }

    out_file = outdir + "/LandmarkSet_moving.am";
    fp = fopen(out_file.c_str(), "w");
    if (fp != NULL) {
      fprintf(fp, "# HyperMesh 3D ASCII 1.0\n\ndefine Markers %lu\n\nParameters {\n   ContentType \"LandmarkSet\",\n   NumSets 2 }\n\n",
              movingPointSet->GetNumberOfPoints());
      fprintf(fp, "Markers { float[3] Coordinates } @1\n");
      fprintf(fp, "Markers { float[3] Coordinates2 } @2\n\n");
      fprintf(fp, "# Data section follows\n\n@1\n");
      for (int i = 0; i < movingPointSet->GetNumberOfPoints(); i++) {
        PointType fPoint;
        movingPointSet->GetPoint(i, &fPoint);
        fprintf(fp, "%f %f %f\n", fPoint[0], fPoint[1], fPoint[2]);
      }
      fprintf(fp, "\n@2\n");
      for (int i = 0; i < movingPointSet->GetNumberOfPoints(); i++) {
        PointType fPoint;
        movingPointSet->GetPoint(i, &fPoint);
        PointType oPoint = transform->TransformPoint(fPoint);
        fprintf(fp, "%f %f %f\n", oPoint[0], oPoint[1], oPoint[2]);
      }
      fprintf(fp, "\n");
      fclose(fp);
    }
  }

  // now find forward (fixed to moving) matching points using closest distances
  // this is using the index of the points in fixedPointSet and movingPointSet after transform
  class MatchInfo {
  public:
    std::vector<int> idxFixed;  // index
    std::vector<int> idxMoving; // index
    std::vector<double> distance;
    std::vector<int> label;
    MatchInfo(int a, int b, double c, int d) {
      idxFixed.push_back(a);
      idxMoving.push_back(b);
      distance.push_back(c);
      label.push_back(d);
    }
    MatchInfo(){};
  };
  std::vector<MatchInfo *> mappedPoints; // index is id in fixed, value is id in moving (after transform)
  std::vector<MatchInfo *> newPoints;
  std::vector<MatchInfo *> missingPoints;
  std::vector<MatchInfo *> mergedPoints; // points that merge together from fixed (t=0) into moving (t=1)
  for (int i = 0; i < fixedPointSet->GetNumberOfPoints(); i++) {
    PointType fPoint;
    fixedPointSet->GetPoint(i, &fPoint);
    double closestDist = std::numeric_limits<double>::max();
    int closestIdx = -1;
    for (int j = 0; j < movingPointSet->GetNumberOfPoints(); j++) {
      PointType f2Point;
      movingPointSet->GetPoint(j, &f2Point);
      PointType oPoint = transform->TransformPoint(f2Point);
      double dist = (oPoint[0] - fPoint[0]) * (oPoint[0] - fPoint[0]) + (oPoint[1] - fPoint[1]) * (oPoint[1] - fPoint[1]) +
                    (oPoint[2] - fPoint[2]) * (oPoint[2] - fPoint[2]);
      if (dist < closestDist) {
        closestIdx = j;
        closestDist = dist;
      }
    }
    if (closestIdx != -1) {
      mappedPoints.push_back(new MatchInfo(i, closestIdx, closestDist, 1));
    }
  }
  // what is the mean distance and variance for all points? We want to use 3*sigma as an outlier criterion
  double mean = 0.0;
  double std = 0.0;
  for (int i = 0; i < mappedPoints.size(); i++) {
    mean += mappedPoints[i]->distance[0];
  }
  if (mappedPoints.size() > 0)
    mean /= mappedPoints.size();
  else {
    fprintf(stderr, "Giving up, not enough points...\n");
  }
  for (int i = 0; i < mappedPoints.size(); i++) {
    std += (mappedPoints[i]->distance[0] - mean) * (mappedPoints[i]->distance[0] - mean);
  }
  if (mappedPoints.size() > 1)
    std /= (mappedPoints.size() - 1);
  else {
    fprintf(stderr, "Giving up, not enough points for std...\n");
  }
  std = sqrt(std);
  fprintf(stdout, "mean and std are: %e %e\n", mean, std);
  // remove outliers
  std::vector<MatchInfo *>::iterator it = mappedPoints.begin();
  int count = 0;
  while (it != mappedPoints.end()) {
    if (fabs((*it)->distance[0] - mean) > 3 * std) {
      fprintf(stdout, "found an outlier at position %d (distance %f > 3*std %f)\n", count, fabs((*it)->distance[0] - mean), 3 * std);
      // move this outlier to the missing points array
      missingPoints.push_back(new MatchInfo((*it)->idxFixed[0], -1, (*it)->distance[0], (*it)->label[0]));
      it = mappedPoints.erase(it);
      count++;
      continue;
    }
    count++;
    ++it;
  }
  // any point in moving (after translation) that is not already in missed + mapped is a new point
  for (int j = 0; j < movingPointSet->GetNumberOfPoints(); j++) {
    bool found = false;
    std::vector<MatchInfo *>::iterator it = mappedPoints.begin();
    while (it != mappedPoints.end()) {
      if ((*it)->idxMoving[0] == j) {
        // fprintf(stdout, "found a point %d %d\n", j, (*it)->idxMoving);
        found = true;
      }
      ++it;
    }
    it = missingPoints.begin();
    while (it != missingPoints.end()) {
      if ((*it)->idxMoving[0] == j) {
        found = true;
      }
      ++it;
    }
    if (!found) {
      // fprintf(stdout, "found a new points: %d", j);
      newPoints.push_back(new MatchInfo(-1, j, -1, 1));
    }
  }
  // if two or more points in fixed map to the same point in moving those are merged points
  for (int i = 0; i < movingPointSet->GetNumberOfPoints(); i++) {
    // for (int i = 0; i < mappedPoints.size(); i++) {
    std::set<int> fixedPoints; // this is for the mappedPoints[i].idxMoving
    int idxMoving = i;
    for (int j = 0; j < mappedPoints.size(); j++) {
      if (mappedPoints[j]->idxMoving[0] == idxMoving)
        fixedPoints.insert(mappedPoints[j]->idxFixed[0]);
    }
    if (fixedPoints.size() > 1) {
      fprintf(stdout, "Found a list of fixed points that map to the same moving point...\n  ");
      std::set<int>::iterator it = fixedPoints.begin();
      std::vector<int> fpp;
      std::vector<double> dists;
      std::vector<int> labels;
      while (it != fixedPoints.end()) {
        fprintf(stdout, "%d, ", (*it));
        fpp.push_back((*it));
        dists.push_back(-1);
        labels.push_back(-1);
        ++it;
      }
      MatchInfo *p = new MatchInfo();
      p->idxFixed = fpp;
      p->idxMoving.push_back(idxMoving);
      p->distance = dists;
      p->label = labels;
      mergedPoints.push_back(p);
      fprintf(stdout, "-> %d\n", idxMoving);
    }
  }

  fprintf(stdout, "Summary: %lu points mapped, %lu points too far away (outliers), %lu new points in moving, %lu merged points\n", mappedPoints.size(),
          missingPoints.size(), newPoints.size(), mergedPoints.size());
  if (1) { // debug the mapping
    std::string out_file = outdir + "/LandmarkSet_matched.am";
    FILE *fp = fopen(out_file.c_str(), "w");
    if (fp != NULL) {
      fprintf(fp, "# HyperMesh 3D ASCII 1.0\n\ndefine Markers %lu\n\nParameters {\n   ContentType \"LandmarkSet\",\n   NumSets 2 }\n\n", mappedPoints.size());
      fprintf(fp, "Markers { float[3] Coordinates } @1\n");
      fprintf(fp, "Markers { float[3] Coordinates2 } @2\n\n");
      fprintf(fp, "# Data section follows\n\n@1\n");
      for (int i = 0; i < mappedPoints.size(); i++) {
        int idx = mappedPoints[i]->idxFixed[0];
        PointType fPoint;
        fixedPointSet->GetPoint(idx, &fPoint);
        fprintf(fp, "%f %f %f\n", fPoint[0], fPoint[1], fPoint[2]);
      }
      fprintf(fp, "\n@2\n");
      for (int i = 0; i < mappedPoints.size(); i++) {
        int idx = mappedPoints[i]->idxMoving[0];
        PointType fPoint;
        movingPointSet->GetPoint(idx, &fPoint);
        PointType oPoint = transform->TransformPoint(fPoint);
        fprintf(fp, "%f %f %f\n", oPoint[0], oPoint[1], oPoint[2]);
      }
      fprintf(fp, "\n");
      fclose(fp);
    }
  }

  // create the final spreadsheet from resultJSON["lesions_fixed"] and resultJSON["lesions_moving"]
  // walk through all the lesions in fixed and moving
  std::vector<std::map<std::string, std::string>> csv;
  for (int i = 0; i < resultJSON["lesions_fixed"].size(); i++) {
    std::map<std::string, std::string> row;
    row.insert(std::pair<std::string, std::string>("PatientID", PatientID));
    row.insert(std::pair<std::string, std::string>("filename", fixed));
    row.insert(std::pair<std::string, std::string>("lesion_id", std::to_string(i)));
    row.insert(std::pair<std::string, std::string>("lesion_id_source", "t0"));
    data_dictionary.insert({"PatientID", "Specified on the command line of MatchPairs and should correspond to the patient identifier."});
    data_dictionary.insert(
        {"filename", "The filename of the label file used as input. Should be the filename of the baseline timepoint (t0) or the followup timepoint (t1)."});
    std::string prefix = "lesion_";
    data_dictionary.insert({"lesion_id", "An id for each lesion. Counting starts with 0."});
    data_dictionary.insert({"lesion_id_source", "Identifies the origin of this lesion. Either t0 (first argument) or t1 (second argument)."});
    std::string type = "";
    for (const auto &it : resultJSON["lesions_fixed"][i].items()) {
      for (const auto &val : it.value().items()) {
        std::string str_val;
        std::ostringstream oss;
        oss << val.value();
        str_val = oss.str();
        if (strlen(val.key().c_str()) == 0) {
          row.insert(std::pair<std::string, std::string>(prefix + it.key(), str_val));
        } else {
          row.insert(std::pair<std::string, std::string>(prefix + it.key() + "_" + val.key(), str_val));
        }
      }
    }
    // does this lesion map to another?
    //  std::vector<MatchInfo *> newPoints;
    bool in_mapped = false;
    int idxMappedTo = -1;
    for (int j = 0; j < mappedPoints.size(); j++) {
      if (mappedPoints[j]->idxFixed[0] == i) {
        in_mapped = true;
        idxMappedTo = mappedPoints[j]->idxMoving[0];
      }
    }
    if (in_mapped) {
      row.insert(std::pair<std::string, std::string>(prefix + "mapped_to_id", std::to_string(idxMappedTo)));
      data_dictionary.insert({prefix + "mapped_to_id", "If a mapped id in the other volume is found after registration this is the id of that lesion."});
      // if they map we would like to know about the volume increase between the two (i, idxMappedTo)
      json b = resultJSON["lesions_fixed"][i]["physical_size"];
      double physical_size_1 = 0.0f;
      for (const auto &val : b.items()) {
        std::string str_val;
        std::ostringstream oss;
        oss << val.value();
        str_val = oss.str();
        physical_size_1 = atof(str_val.c_str());
      }
      double physical_size_2 = 0.0f;
      b = resultJSON["lesions_moving"][idxMappedTo]["physical_size"];
      for (const auto &val : b.items()) {
        std::string str_val;
        std::ostringstream oss;
        oss << val.value();
        str_val = oss.str();
        physical_size_2 = atof(str_val.c_str());
      }
      double change = 0.0f;
      if (physical_size_1 > 0)
        change = physical_size_2 / physical_size_1;
      row.insert(std::pair<std::string, std::string>(prefix + "relative_size_change", std::to_string(change)));
      type = "mapped";
      data_dictionary.insert({prefix + "relative_size_change",
                              "The physical size of the lesion (number of voxel) at t1 divided by the size of the lesion at t0. Values larger "
                              "than 1 indicate volume increase. Values smaller than 0 indicate that a lesion is shrinking from t0 to t1."});
    } else {
      row.insert(std::pair<std::string, std::string>(prefix + "mapped_to_id", "")); // if not its an orphan
    }
    // does this lesion merge with another
    //    std::vector<MatchInfo *> mergedPoints; // points that merge together from fixed (t=0) into moving (t=1)
    idxMappedTo = -1;
    for (int j = 0; j < mergedPoints.size(); j++) {
      for (int k = 0; k < mergedPoints[j]->idxFixed.size(); k++) {
        if (mergedPoints[j]->idxFixed[k] == i) {
          row.insert(std::pair<std::string, std::string>(prefix + "merged_to", std::to_string(mergedPoints[j]->idxMoving[0]))); // if not its an orphan
        }
      }
    }
    data_dictionary.insert({prefix + "merged_to", "If two or more lesions in t0 are mapped to the same lesion in t1 they could as merged lesions. This column "
                                                  "contains the id of the target lesion in t1."});

    bool is_missing = false;
    for (int j = 0; j < missingPoints.size(); j++) {
      if (missingPoints[j]->idxFixed[0] == i) {
        row.insert(std::pair<std::string, std::string>(prefix + "removed_point", "yes"));
        type = "missing_in_t1";
        is_missing = true;
      }
    }
    if (!is_missing) {
      row.insert(std::pair<std::string, std::string>(prefix + "removed_point", "no"));
    }
    data_dictionary.insert({prefix + "removed_point",
                            "If a lesion in t0 does not have a corresponding lesion in t1 it counts as removed. This value will be yes for such a lesion."});
    row.insert(std::pair<std::string, std::string>(prefix + "type", type));
    data_dictionary.insert({prefix + "type", "The type of a lesion can be mapped, missing_in_t1, or new."});
    row.insert(std::pair<std::string, std::string>(prefix + "num_lesions", std::to_string(resultJSON["lesions_fixed"].size())));
    data_dictionary.insert({prefix + "num_lesions", "The number of lesions in this volume. This is a constant value for all lesion entries."});
    row.insert(std::pair<std::string, std::string>(prefix + "num_mapped", std::to_string(mappedPoints.size())));
    data_dictionary.insert({prefix + "num_mapped", "The number of mapped lesions in this volume. This is a constant value for all lesion entries."});
    row.insert(std::pair<std::string, std::string>(prefix + "num_missing", std::to_string(missingPoints.size())));
    data_dictionary.insert({prefix + "num_missing", "The number of missing lesions in the t0 volume. This is a constant value for all lesion entries in t0."});
    row.insert(std::pair<std::string, std::string>(prefix + "num_merged", std::to_string(mergedPoints.size())));
    data_dictionary.insert(
        {prefix + "num_merged",
         "As lesions grow they might merge with neighbors. This entry is the total number of lesions in t0 that do not map to a unique lesion "
         "in t1. The value is the same for all lesions in t0."});
    for (const auto &it : resultJSON["total_lesion_size_fixed"]) {
      std::string str_val;
      std::ostringstream oss;
      oss << it;
      str_val = oss.str();
      row.insert(std::pair<std::string, std::string>(prefix + "num_total_voxel", str_val));
    }
    data_dictionary.insert(
        {prefix + "num_total_voxel", "The sum of all voxel in this volume that are part of a lesion. The value is repeated for all entries of the volume."});

    csv.push_back(row);
  }
  for (int i = 0; i < resultJSON["lesions_moving"].size(); i++) {
    std::map<std::string, std::string> row;
    row.insert(std::pair<std::string, std::string>("PatientID", PatientID));
    row.insert(std::pair<std::string, std::string>("filename", fixed));
    row.insert(std::pair<std::string, std::string>("lesion_id", std::to_string(i)));
    row.insert(std::pair<std::string, std::string>("lesion_id_source", "t1"));
    std::string prefix = "lesion_";
    std::string type = "mapped";
    for (const auto &it : resultJSON["lesions_moving"][i].items()) {
      for (const auto &val : it.value().items()) {
        std::string str_val;
        std::ostringstream oss;
        oss << val.value();
        str_val = oss.str();
        if (strlen(val.key().c_str()) == 0) {
          row.insert(std::pair<std::string, std::string>(prefix + it.key(), str_val));
        } else {
          row.insert(std::pair<std::string, std::string>(prefix + it.key() + "_" + val.key(), str_val));
        }
      }
    }
    //  std::vector<MatchInfo *> missingPoints;  points that don't have a matching entry in moving
    bool is_new = false;
    for (int j = 0; j < newPoints.size(); j++) {
      if (newPoints[j]->idxMoving[0] == i) {
        row.insert(std::pair<std::string, std::string>(prefix + "new_point", "yes"));
        is_new = true;
        type = "new";
      }
    }
    if (!is_new) {
      row.insert(std::pair<std::string, std::string>(prefix + "new_point", "no"));
    }
    row.insert(std::pair<std::string, std::string>(prefix + "type", type));
    row.insert(std::pair<std::string, std::string>(prefix + "num_lesions", std::to_string(resultJSON["lesions_moving"].size())));
    for (const auto &it : resultJSON["total_lesion_size_moving"]) {
      std::string str_val;
      std::ostringstream oss;
      oss << it;
      str_val = oss.str();
      row.insert(std::pair<std::string, std::string>(prefix + "num_total_voxel", str_val));
    }
    // add the mapped_to_id for this moving lesion
    bool in_mapped = false;
    int idxMappedTo = -1;
    for (int j = 0; j < mappedPoints.size(); j++) {
      if (mappedPoints[j]->idxMoving[0] == i) {
        in_mapped = true;
        idxMappedTo = mappedPoints[j]->idxFixed[0];
      }
    }
    if (in_mapped) {
      row.insert(std::pair<std::string, std::string>(prefix + "mapped_to_id", std::to_string(idxMappedTo)));
    }

    csv.push_back(row);
  }

  // print out the table
  std::set<std::string> allKeys;
  for (int i = 0; i < csv.size(); i++) {
    std::map<std::string, std::string>::iterator it = csv[i].begin();
    while (it != csv[i].end()) {
      allKeys.insert(it->first);
      ++it;
    }
  }
  std::string out_file = outdir + "/stats.csv";
  FILE *fp = fopen(out_file.c_str(), "w");
  if (fp != NULL) {
    std::vector<std::string> header;
    header.assign(allKeys.begin(), allKeys.end());
    std::sort(header.begin(), header.end());
    for (int i = 0; i < header.size(); i++) {
      fprintf(fp, "%s", header[i].c_str());
      if (i < header.size() - 1) {
        fprintf(fp, ", ");
      }
    }
    fprintf(fp, "\n");
    for (int i = 0; i < csv.size(); i++) {
      for (int j = 0; j < header.size(); j++) {
        std::string val = "";
        std::map<std::string, std::string>::iterator it = csv[i].find(header[j]);
        if (it != csv[i].end()) {
          val = it->second;
        }
        fprintf(fp, "%s", val.c_str());
        if (j < header.size() - 1)
          fprintf(fp, ",");
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }

  // write the data dictionary for stats.csv
  out_file = outdir + "/stats_data_dictionary.csv";
  fp = fopen(out_file.c_str(), "w");
  if (fp != NULL) {
    for (const auto &pair : data_dictionary) {
      std::string key;
      std::string value;
      std::ostringstream oss_key;
      oss_key << pair.first;
      key = oss_key.str();
      std::ostringstream oss_value;
      oss_value << pair.second;
      value = oss_value.str();
      fprintf(fp, "\"%s\",\"%s\"\n", key.c_str(), value.c_str());
    }
    fclose(fp);
  }

  // it would be good to collect summary information based on representative_input_values field values
  // for the relative size change (average value)
  std::map<int, std::vector<double>> repr2relativeSizeChange = groupBy(csv, "lesion_representative_input_value", "lesion_relative_size_change");
  json erg = robustAnalysis(repr2relativeSizeChange);
  resultJSON["relative_size_change_by_representative_input_value"] = erg;

  std::map<int, std::vector<double>> repr2physicalSizeChange = groupBy(csv, "lesion_representative_input_value", "lesion_physical_size");
  erg = robustAnalysis(repr2physicalSizeChange);
  resultJSON["physical_size_change_by_representative_input_value"] = erg;

  std::map<int, std::vector<double>> repr2roundnessChange = groupBy(csv, "lesion_representative_input_value", "lesion_roundness");
  erg = robustAnalysis(repr2roundnessChange);
  resultJSON["roundness_change_by_representative_input_value"] = erg;

  std::map<int, std::vector<double>> repr2flatnessChange = groupBy(csv, "lesion_representative_input_value", "lesion_flatness");
  erg = robustAnalysis(repr2flatnessChange);
  resultJSON["flatness_change_by_representative_input_value"] = erg;

  std::map<int, std::vector<double>> repr2elongationChange = groupBy(csv, "lesion_representative_input_value", "lesion_elongation");
  erg = robustAnalysis(repr2elongationChange);
  resultJSON["elongation_change_by_representative_input_value"] = erg;

  std::map<int, std::vector<double>> repr2perimeterChange = groupBy(csv, "lesion_representative_input_value", "lesion_perimeter");
  erg = robustAnalysis(repr2perimeterChange);
  resultJSON["perimeter_change_by_representative_input_value"] = erg;

  std::map<int, std::vector<double>> repr2feretDiameterChange = groupBy(csv, "lesion_representative_input_value", "lesion_feret_diameter");
  erg = robustAnalysis(repr2feretDiameterChange);
  resultJSON["feret_diameter_change_by_representative_input_value"] = erg;

  std::ostringstream o;
  std::string si(outdir + "/summary");
  si.erase(std::remove(si.begin(), si.end(), '\"'), si.end());
  lastdot = si.find_last_of(".");
  if (lastdot == std::string::npos)
    si = si + ".json";
  else
    si = si.substr(0, lastdot) + ".json";

  o << si;
  // resultJSON["z_comment"] =
  //    std::string("jq -r '.lesions | map(.filename), map(.id), map(.num_voxel), map(.flatness), map(.roundness), map(.elongation) | @csv' ") + o.str();

  std::ofstream out(o.str());
  std::string res = resultJSON.dump(4) + "\n";
  out << res;
  out.close();

  fprintf(stdout, "%s", res.c_str());

  return EXIT_SUCCESS;
}
