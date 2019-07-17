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

//#include "itkMetaDataDictionary.h"
#include "json.hpp"
#include "metaCommand.h"
#include <boost/filesystem.hpp>
#include <map>

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

json resultJSON;

int main(int argc, char *argv[]) {

  itk::MultiThreaderBase::SetGlobalMaximumNumberOfThreads(4);

  MetaCommand command;
  command.SetAuthor("Hauke Bartsch");
  command.SetDescription("Simple ConnectedComponents based on itk.");
  command.AddField("infile", "Input mask", MetaCommand::STRING, true);
  command.AddField("outdir", "Output masks directory", MetaCommand::STRING, true);

  command.SetOption("Threshold", "t", false, "Specify the threshold applied to the input to create a mask (0).");
  command.AddOptionField("Threshold", "threshold", MetaCommand::FLOAT, true);

  command.SetOption("minPixel", "m", false, "Specify the minimum number of voxel in a lesion (1).");
  command.AddOptionField("minPixel", "minPixel", MetaCommand::INT, true);

  command.SetOption("Verbose", "v", false, "Print more verbose output");

  if (!command.Parse(argc, argv)) {
    return 1;
  }

  std::string input = command.GetValueAsString("infile");
  std::string outdir = command.GetValueAsString("outdir");

  int minPixel = 1;
  if (command.GetOptionWasSet("minPixel"))
    minPixel = command.GetValueAsInt("minPixel", "minPixel");

  float threshold = 0.00001; // > 0
  if (command.GetOptionWasSet("Threshold"))
    threshold = command.GetValueAsFloat("Threshold", "threshold");

  bool verbose = false;
  if (command.GetOptionWasSet("Verbose"))
    verbose = true;

  // store information in the result json file
  resultJSON["command_line"] = json::array();
  for (int i = 0; i < argc; i++) {
    resultJSON["command_line"].push_back(std::string(argv[i]));
  }
  path p(input);
  std::string fn = p.filename().string();
  size_t lastdot = fn.find_last_of(".");
  std::string output_filename;
  if (lastdot == std::string::npos)
    output_filename = fn + "_label.nii";
  else
    output_filename = fn.substr(0, lastdot) + "_label.nii";

  resultJSON["output_labels"] = outdir + "/" + output_filename;

  constexpr unsigned int ImageDimension = 3;
  using PixelType = float;
  typedef itk::Image<float, ImageDimension> ImageType;
  typedef itk::ImageFileReader<ImageType> ImageReaderType;
  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName(input);
  imageReader->Update();

  using OutputPixelType = unsigned short;
  using OutputImageType = itk::Image<OutputPixelType, ImageDimension>;
  using OutputMaskType = itk::Image<unsigned char, ImageDimension>;
  using FilterType = itk::BinaryThresholdImageFilter<ImageType, OutputImageType>;
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput(imageReader->GetOutput());
  filter->SetOutsideValue(0);
  filter->SetInsideValue(1);
  filter->SetLowerThreshold(threshold);
  filter->SetUpperThreshold(100); // max value in volume
  filter->Update();

  typedef itk::ConnectedComponentImageFilter<OutputImageType, OutputImageType> ConnectedComponentImageFilterType;
  ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
  connected->SetBackgroundValue(0);
  connected->SetInput(filter->GetOutput());
  connected->Update();
  OutputImageType::Pointer con = connected->GetOutput();
  con->SetOrigin(imageReader->GetOutput()->GetOrigin());
  con->SetSpacing(imageReader->GetOutput()->GetSpacing());
  con->SetDirection(imageReader->GetOutput()->GetDirection());

  if (1) { // save the connected components image as a single volume
    typedef itk::ImageFileWriter<OutputImageType> WriterType;
    WriterType::Pointer writer = WriterType::New();
    // check if that directory exists, create before writing
    writer->SetFileName(resultJSON["output_labels"]);
    writer->SetInput(con);

    std::cout << "Writing all segments as a single file " << std::endl;
    std::cout << resultJSON["output_labels"] << std::endl << std::endl;
    resultJSON["output_all_lesions"] = resultJSON["output_labels"];
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
  LabelType::Pointer label = LabelType::New();
  label->SetInput(connected->GetOutput());
  label->SetComputePerimeter(true);
  label->Update();

  LabelMapType *labelMap = label->GetOutput();
  if (labelMap->GetNumberOfLabelObjects() == 0) {
    // error case
    fprintf(stderr, "Error: Could not find any lesions using the current set of thresholds\n");
  }
  resultJSON["voxel_size"] = json::array();
  resultJSON["voxel_size"].push_back(imageReader->GetOutput()->GetSpacing()[0]);
  resultJSON["voxel_size"].push_back(imageReader->GetOutput()->GetSpacing()[1]);
  resultJSON["voxel_size"].push_back(imageReader->GetOutput()->GetSpacing()[2]);

  resultJSON["lesions"] = json::array();
  int counter = 0;
  size_t totalVolume = 0;
  for (unsigned int n = 0; n < labelMap->GetNumberOfLabelObjects(); ++n) {
    ShapeLabelObjectType *labelObject = labelMap->GetNthLabelObject(n);
    if (labelObject->GetNumberOfPixels() < minPixel)
      continue; // ignore this region
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

    OutputMaskType::Pointer mask = OutputMaskType::New();
    mask->SetRegions(region);
    mask->Allocate();
    mask->FillBuffer(0); // density for air
    mask->SetOrigin(imageReader->GetOutput()->GetOrigin());
    mask->SetSpacing(imageReader->GetOutput()->GetSpacing());
    mask->SetDirection(imageReader->GetOutput()->GetDirection());
    itk::ImageRegionIterator<OutputImageType> imageIterator(connected->GetOutput(), region);
    itk::ImageRegionIterator<OutputMaskType> maskIterator(mask, region);
    while (!imageIterator.IsAtEnd() && !maskIterator.IsAtEnd()) {
      if (imageIterator.Get() == labelObject->GetLabel()) {
        maskIterator.Set(1);
      }
      ++imageIterator;
      ++maskIterator;
    }
    // and safe that volume now
    if (1) { // save the connected components image as a single volume
      typedef itk::ImageFileWriter<OutputMaskType> WriterType;
      WriterType::Pointer writer = WriterType::New();
      // check if that directory exists, create before writing
      std::string fn = resultJSON["output_labels"];
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

      std::cout << "Writing lesion id " << counter << std::endl;
      std::cout << " to " << filename << std::endl;
      lesion["filename"] = filename;
      try {
        writer->Update();
      } catch (itk::ExceptionObject &ex) {
        std::cout << ex << std::endl;
        return EXIT_FAILURE;
      }
    }
    resultJSON["lesions"].push_back(lesion);
    counter++;
  }
  resultJSON["num_lesions"] = counter;
  resultJSON["total_lesion_size"] = totalVolume;

  std::string res = resultJSON.dump(4) + "\n";
  std::ostringstream o;
  std::string si(resultJSON["output_labels"]);
  si.erase(std::remove(si.begin(), si.end(), '\"'), si.end());
  lastdot = si.find_last_of(".");
  if (lastdot == std::string::npos)
    si = si + ".json";
  else
    si = si.substr(0, lastdot) + ".json";

  o << si;
  resultJSON["z_comment"] =
      std::string("jq -r '.lesions | map(.filename), map(.id), map(.num_voxel), map(.flatness), map(.roundness), map(.elongation) | @csv' ") + o.str();

  std::ofstream out(o.str());
  out << res;
  out.close();

  fprintf(stdout, "%s", res.c_str());

  return EXIT_SUCCESS;
}
