// ITK includes
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkResampleImageFilter.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkRealTimeClock.h>

// TCLAP includes
#include <tclap/ArgException.h>
#include <tclap/CmdLine.h>
#include <tclap/SwitchArg.h>
#include <tclap/ValueArg.h>

// STD includes
#include <cstdlib>

// =========================================================================
// Arguments structure
// =========================================================================
struct Arguments {
  enum DataType { _short, _int, _float };
  enum InterpolationType{ nearest, linear, cubic};
  std::string inputFileName;
  std::string outputFileName;
  DataType dataType;
  InterpolationType interpolationType;
  bool isUnsigned;
  unsigned short int x, y, z;
};

// =========================================================================
// DoIt Lippincott function
// =========================================================================
template <class T> int DoIt(const Arguments &arguments, T)
{
  // =========================================================================
  // Datatype definitions
  // =========================================================================
  using ImageType = itk::Image<T, 3>;
  using ImageReaderType = itk::ImageFileReader<ImageType>;
  using ImageWriterType = itk::ImageFileWriter<ImageType>;
  using ResampleImageFilterType = itk::ResampleImageFilter<ImageType, ImageType>;
  using InterpolatorType = typename ResampleImageFilterType::InterpolatorType;
  using BSplineInterpolatorType = itk::BSplineInterpolateImageFunction<ImageType, double, double>;
  using LinearInterpolatorType = itk::LinearInterpolateImageFunction<ImageType, double>;
  using NearestNeighborInterpolatorType = itk::NearestNeighborInterpolateImageFunction<ImageType, double>;

  // =========================================================================
  // Datatype definitions
  // =========================================================================
  auto imageReader = ImageReaderType::New();
  imageReader->SetFileName(arguments.inputFileName);
  imageReader->Update();

  // =========================================================================
  // Resample
  // =========================================================================
  auto imageResampleFilter = ResampleImageFilterType::New();


  auto inputImage = imageReader->GetOutput();
  auto outputOrigin = inputImage->GetOrigin();
  auto outputDirection= inputImage->GetDirection();
  auto inputSizePixels = inputImage->GetLargestPossibleRegion().GetSize();
  auto inputSpacing = inputImage->GetSpacing();
  double inputSize[3] = {inputSpacing[0] * inputSizePixels[0],
                         inputSpacing[1] * inputSizePixels[1],
                         inputSpacing[2] * inputSizePixels[2]};

  typename ImageType::SpacingType outputSpacing;
  outputSpacing[0] = inputSize[0] / static_cast<float>(arguments.x);
  outputSpacing[1] = inputSize[1] / static_cast<float>(arguments.y);
  outputSpacing[2] = inputSize[2] / static_cast<float>(arguments.z);

  typename ImageType::SizeType outputSize = {arguments.x, arguments.y, arguments.z};

  typename InterpolatorType::Pointer interpolator;
  switch(arguments.interpolationType)
    {
    case Arguments::InterpolationType::linear:
      interpolator = LinearInterpolatorType::New();
      break;

    case Arguments::InterpolationType::nearest:
      interpolator = NearestNeighborInterpolatorType::New();
      break;

    case Arguments::InterpolationType::cubic:
      auto inter = BSplineInterpolatorType::New();
      inter->SetSplineOrder(3);
      interpolator = inter;
    }

  auto realTimeClock = itk::RealTimeClock::New();
  imageResampleFilter->SetInput(imageReader->GetOutput());
  imageResampleFilter->SetOutputOrigin(outputOrigin);
  imageResampleFilter->SetOutputSpacing(outputSpacing);
  imageResampleFilter->SetOutputDirection(outputDirection);
  imageResampleFilter->SetInterpolator(interpolator);
  imageResampleFilter->SetSize(outputSize);
  auto start = realTimeClock->GetRealTimeStamp();
  imageResampleFilter->Update();
  auto end = realTimeClock->GetRealTimeStamp();

  std::cout << "Time (ms):" << (end - start).GetTimeInMilliSeconds() << std::endl;

  // =========================================================================
  // Write out the resampled image
  // =========================================================================
  auto imageWriter = ImageWriterType::New();
  imageWriter->SetFileName(arguments.outputFileName);
  imageWriter->SetInput(imageResampleFilter->GetOutput());
  imageWriter->Write();

  return EXIT_SUCCESS;
}

// =========================================================================
// Entry Point
// =========================================================================
int main(int argc, char **argv)
{
  // =========================================================================
  // Parse arguments
  // =========================================================================

  Arguments arguments;

  try
    {

    TCLAP::CmdLine cmd("itkImageResample");

    TCLAP::ValueArg<std::string> inputInput("i", "input", "Input Image", true, "None", "string");
    TCLAP::ValueArg<std::string> outputInput("o", "output", "Output Image", true, "None", "string");
    TCLAP::ValueArg<unsigned short int> datatypeInput("d", "datatype", "Datatype: (0) short, (1) int or (2) float", true, 0, "unsigned short int");
    TCLAP::ValueArg<unsigned short int> interpolationInput("l", "interpolation", "(0) nearest neighbours , (1) linear (2) cubic bspline", true, 0, "unsigned short int");
    TCLAP::ValueArg<unsigned short int> xInput("x", "size_x", "New size (voxels) x-axis", true, 0, "unsigned short int");
    TCLAP::ValueArg<unsigned short int> yInput("y", "size_y", "New size (voxels) y-axis", true, 0, "unsigned short int");
    TCLAP::ValueArg<unsigned short int> zInput("z", "size_z", "New size (voxels) z-axis", true, 0, "unsigned short int");
    TCLAP::SwitchArg unsignedInput("u", "unsigned", "Unsigned values", false);

    cmd.add(inputInput);
    cmd.add(outputInput);
    cmd.add(datatypeInput);
    cmd.add(interpolationInput);
    cmd.add(unsignedInput);
    cmd.add(xInput);
    cmd.add(yInput);
    cmd.add(zInput);

    cmd.parse(argc, argv);

    arguments.inputFileName = inputInput.getValue();
    arguments.outputFileName = outputInput.getValue();
    arguments.interpolationType = static_cast<Arguments::InterpolationType>(interpolationInput.getValue());
    arguments.x = xInput.getValue();
    arguments.y = yInput.getValue();
    arguments.z = zInput.getValue();
    arguments.dataType = static_cast<Arguments::DataType>(datatypeInput.getValue());
    arguments.isUnsigned = unsignedInput.getValue();

    if (arguments.dataType == Arguments::DataType::_float &&
        arguments.isUnsigned) {
      std::cerr << "Unsigned flag cannot be selected when the datatype is float"
                << std::endl;
      return EXIT_FAILURE;
    }

  } catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
  }

  // =========================================================================
  // Call the right DoIt depending on the input arguments
  // =========================================================================
  switch (arguments.dataType)
    {
    case Arguments::DataType::_short:
      return DoIt(arguments, arguments.isUnsigned ? static_cast<unsigned short int>(0)
                                                  : static_cast<short int>(0));

    case Arguments::DataType::_int:
      return DoIt(arguments, arguments.isUnsigned ? static_cast<unsigned int>(0)
                                                  : static_cast<int>(0));

    case Arguments::DataType::_float:
      return DoIt(arguments, static_cast<float>(0));
    }

  return EXIT_SUCCESS;
}
