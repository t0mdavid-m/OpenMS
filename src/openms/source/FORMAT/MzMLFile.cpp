// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Chris Bielow, Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>

#include <OpenMS/FORMAT/HANDLERS/MzMLHandler.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/CVMappingFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>
#include <OpenMS/FORMAT/VALIDATORS/MzMLValidator.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataTransformingConsumer.h>
#include <OpenMS/FORMAT/HANDLERS/IndexedMzMLDecoder.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <exception>
#include <fstream>
#include <sstream>

namespace OpenMS
{

  MzMLFile::MzMLFile() :
    XMLFile("/SCHEMAS/mzML_1_10.xsd", "1.1.0"),
    indexed_schema_location_("/SCHEMAS/mzML_idx_1_10.xsd")
  {
  }

  MzMLFile::~MzMLFile() = default;

  PeakFileOptions& MzMLFile::getOptions()
  {
    return options_;
  }

  const PeakFileOptions& MzMLFile::getOptions() const
  {
    return options_;
  }

  void MzMLFile::setOptions(const PeakFileOptions& options)
  {
    options_ = options;
  }

  bool MzMLFile::hasIndex(const String& filename)
  {
    const std::streampos NOT_FOUND {-1};
    return NOT_FOUND != IndexedMzMLDecoder().findIndexListOffset(filename);
  }

  // reimplemented in order to handle index MzML
  bool MzMLFile::isValid(const String& filename, std::ostream& os)
  {
    //determine if this is indexed mzML or not
    bool indexed = false;
    TextFile file(filename, true, 4);
    String s;
    s.concatenate(file.begin(), file.end());
    if (s.hasSubstring("<indexedmzML"))
    {
      indexed = true;
    }
    // find the corresponding schema
    String current_location;
    if (indexed)
    {
      current_location = File::find(indexed_schema_location_);
    }
    else
    {
      current_location = File::find(schema_location_);
    }

    return XMLValidator().isValid(filename, current_location, os);
  }

  bool MzMLFile::isSemanticallyValid(const String& filename, StringList& errors, StringList& warnings)
  {
    // load mapping
    CVMappings mapping;
    CVMappingFile().load(File::find("/MAPPING/ms-mapping.xml"), mapping);

    // validate
    Internal::MzMLValidator v(mapping, ControlledVocabulary::getPSIMSCV());
    bool result = v.validate(filename, errors, warnings);

    return result;
  }

  void MzMLFile::loadSize(const String& filename, Size& scount, Size& ccount)
  {
    PeakMap dummy;
    Internal::MzMLHandler handler(dummy, filename, getVersion(), *this);
    handler.setOptions(options_);
    if (options_.hasFilters())
    {
      handler.setLoadDetail(Internal::XMLHandler::LD_COUNTS_WITHOPTIONS);
    }
    else
    { // no filters where specified. Just take the 'counts' attributes from the mzML file and end parsing
      handler.setLoadDetail(Internal::XMLHandler::LD_RAWCOUNTS);
    }
    
    safeParse_(filename, &handler);
    handler.getCounts(scount, ccount);
  }

  void MzMLFile::safeParse_(const String& filename, Internal::XMLHandler* handler)
  {
    try
    {
      parse_(filename, handler);
    }
    catch (Exception::BaseException& e)
    {
      String expr;
      expr += e.getFile();
      expr += "@";
      expr += e.getLine();
      expr += "-";
      expr += e.getFunction();
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, expr, String("- due to that error of type ") + e.getName());
    }
  }

  void MzMLFile::loadBuffer(const std::string& buffer, PeakMap& map)
  {
    map.reset();

    Internal::MzMLHandler handler(map, "memory", getVersion(), *this);
    handler.setOptions(options_);
    parseBuffer_(buffer, &handler);
    map.updateRanges();
  }

  void MzMLFile::load(const String& filename, PeakMap& map)
  {
    map.reset();

    //set DocumentIdentifier
    map.setLoadedFileType(filename);
    map.setLoadedFilePath(filename);

    Internal::MzMLHandler handler(map, filename, getVersion(), *this);
    handler.setOptions(options_);
    safeParse_(filename, &handler);
    map.updateRanges();
  }

  namespace
  {
    bool loadIndexedParallel_(const std::string& file, const PeakFileOptions& options, MSExperiment& map, int threads,
                              const std::function<void(Size, MSExperiment&)>& chunk_filter)
    {
      std::ifstream in(file, std::ios::binary);
      if (! in) { return false; }
      in.seekg(0, std::ios::end);
      const std::streamoff file_size = in.tellg();
      if (file_size < 64) { return false; }
      // Only plain (uncompressed) indexed mzML: non-indexed or compressed inputs fall back silently.
      char magic[2] = {0, 0};
      in.seekg(0);
      in.read(magic, 2);
      if ((magic[0] == '\x1f' && magic[1] == '\x8b') || (magic[0] == 'B' && magic[1] == 'Z') || (magic[0] == 'P' && magic[1] == 'K')) { return false; }
      const std::streamoff tail_size = std::min<std::streamoff>(file_size, 4096);
      std::string tail((Size)tail_size, '\0');
      in.seekg(file_size - tail_size);
      in.read(tail.data(), tail_size);
      if (in.gcount() != tail_size || tail.find("<indexListOffset") == std::string::npos) { return false; }

      IndexedMzMLDecoder decoder;
      const std::streampos index_offset = decoder.findIndexListOffset(file);
      if (index_offset == std::streampos(-1)) { return false; }
      IndexedMzMLDecoder::OffsetVector spectrum_offsets, chromatogram_offsets;
      if (decoder.parseOffsets(file, index_offset, spectrum_offsets, chromatogram_offsets) != 0 || spectrum_offsets.size() < 2) { return false; }
      for (Size i = 1; i < spectrum_offsets.size(); ++i)
      {
        if (spectrum_offsets[i].second <= spectrum_offsets[i - 1].second) { return false; } // index not in document order
      }
      if (spectrum_offsets.back().second >= file_size) { return false; }
      const std::streamoff first = spectrum_offsets.front().second;
      std::string header((Size)first, '\0');
      in.seekg(0);
      in.read(header.data(), first);
      if (in.gcount() != first) { return false; }
      // The header must end inside the spectrum list opening and contain no spectrum (a stale or partial index
      // would otherwise duplicate spectra in every chunk).
      const Size list_pos = header.rfind("<spectrumList");
      if (list_pos == std::string::npos || header.find("<spectrum ", list_pos + 1) != std::string::npos) { return false; }
      const bool indexed = header.find("<indexedmzML") != std::string::npos;
      const std::string closing = std::string("\n</spectrumList>\n</run>\n</mzML>\n") + (indexed ? "</indexedmzML>\n" : "");
      // header with the spectrum list count replaced, so each chunk handler reserves only its own spectra
      const Size count_pos = header.find("count=\"", list_pos);
      const Size count_end = count_pos == std::string::npos ? std::string::npos : header.find('"', count_pos + 7);
      if (count_end == std::string::npos || count_end > header.find('>', list_pos)) { return false; }
      const auto chunk_header = [&](Size n) { return header.substr(0, count_pos + 7) + std::to_string(n) + header.substr(count_end); };

      // Byte-balanced chunks of whole <spectrum> elements: at least 4 per thread and at most ~64 MiB of text
      // each, which bounds the raw text in flight. The last chunk runs to the end of the file.
      const std::streamoff body = file_size - first;
      // at most one thread per 2 MiB of spectra: more threads only add chunk overhead (79 MB file: 32 threads 0.19 s, 192 threads 0.36 s)
      threads = (int)std::min<std::streamoff>(std::max(1, threads), body / (2 << 20) + 1);
      const Size chunk_count = std::min<Size>(spectrum_offsets.size(), std::max<Size>((Size)std::max(1, threads) * 4, (Size)(body >> 26) + 1));
      const double target = double(body) / chunk_count;
      std::vector<Size> starts {0};
      for (Size i = 1; i < spectrum_offsets.size() && starts.size() < chunk_count; ++i)
      {
        if (double(spectrum_offsets[i].second - first) >= target * starts.size()) { starts.push_back(i); }
      }
      const Size parts_n = starts.size();
      std::vector<MSExperiment> parts(parts_n);
      std::vector<std::exception_ptr> errors(parts_n);
      const auto parse = [&](Size c) {
        const Size first_spectrum = starts[c];
        const Size end_spectrum = c + 1 < parts_n ? starts[c + 1] : spectrum_offsets.size();
        const std::streamoff begin = spectrum_offsets[first_spectrum].second;
        const std::streamoff end = c + 1 < parts_n ? (std::streamoff)spectrum_offsets[end_spectrum].second : file_size;
        std::ifstream f(file, std::ios::binary);
        std::string buffer = chunk_header(end_spectrum - first_spectrum);
        const Size head = buffer.size();
        buffer.resize(head + (Size)(end - begin));
        f.seekg(begin);
        f.read(buffer.data() + head, end - begin);
        if (f.gcount() != end - begin || buffer.compare(head, 9, "<spectrum") != 0)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, file, "index offset does not point to a <spectrum> element");
        }
        if (c + 1 < parts_n) { buffer += closing; }
        MzMLFile mzml;
        mzml.getOptions() = options;
        mzml.setLogType(ProgressLogger::NONE);
        mzml.loadBuffer(buffer, parts[c]);
        if (chunk_filter) { chunk_filter(first_spectrum, parts[c]); }
      };
      {
        // Parse the header alone first: initialises Xerces and the CV singletons before the threads start.
        MSExperiment warm_up;
        MzMLFile mzml;
        mzml.getOptions() = options;
        mzml.setLogType(ProgressLogger::NONE);
        mzml.loadBuffer(chunk_header(0) + closing, warm_up);
      }
#pragma omp parallel for schedule(dynamic, 1) num_threads(std::max(1, threads))
      for (SignedSize c = 0; c < (SignedSize)parts_n; ++c)
      {
        try
        {
          parse((Size)c);
        }
        catch (...)
        {
          errors[c] = std::current_exception();
        }
      }
      for (const auto& e : errors)
      {
        if (e) { std::rethrow_exception(e); }
      }
      // Each chunk handler creates its own DataProcessing objects. The serial loader shares one object per
      // dataProcessing element, and the writer compares these pointers, so equal ones are shared again.
      std::vector<DataProcessingPtr> shared;
      const auto share = [&shared](std::vector<DataProcessingPtr>& dps) {
        for (auto& dp : dps)
        {
          auto it = std::find_if(shared.begin(), shared.end(), [&dp](const DataProcessingPtr& s) { return s == dp || *s == *dp; });
          if (it == shared.end()) { shared.push_back(dp); }
          else { dp = *it; }
        }
      };
      map = std::move(parts[0]);
      map.reserveSpaceSpectra(spectrum_offsets.size());
      for (auto& spectrum : map.getSpectra()) { share(spectrum.getDataProcessing()); }
      for (auto& chromatogram : map.getChromatograms()) { share(chromatogram.getDataProcessing()); }
      for (Size c = 1; c < parts_n; ++c)
      {
        for (auto& spectrum : parts[c].getSpectra())
        {
          share(spectrum.getDataProcessing());
          map.addSpectrum(std::move(spectrum));
        }
        for (auto& chromatogram : parts[c].getChromatograms())
        {
          share(chromatogram.getDataProcessing());
          map.addChromatogram(std::move(chromatogram));
        }
        parts[c] = MSExperiment();
      }
      map.setLoadedFileType(file);
      map.setLoadedFilePath(file);
      map.updateRanges();
      return true;
    }
  }

  bool MzMLFile::loadIndexedParallel(const String& filename, PeakMap& map, int threads, const std::function<void(Size, PeakMap&)>& chunk_filter)
  {
    try
    {
      if (loadIndexedParallel_(filename, options_, map, threads, chunk_filter)) { return true; }
    }
    catch (...)
    {
    }
    map.clear(true);
    return false;
  }

  void MzMLFile::store(const String& filename, const PeakMap& map) const
  {
    Internal::MzMLHandler handler(map, filename, getVersion(), *this);
    handler.setOptions(options_);
    save_(filename, &handler);
  }

  void MzMLFile::storeBuffer(std::string& output, const PeakMap& map) const
  {
    Internal::MzMLHandler handler(map, "dummy", getVersion(), *this);
    handler.setOptions(options_);
    {
      std::stringstream os;

      //set high precision for writing of floating point numbers
      os.precision(writtenDigits(double()));

      // write data and close stream
      handler.writeTo(os);
      output = os.str();
    }
  }

  void MzMLFile::transform(const String& filename_in, Interfaces::IMSDataConsumer* consumer, bool skip_full_count, bool skip_first_pass)
  {
    // First pass through the file -> get the meta-data and hand it to the consumer
    if (!skip_first_pass) transformFirstPass_(filename_in, consumer, skip_full_count);

    // Second pass through the data, now read the spectra!
    {
      PeakMap dummy;
      Internal::MzMLHandler handler(dummy, filename_in, getVersion(), *this);
      handler.setOptions(options_);
      handler.setMSDataConsumer(consumer);
      safeParse_(filename_in, &handler);
    }
  }

  void MzMLFile::transform(const String& filename_in, Interfaces::IMSDataConsumer* consumer, PeakMap& map, bool skip_full_count, bool skip_first_pass)
  {
    // First pass through the file -> get the meta-data and hand it to the consumer
    if (!skip_first_pass)
    {
      transformFirstPass_(filename_in, consumer, skip_full_count);
    }
    // Second pass through the data, now read the spectra!
    {
      PeakFileOptions tmp_options(options_);
      Internal::MzMLHandler handler(map, filename_in, getVersion(), *this);
      tmp_options.setAlwaysAppendData(true);
      handler.setOptions(tmp_options);
      handler.setMSDataConsumer(consumer);

      safeParse_(filename_in, &handler);
    }
  }

  void MzMLFile::transformFirstPass_(const String& filename_in, Interfaces::IMSDataConsumer* consumer, bool skip_full_count)
  {
    // Create temporary objects and counters
    PeakFileOptions tmp_options(options_);
    Size scount = 0, ccount = 0;
    PeakMap experimental_settings;
    Internal::MzMLHandler handler(experimental_settings, filename_in, getVersion(), *this);

    // set temporary options for handler
    tmp_options.setMetadataOnly( skip_full_count );
    handler.setOptions(tmp_options);
    handler.setLoadDetail(Internal::XMLHandler::LD_RAWCOUNTS);

    safeParse_(filename_in, &handler);

    // After parsing, collect information
    handler.getCounts(scount, ccount);
    consumer->setExpectedSize(scount, ccount);
    consumer->setExperimentalSettings(experimental_settings);
  }

  std::map<UInt, MzMLFile::SpecInfo> MzMLFile::getCentroidInfo(const String& filename, const Size first_n_spectra_only)
  {
    bool oldoption = options_.getFillData();
    options_.setFillData(true); // we want the data as well (to allow estimation from data if metadata is missing)
    std::map<UInt, SpecInfo> ret;
    Size first_n_spectra_only_remaining = first_n_spectra_only;
    auto f = [&ret, &first_n_spectra_only_remaining](const MSSpectrum& s)
    {
        UInt lvl = s.getMSLevel();
        switch (s.getType(true))
        {
          case (MSSpectrum::SpectrumType::CENTROID):
            ++ret[lvl].count_centroided;
            --first_n_spectra_only_remaining;
            break;
          case (MSSpectrum::SpectrumType::PROFILE):
            ++ret[lvl].count_profile;
            --first_n_spectra_only_remaining;
            break;
          case (MSSpectrum::SpectrumType::UNKNOWN):  // this can only happen for spectra with very few peaks (or completely empty spectra)
            ++ret[lvl].count_unknown;
            break;
          default:
            // make sure we did not forget a case
            throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
        }

        // stop parsing after 10 or so spectra
        if (first_n_spectra_only_remaining == 0)
        {
          throw Internal::XMLHandler::EndParsingSoftly(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
        }
    };
    MSDataTransformingConsumer c;
    c.setSpectraProcessingFunc(f);
    transform(filename, &c, true, true); // no first pass

    // restore old state
    options_.setFillData(oldoption);

    return ret;
  }

} // namespace OpenMS
