/**
 * Exons.hpp
 *
 * Exons for splicegraph
 *
 * Copyright 2020 <University of Pennsylvania>
 */
#ifndef MAJIQ_ANNOTATEDTRANSCRIPTS_HPP
#define MAJIQ_ANNOTATEDTRANSCRIPTS_HPP

#include <algorithm>
#include <functional>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Contigs.hpp"
#include "GeneRegion.hpp"
#include "Genes.hpp"
#include "Interval.hpp"
#include "Regions.hpp"
#include "checksum.hpp"
#include "ConnectionData.hpp"

namespace majiq {

using AnnotatedTranscriptData_t = std::pair<geneid_t, std::vector<ClosedInterval>>;

struct AnnotatedTranscript : detail::GeneRegion<ExonIntervalT, AnnotatedTranscriptData_t> {
 public:
  using BaseT = detail::GeneRegion<ExonIntervalT, AnnotatedTranscriptData_t>;

    inline geneid_t& transcript_id() noexcept { return data.first; }
    inline const geneid_t& transcript_id() const noexcept { return data.first; }
    inline std::vector<ClosedInterval>& exons() noexcept { return data.second; }
    inline const std::vector<ClosedInterval>& exons() const noexcept { return data.second; }

  // constructors
  AnnotatedTranscript(KnownGene _gene, ExonIntervalT _coordinates)
      : BaseT{_gene, _coordinates, std::pair("N/A", std::vector<ClosedInterval>())} {}
  AnnotatedTranscript(KnownGene _gene, ExonIntervalT _coordinates, geneid_t _transcript_id, std::vector<ClosedInterval> exons)
      : BaseT{_gene, _coordinates, std::pair(_transcript_id, exons)} {}
  AnnotatedTranscript(KnownGene _gene, ExonIntervalT _coordinates, const std::pair<geneid_t, std::vector<ClosedInterval>> &_transcript_exons)
      : BaseT{_gene, _coordinates, _transcript_exons} {}


};

struct AnnotatedTranscriptsEx {
    std::vector<size_t> exon_idx_start;

    AnnotatedTranscriptsEx(std::vector<size_t>& _exon_idx_start) : exon_idx_start{_exon_idx_start} { }
};

class AnnotatedTranscripts : public detail::Regions<AnnotatedTranscript, true, 0, true>,
    public std::enable_shared_from_this<AnnotatedTranscripts>{

    using BaseT = detail::Regions<AnnotatedTranscript, true, 0, true>;

 public:

    std::vector<size_t> exon_idx_start;
    std::vector<size_t> exon_idx_end;
    std::vector<position_t> exon_start;
    std::vector<position_t> exon_end;

  AnnotatedTranscripts(const std::shared_ptr<Genes>& genes, std::vector<AnnotatedTranscript>&& x)
      : BaseT{genes, std::move(x)},
    exon_idx_start{__exon_idx_start()}, exon_idx_end{__exon_idx_end()},
    exon_start{__exon_start()}, exon_end{__exon_end()} {
    if (parents() == nullptr) {
      throw std::invalid_argument("AnnotatedTranscripts cannot have null genes");
    }
  }
  AnnotatedTranscripts(const std::shared_ptr<Genes>& genes, std::vector<AnnotatedTranscript>&& x,
      const std::vector<size_t>& _exon_idx_start, const std::vector<size_t>& _exon_idx_end,
      const std::vector<position_t>& _exon_start, const std::vector<position_t>& _exon_end)
      : BaseT{genes, std::move(x)},
    exon_idx_start{_exon_idx_start}, exon_idx_end {_exon_idx_end},
    exon_start{_exon_start}, exon_end{_exon_end} {}


  std::vector<geneid_t> transcript_ids() const {
      std::vector<geneid_t> result{size()};
      std::transform(begin(), end(), result.begin(),
                     [](const AnnotatedTranscript& a) { return a.transcript_id(); });
      return result;
  }

  std::vector<size_t> __exon_idx_start() const {
      std::vector<size_t> result;
      result.reserve(size());
      size_t offset{0};
      for (int i = 0; i < size(); i++ ) {
          result.push_back(offset);
          offset += elements_[i].exons().size();
      }
      return result;
  }

  std::vector<size_t> __exon_idx_end() const {
      std::vector<size_t> result;
      result.reserve(size());
      size_t offset{0};
      for (int i = 0; i < size(); ++i ) {
          offset += elements_[i].exons().size();
          result.push_back(offset);
      }
      return result;
  }

  std::vector<position_t> __exon_start() const {
      std::vector<position_t> result;
      for (int i = 0; i < size(); ++i ) {
          for (const auto& exon: elements_[i].exons()) {
              result.push_back(exon.start);
          }
      }
      return result;
  }

  std::vector<position_t> __exon_end() const {
      std::vector<position_t> result;
      for (int i = 0; i < size(); ++i ) {
          for (const auto& exon: elements_[i].exons()) {
              result.push_back(exon.end);
          }
      }
      return result;
  }

};


}  // namespace majiq

#endif  // MAJIQ_EXONS_HPP
